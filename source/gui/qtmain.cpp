/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * QT threads
 *
 ******************************************************************************/

#include "mainwindow.h"
#include "qtmain.h"
#include "customctrl.h"

using namespace std;

// execute python script
// from pymain.cpp
extern void runScript(vector<string>& args);

namespace Manta {
    
GuiThread* gGuiThread = NULL;    
MainThread* gMainThread = NULL;    
    
MainThread::MainThread(vector<string>& args) : mArgs(args), mFinished(false) {
}

void MainThread::run() {
    runScript(mArgs);
}

void MainThread::sendAndWait(int e) {
    mMutex.lock();
    emit sendToGui(e);
    while(!mWait.wait(&mMutex, 250))
        if (gGuiThread->getWindow()->closeRequest()) {
            mMutex.unlock();
            throw Error("User interrupt");    
        }
    mMutex.unlock();
}

void MainThread::send(int e) {
    emit sendToGui(e);    
}

void MainThread::killMe() {
    if (!mFinished) {
        wait(1000);
        if (!mFinished) {
            cout << "worker thread still running, terminate" << endl;
            terminate();
            return;
        }
    }    
    wait();
}

void MainThread::wakeUp() {
    mMutex.lock();
    mWait.wakeAll();
    mMutex.unlock();
}
    
GuiThread::GuiThread(QApplication& app) : mApp(app), mWnd() {    
}

void GuiThread::sendEvent(int e) {
    mApp.postEvent(&mWnd, new QEvent((QEvent::Type)e));
}

void GuiThread::exitApp() {
    mApp.exit(1);
}

void guiMain(int argc, char* argv[]) {
    QApplication app(argc, argv);
    
	// parse arguments
	vector<string> args;
    for (int i=1;i<argc;i++) args.push_back(argv[i]);
	
	// Show file dialog if no argument is present
	if (argc <= 1) {
		QString filename = QFileDialog::getOpenFileName(0, "Open scene file", "", "Python scene files (*.py)");
		args.push_back(filename.toAscii().data());
	}
	
    GuiThread gui(app);
    MainThread worker(args);
    
    gGuiThread = &gui;
    gMainThread = &worker;
    
    // connect thread wakeup and termination signals
    QObject::connect(&worker, SIGNAL(sendToGui(int)), &gui, SLOT(sendEvent(int)));
    QObject::connect(gui.getWindow(), SIGNAL(wakeMain()), &worker, SLOT(wakeUp()));
    QObject::connect(gui.getWindow(), SIGNAL(killMain()), &worker, SLOT(killMe()));
    QObject::connect(gui.getWindow(), SIGNAL(exitApp()), &gui, SLOT(exitApp()));
    app.setQuitOnLastWindowClosed(true);
        
	// Start main program threads
    worker.start();
    app.exec();
}

void guiWaitFinish() {
    gMainThread->setFinished();    
    gMainThread->send((int)MainWnd::EventInstantKill);
    /*
    if (gGuiThread->getWindow()->closeRequest())
        return;
    
    gMainThread->sendAndWait((int)MainWnd::EventFinalUpdate);
    gGuiThread->getWindow()->pause();
    while (gGuiThread->getWindow()->pauseRequest())
        gMainThread->threadSleep(10);    */
}

//******************************************************************************
// Python adapter class


// external callback functions 
void updateQtGui(bool full, int frame, const string& curPlugin) {    
    if (!gGuiThread->getWindow()->isVisible()) return;
    if (gGuiThread->getWindow()->closeRequest()) throw Error("User interrupt");    
    
    if (full && frame >= 0) gGuiThread->getWindow()->setFrame(frame);
    gMainThread->sendAndWait(full ? (int)MainWnd::EventFullUpdate : (int)MainWnd::EventStepUpdate);
    
    if (gGuiThread->getWindow()->pauseRequest()) {
        if (!curPlugin.empty()) {
            cout << "Step: " << curPlugin <<  endl;
        }
        gGuiThread->getWindow()->setPauseStatus(true);
        while (gGuiThread->getWindow()->pauseRequest()) {            
            gMainThread->threadSleep(10);
        }
        if (gGuiThread->getWindow()->closeRequest()) throw Error("User interrupt");
        gGuiThread->getWindow()->setPauseStatus(false);
    }
    gGuiThread->getWindow()->stepReset(full);
}

} //namespace
