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
#ifndef _QTMAIN_H_
#define _QTMAIN_H_

#include <QThread>
#include <QtGui/QApplication>
#include <vector>
#include <QMutex>
#include <QWaitCondition>
#include "mainwindow.h"
#include "pclass.h"

namespace Manta {    

//! encapsulates GUI thread
class GuiThread : public QObject {
    Q_OBJECT        
public:
    
    GuiThread(QApplication& app);
    
    //! obtain window handle
    inline MainWnd* getWindow() { return &mWnd; }
    
public slots:
    void sendEvent(int e);
    void exitApp();
    
protected:
    QApplication& mApp;
    MainWnd mWnd;
};

//! encapsulates working/python thread
class MainThread : public QThread {
    Q_OBJECT        
public:
    MainThread(std::vector<std::string>& args);
    
    //! send event to GUI and wait for completion
    void sendAndWait(int e);
    void send(int e);
    
    //! sleep for given number of milliseconds
    inline void threadSleep(int msec) { msleep(msec); }
    inline bool isFinished() { return mFinished; }
    inline void setFinished() { mFinished = true; }
    
public slots:
    void wakeUp();
    void killMe();
    
signals:
    void sendToGui(int event);
    
protected:
    QMutex mMutex;
    QWaitCondition mWait;
    bool mFinished;
    std::vector<std::string> mArgs;
    void run();    
};
  
} // namespace

#endif