/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * QT main window
 *
 ******************************************************************************/

#include "mainwindow.h"
#include "qtmain.h"

#include <QtGui/QLabel>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QAction>
#include <QtOpenGL>
#include <sstream>
#include "meshpainter.h"
#include "particlepainter.h"

using namespace std;

namespace Manta {

MainWnd::MainWnd() : mPaused(false), mRequestPause(false), mRequestClose(false), mStep(0), QMainWindow(0)
{
    // Frame info label
    mInfo = new QLabel;
    setFrame(0);
    
    // register GL widget
    mGlWidget = new GLWidget();
    setCentralWidget(mGlWidget);  
        
    // register grid painters
    mPainterLayout = new QVBoxLayout;    
    mPainterLayout->setAlignment(Qt::AlignTop);
    mPainterLayout->addWidget(mInfo);
    GridPainter<int>* intPainter = new GridPainter<int>(NULL, this);     
    mPainter.push_back(new GridPainter<Real>((FlagGrid**)intPainter->getGridPtr(), this));    
    mPainter.push_back(new GridPainter<Vec3>(NULL, this));    
    mPainter.push_back(intPainter);
    MeshPainter* ptr = new MeshPainter(this);
    mPainter.push_back(new ParticlePainter(intPainter, this));
    mPainter.push_back(ptr);    
    connect(this, SIGNAL(setBackgroundMesh(Mesh*)), ptr, SLOT(setBackgroundMesh(Mesh*)));
    for (int i=0; i<(int)mPainter.size(); i++) {
        connect(mGlWidget, SIGNAL(paintSub()), mPainter[i], SLOT(paint()));
        connect(mGlWidget, SIGNAL(painterEvent(int, int)), mPainter[i], SLOT(doEvent(int, int)));
        connect(this, SIGNAL(painterEvent(int, int)), mPainter[i], SLOT(doEvent(int, int)));
        connect(mPainter[i], SIGNAL(setViewport(const Vec3i&)), mGlWidget, SLOT(setViewport(const Vec3i&)));
        mPainter[i]->attachWidget(mPainterLayout);
    }
    
    // docking widget for painters
    QDockWidget* painterDock = new QDockWidget("Info", this);
    QWidget* painterProxy = new QWidget;
    painterProxy->setLayout(mPainterLayout);    
    painterDock->setWidget(painterProxy);
    painterDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    addDockWidget(Qt::RightDockWidgetArea, painterDock);
     
    // Top toolbar
    QToolBar* toolbar = addToolBar("Control");
    toolbar->setAllowedAreas(Qt::TopToolBarArea);
    toolbar->setMovable(false);    
    mAcPlay = toolbar->addAction(QIcon(":/play.png"),"Play");
    mAcPlay->setStatusTip("Continue simulation");
    connect(mAcPlay, SIGNAL(triggered()), SLOT(play()));
    mAcPause = toolbar->addAction(QIcon(":/pause.png"),"Pause");
    mAcPause->setStatusTip("Pause simulation");
    connect(mAcPause, SIGNAL(triggered()), SLOT(pause()));
    emit play();
    
    // build menu
    /*QAction* a = new QAction(this);
    a->setText( "Quit" );
    connect(a, SIGNAL(triggered()), SLOT(close()) );
    menuBar()->addMenu( "File" )->addAction( a );        */
    
    mGlWidget->setFocus();
    this->raise();
    this->activateWindow();
}

void MainWnd::addControl(void* ctrl) {
    CustomControl* control = (CustomControl*) ctrl;
    mCtrls.push_back(control);
    control->init(mPainterLayout);
}

void MainWnd::setFrame(int f) {
    std::stringstream s;
    s << "Simulation frame " << f;
    mInfo->setText(s.str().c_str());
}

void MainWnd::setPauseStatus(bool v)
{
    mPaused = v; 
}

bool MainWnd::event(QEvent* e) {
    if (e->type() == (QEvent::Type)EventGuiShow) {
        if (!mRequestClose) {
            this->show();
            emit painterEvent(Painter::UpdateFull);
            mGlWidget->updateGL();
        }
        emit wakeMain();
        return true;
    }
    else if (e->type() == (QEvent::Type)EventFullUpdate) {        
        if (!mRequestClose) {
            emit painterEvent(Painter::UpdateFull);
            mGlWidget->updateGL();
        }
        emit wakeMain();
        return true;
    }
    else if (e->type() == (QEvent::Type)EventStepUpdate) {        
        if (!mRequestClose) {
            if (mRequestPause)
                emit painterEvent(Painter::UpdateFull);
            else
                emit painterEvent(Painter::UpdateStep);
            mGlWidget->updateGL();
        }
        emit wakeMain();
        return true;
    }
    else if (e->type() == (QEvent::Type)EventFinalUpdate) {        
        if (!mRequestClose) {
            emit painterEvent(Painter::UpdateFull);
            mGlWidget->updateGL();
        }
        mRequestClose = true;
        emit wakeMain();
        return true;
    }
    else if (e->type() == (QEvent::Type)EventInstantKill) {        
        emit killMain();
        emit exitApp();
        return true;
    }
    
    return QMainWindow::event(e);
}

void MainWnd::keyPressEvent(QKeyEvent* e) {
    if (e->key() == Qt::Key_Escape) {
        mRequestClose = true;
        emit killMain();
        this->close();
    } else if (e->key() == Qt::Key_Space) {
        if (mRequestClose) {
            emit killMain();
            this->close();
        } else {
            emit painterEvent(mPaused ? Painter::UpdateFull : Painter::UpdateRequest);
            mGlWidget->updateGL();
        }
    } else if (e->key() == Qt::Key_P) {
        if (mRequestClose) {
            emit killMain();
            this->close();
        } else if (mRequestPause)
            emit play();
        else    
            emit pause();
    } else if (e->key() == Qt::Key_L) {
        if (mRequestClose) {
            emit killMain();
            this->close();
        } else if (mRequestPause) {
            mRequestPause = false;
            mStep = (e->modifiers() & Qt::ShiftModifier) ? 1 : 2;                
        } else
            emit pause();
    } else
        QMainWindow::keyPressEvent(e);
}

void MainWnd::pause() {
    mRequestPause = true;
    mAcPlay->setEnabled(true);
    mAcPause->setEnabled(false);    
}

void MainWnd::play() {
    mRequestPause = false;
    mAcPlay->setEnabled(false);
    mAcPause->setEnabled(true);    
}

void MainWnd::step() {
    mStep = 2;
    mRequestPause = false;
}

MainWnd::~MainWnd() {
}

void MainWnd::screenshot(QString file) {
    mGlWidget->screenshot(file);
}


}
