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

#ifndef _MAINWINDOW_H_
#define _MAINWINDOW_H_

#include <QMainWindow>
#include "glwidget.h"
#include "customctrl.h"
#include "painter.h"
#include <vector>

namespace Manta {
class Mesh;    
	
class MainWnd : public QMainWindow
{
Q_OBJECT
public:    
	enum EventType { EventFullUpdate = QEvent::User, EventGuiShow, 
		EventStepUpdate, EventFinalUpdate, EventInstantKill, EventSet2DCam };
	
	MainWnd();
	virtual ~MainWnd();
	bool event(QEvent* e);
	void keyPressEvent(QKeyEvent* e);
	void keyReleaseEvent(QKeyEvent* e);
	inline bool pauseRequest() { return mRequestPause && !mRequestClose;  }
	inline bool closeRequest() { return mRequestClose;  }
	void setPauseStatus(bool v);
	void stepReset(bool fullUpdate) { if (mStep == 1 || (mStep == 2 && fullUpdate)) {mRequestPause = true; mStep = 0;} }
	void requestClose() { mRequestClose =true; }
	void setFrame(int f);
	void setBackground(Mesh *m) { emit setBackgroundMesh(m); }
	
public slots:
	void pause();
	void play();
	void step();
	void addControl(void* ctrl);
	void screenshot(QString file);
	void clickLine(QPoint pos, float p0, float p1,float p2, float q0, float q1, float q2);
	
signals:
	void painterEvent(int e, int param=0);    
	void wakeMain();
	void setBackgroundMesh(Mesh* bgr);
	void killMain();
	void exitApp();
	
protected:
	bool mPaused, mRequestPause, mRequestClose;
	int mStep;
	GLWidget* mGlWidget;
	QAction* mAcPlay, *mAcPause;
	std::vector<Painter*> mPainter;
	std::vector<CustomControl*> mCtrls;
	QLabel* mInfo;
	QVBoxLayout* mPainterLayout;

    QGraphicsScene *scene; // NT_DEBUG
    QGraphicsView *view;
    QGraphicsPixmapItem *item; // NT_DEBUG
};

}

#endif
