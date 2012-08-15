/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Base class for objects painting into the GL widget
 *
 ******************************************************************************/

#ifndef _PAINTER_H_
#define _PAINTER_H_

#include <QWidget>
#include <QLabel>
#include <map>
#include "grid.h"

namespace Manta {
    
// forward decl.
class PbClass;
    
//! Base class for all painter
/*! Derived classes have to implement paint, doEvent */
class Painter : public QObject {
    Q_OBJECT
public:
    enum PainterEvent { EventNone = 0, UpdateRequest, UpdateFull, UpdateStep,
                        EventScaleVecUp, EventScaleVecDown, EventScaleRealUp, EventScaleRealDown, EventChangePlane, 
                        EventSetPlane, EventSetDim, EventNextInt, EventNextReal, EventNextVec, EventToggleVels, EventToggleCentered,
                        EventNextMesh, EventMeshMode, EventToggleGridDisplay, EventScaleMeshUp, EventScaleMeshDown, EventMeshColorMode,
                        EventNextSystem, EventToggleParticles, EventToggleBackgroundMesh, EventSetMax};
    
    Painter(QWidget* par = 0) : QObject(par) {}
    virtual ~Painter() {}
    
    virtual void attachWidget(QLayout* layout) {}
signals:
    void setViewport(const Vec3i& gridsize);
    
public slots:
    virtual void paint() = 0;
    virtual void doEvent(int e, int param=0) = 0;
};

//! Base clas for all painters that require access to a locked PbClass
/*! Derived classes have to implement paint, update, getID, processKeyEvent. doEvent is handled in this class */
class LockedObjPainter : public Painter {
    Q_OBJECT
public:
    LockedObjPainter(QWidget* par = 0) : Painter(par), mRequestUpdate(false), mObjIndex(-1), mObject(NULL) {}

    void doEvent(int e, int param=0); // don't overload, use processKeyEvent and update instead
    
protected:
    void nextObject();
    virtual std::string getID() = 0;
    virtual void update() = 0;
    virtual void processKeyEvent(PainterEvent e, int param) = 0;
    
    bool mRequestUpdate;
    PbClass* mObject;
    int mObjIndex;
};

//! Painter object for int,Real,Vec3 grids
template<class T>
class GridPainter : public LockedObjPainter {
public:
    GridPainter(FlagGrid** flags = NULL, QWidget* par = 0);
    ~GridPainter();
    
    void paint();
    void attachWidget(QLayout* layout);
    Grid<T>** getGridPtr() { return &mLocalGrid; }
    int getPlane() { return mPlane; }
    int getDim() { return mDim; }
    int getMax() { return mMax; }
protected:
    std::string getID();
    Real getScale();
    void update();
    void updateText();
    void processKeyEvent(PainterEvent e, int param);
    void processSpecificKeyEvent(PainterEvent e, int param);
    //void paintGridLines(bool lines, bool box);
    
    Real mMaxVal;
    int mDim, mPlane, mMax;
    Grid<T>* mLocalGrid;
    FlagGrid** mFlags;
    QLabel* mInfo;
    bool mHide, mCentered, mHideLocal;
    std::map<PbClass*, Real> mValScale;
};

}

#endif