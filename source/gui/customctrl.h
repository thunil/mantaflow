/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * GUI extension from python
 *
 ******************************************************************************/

#ifndef _CUSTOMCTRL_H__
#define _CUSTOMCTRL_H__

#include <QSlider>
#include <QLabel>
#include <QLayout>
#include "pclass.h"

namespace Manta {

// fwd decl.
class Mesh;
class GuiThread;
class MainThread;
    
//! Interface for python declared controls
PYTHON class CustomControl : public PbClass {
public:
    PYTHON CustomControl();
    
    virtual void init(QLayout* layout) {};

protected:
};

//! Slider with attached text display
class TextSlider : public QSlider {
    Q_OBJECT
public:
    TextSlider(const std::string& name, float val, float min, float max);
    void attach(QLayout* layout);
    void set(float v);
    float get();
    
public slots:
    void update(int v);
        
protected:
    float mMin, mMax, mScale;
    QLabel* mLabel;    
    QString mSName;    
};
    
//! Links a slider control
PYTHON(name=Slider)
class CustomSlider : public CustomControl  {
public:
    PYTHON CustomSlider(std::string text, float val, float min, float max);
    virtual void init(QLayout* layout);
    
    PYTHON float get();
    PYTHON void set(float v);
    
protected:
    float mMin, mMax, mVal;
    std::string mSName;
    TextSlider* mSlider;
};
    

//! GUI adapter class to call from Python
PYTHON class Gui : public PbClass {
public:
    PYTHON Gui();
    
    PYTHON void setBackgroundMesh(Mesh* m);
    PYTHON void show();
    PYTHON void update();
    PYTHON void pause();
    PYTHON PbClass* addControl(PbType t);
    PYTHON void screenshot(std::string filename);
    
protected:
    GuiThread* mGuiPtr;
    MainThread* mMainPtr;
};
    
} // namespace

#endif

