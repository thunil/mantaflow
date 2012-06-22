/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Painting particle systems
 *
 ******************************************************************************/

#include "particlepainter.h"
#include <QtOpenGL>
#include <sstream>
#include <iomanip>
#ifdef MESHCODE
    #include "vortexpart.h"
#endif
#include "test.h"

using namespace std;

namespace Manta {

ParticlePainter::ParticlePainter(QWidget* par) 
    : LockedObjPainter(par), mHide(false), mLocal(0)
{    
    mInfo = new QLabel();
}

ParticlePainter::~ParticlePainter() {
    if (mLocal)
        delete mLocal;
}

void ParticlePainter::attachWidget(QLayout* layout) {
    layout->addWidget(mInfo);
}

void ParticlePainter::update() {
    ParticleBase* src = (ParticleBase*) mObject;
    
    // always reallocate
    if (mLocal) 
        delete mLocal;
    
    mLocal = src->clone();
    
    updateText();    
}

string ParticlePainter::getID() { return "ParticleBase"; }


void ParticlePainter::processKeyEvent(PainterEvent e, int param) {
    if (e == EventNextSystem)
        nextObject();
    else if (e == EventToggleParticles) {
        mHide = !mHide;
    }
    else return;
        
    updateText();
}

void ParticlePainter::updateText() {
    stringstream s;
    
    if (mObject && !mHide) {
        s << "ParticleSystem '" << mLocal->getName() << "' [" << mLocal->size() << " parts]" << endl;
    }
    mInfo->setText(s.str().c_str());    
}


static inline void glVertex(const Vec3& v, Real dx) {
    glVertex3f(v.x * dx, v.y * dx, v.z * dx);
}
/*
static inline void glColor(const Vec3& color) {
    glColor3f(std::max(0.0f,std::min(1.0f,color.x)), std::max(0.0f,std::min(1.0f,color.y)), std::max(0.0f,std::min(1.0f,color.z)));
}*/

void ParticlePainter::paint() {
    if (!mObject || mHide) return;
    float dx = mLocal->getParent()->getDx();
    
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    
//#ifdef MESHCODE
    // draw vertex points
    if(mLocal->getType() == ParticleBase::VORTEX) {
        VortexParticleSystem* vp = (VortexParticleSystem*) mLocal;
        glColor3f(1,1,0);
        glPointSize(3.0);
        for(int i=0; i<(int)mLocal->size(); i++) {
            Vec3 pos = (*vp)[i].pos;
            
            if ((*vp)[i].vort > 0)
                glColor3f(1,0.2,0);
            else
                glColor3f(0,0.2,1);
        
            glBegin(GL_POINTS);
            glVertex(pos, dx);
            glEnd();
        }        
        glPointSize(1.0);
    }
    if(mLocal->getType() == ParticleBase::TRACER) {
        TracerParticleSystem* vp = (TracerParticleSystem*) mLocal;
        glPointSize(3.0);
        for(int i=0; i<(int)mLocal->size(); i++) {
            Vec3 pos = (*vp)[i].pos;
            
            glColor3f((*vp)[i].color.x,(*vp)[i].color.y,(*vp)[i].color.z);
            glBegin(GL_POINTS);
            glVertex(pos, dx);
            glEnd();
        }        
        glPointSize(1.0);
    }
//#endif
}

} // namespace

