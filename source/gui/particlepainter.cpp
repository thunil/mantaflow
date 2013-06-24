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
#include "vortexpart.h"
#include "vortexfilament.h"
#include "flip.h"

using namespace std;

namespace Manta {

ParticlePainter::ParticlePainter(GridPainter<int>* gridRef, QWidget* par) 
    : LockedObjPainter(par), mHide(false), mLocal(0), mGridRef(gridRef)
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
        s << mLocal->infoString() << endl;
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
    
    Real scale = 0.4;
    
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    
    // obtain current plane
    int dim = mGridRef->getDim();
    Real factor = mGridRef->getMax() / mLocal->getParent()->getGridSize()[dim];
    int plane = factor * mGridRef->getPlane();
    
    // draw points
    if(mLocal->getType() == ParticleBase::VORTEX) {
        VortexParticleSystem* vp = (VortexParticleSystem*) mLocal;
        glColor3f(1,1,0);
        for(int i=0; i<vp->size(); i++) {
            if (vp->isActive(i)) {
                Vec3 pos = (*vp)[i].pos;
            
                glPointSize((*vp)[i].sigma);

                glBegin(GL_POINTS);
                glVertex(pos, dx);
                glEnd();
            }
        }        
    } else if (mLocal->getType() == ParticleBase::FLIP) {
        FlipSystem* fp = (FlipSystem*) mLocal;
        glColor3f(0,1,1);
        glPointSize(1.0);
        glBegin(GL_LINES);
            
        for(int i=0; i<fp->size(); i++) {
            if (fp->isActive(i)) {
                Vec3 pos = (*fp)[i].pos;
                Vec3 vel = (*fp)[i].vel;
            
                if (pos[dim] >= plane && pos[dim] <= plane + 1.0f) {
                    glVertex(pos, dx);
                    glVertex(pos + vel * scale, dx);                    
                }
            }
        }   
        glEnd();
        glPointSize(1.0);
        glBegin(GL_POINTS);
            
        for(int i=0; i<fp->size(); i++) {
            if (fp->isActive(i)) {
                Vec3 pos = (*fp)[i].pos;
                if (pos[dim] >= plane && pos[dim] <= plane + 1.0f)
                    glVertex(pos, dx);
            }
        }   
        glEnd();
    } else if (mLocal->getType() == ParticleBase::FILAMENT) {
        VortexFilamentSystem* fp = (VortexFilamentSystem*) mLocal;
        glColor3f(1,1,0);
            
        for(int i=0; i<fp->segSize(); i++) {
            if (!fp->isSegActive(i)) continue;
            const VortexRing& r = fp->seg(i);
            
            glPointSize(1.0);
            glBegin(GL_LINES);
            for(int j=0; j<r.size(); j++) {
                glVertex( (*fp)[r.idx0(j)].pos, dx);
                glVertex( (*fp)[r.idx1(j)].pos, dx);
            }
            glEnd();
            
            /*glPointSize(3.0);
            glBegin(GL_POINTS);
            glVertex((*fp)[r.idx0(0)].pos,dx);
            glEnd();        */
        }   
    } else if(mLocal->getType() == ParticleBase::PARTICLE) {
        TracerParticleSystem* vp = (TracerParticleSystem*) mLocal;
        glPointSize(0.5);
        glColor3f(0,1,0);
        glBegin(GL_POINTS);
        for(int i=0; i<(int)vp->size(); i++) {
            Vec3 pos = (*vp)[i].pos;
            
            glVertex(pos, dx);
            
        }   
        glEnd();
        
    }
    glPointSize(1.0);
    /*if(mLocal->getType() == ParticleBase::TRACER) {
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
    }*/
//#endif
}

} // namespace

