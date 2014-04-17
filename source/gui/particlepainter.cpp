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

#include <ctime>
#include "particlepainter.h"
#include <sstream>
#include <iomanip>
#include <QtOpenGL>
#include "vortexpart.h"
#include "vortexfilament.h"
#include "turbulencepart.h"

using namespace std;

namespace Manta {

ParticlePainter::ParticlePainter(GridPainter<int>* gridRef, QWidget* par) 
    : LockedObjPainter(par), mMode(PaintVel), mLocal(0), mGridRef(gridRef),
	mLastPdata(-1), mHavePdata(false), mMaxVal(0.)
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

Real ParticlePainter::getScale() {
    if (!mObject) return 0;
    
    if (mValScale.find(mObject) == mValScale.end()) {
        Real s = 1.0;
        //if (mLocalGrid->getType() & GridBase::TypeVec3) s = 0.4;
        mValScale[mObject] = s;
    }
    return mValScale[mObject];
    
}

void ParticlePainter::processKeyEvent(PainterEvent e, int param) {
    if (e == EventNextSystem)
        nextObject();
    else if (e == EventScalePdataDown && mObject)
        mValScale[mObject] = getScale() * 0.5;
    else if (e == EventScalePdataUp && mObject)
        mValScale[mObject] = getScale() * 2.0;
    else if (e == EventToggleParticles) {
        mMode++;  // apply modulo later depending on particle system
		//if(mMode>PaintVel) mMode=PaintOff;
    }
    else return;
        
    updateText();
}

void ParticlePainter::updateText() {
    ostringstream s;
    
    if (mObject && !(mMode==PaintOff) ) {
        s << mLocal->infoString() << endl;
		s << mPdataInfo;
		if(mHavePdata) {
        	s << "-> Max " << fixed << setprecision(2) << mMaxVal << "  Scale " << getScale() << endl;
		}
    }
    mInfo->setText( s.str().c_str() );    
}


static inline void glVertex(const Vec3& v, Real dx) {
    glVertex3f(v.x * dx, v.y * dx, v.z * dx);
}

static inline void glColor(const Vec3& color) {
    glColor3f(std::max(0.0f,std::min(1.0f,color.x)), std::max(0.0f,std::min(1.0f,color.y)), std::max(0.0f,std::min(1.0f,color.z)));
}

void ParticlePainter::paint() {
    if (!mObject) return;
	if (mMode == PaintOff) return;
    float dx = mLocal->getParent()->getDx();
	mHavePdata = false;
	mMaxVal = 0.;
    
    Real scale = getScale(); 
    
    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST); // disable depth test for particles, clashes with display plane for regular ones
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
    } else if(mLocal->getType() == ParticleBase::TURBULENCE) {
        TurbulenceParticleSystem* vp = (TurbulenceParticleSystem*) mLocal;
        glPointSize(2.5);
        glColor3f(0,1,0);
        glBegin(GL_POINTS);
        for(int i=0; i<(int)vp->size(); i++) {
            Vec3 pos = (*vp)[i].pos;
            glColor((*vp)[i].color);
            glVertex(pos, dx);
            
        }   
        glEnd();
        
    } else if(mLocal->getType() == ParticleBase::PARTICLE) {
        BasicParticleSystem* bp = (BasicParticleSystem*) mLocal;

		// draw other particle data, if available
		int pdataId = mMode % (bp->getNumPdata() + 2);
		std::ostringstream infoStr;
		bool drewPoints = false;

		if( pdataId==0 ) {
			// dont draw any points
			infoStr << "Off\n";
			drewPoints = true;
		} else if( pdataId==1 ) {
			// dont draw data, only flags with center below
			infoStr << "Drawing center & flags\n";
		} else if (bp->getNumPdata() > 0)  {
			int pdNum = pdataId-2; // start at 0
			ParticleDataBase* pdb = bp->getPdata(pdNum);

			switch (pdb->getType() ) {

			case ParticleDataBase::DATA_REAL: {
				ParticleDataImpl<Real>* pdi = dynamic_cast<ParticleDataImpl<Real>*>(pdb);
				if(!pdi) break;
				mHavePdata = true;
				drewPoints = true;
				glPointSize(1.5);
				glBegin(GL_POINTS); 
				for(int i=0; i<(int)bp->size(); i++) {
            		if (!bp->isActive(i)) continue;
					Vec3 pos = (*bp)[i].pos; 
					if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
					mMaxVal = std::max( pdi->get(i), mMaxVal );
					Real val = pdi->get(i) * scale;
					glColor3f(0,val,0);
					glVertex(pos, dx); 
				}   
				glEnd();
				infoStr << "Pdata '"<<pdi->getName()<<"' #"<<pdNum<<", real\n";
				} break;

			case ParticleDataBase::DATA_INT: {
				ParticleDataImpl<int>* pdi = dynamic_cast<ParticleDataImpl<int>*>(pdb);
				if(!pdi) break;
				mHavePdata = true;
				drewPoints = true;
				glPointSize(1.5);
				glBegin(GL_POINTS); 
				for(int i=0; i<(int)bp->size(); i++) {
            		if (!bp->isActive(i)) continue;
					Vec3 pos = (*bp)[i].pos; 
					if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
					Real val = pdi->get(i);
					mMaxVal = std::max( val, mMaxVal );
					val *= scale;
					glColor3f(0,val,0);
					glVertex(pos, dx); 
				}   
				glEnd();
				infoStr << "Pdata '"<<pdi->getName()<<"' #"<<pdNum<<", int\n";
				} break;

			case ParticleDataBase::DATA_VEC3: {
				ParticleDataImpl<Vec3>* pdi = dynamic_cast<ParticleDataImpl<Vec3>*>(pdb);
				if(!pdi) break;
				mHavePdata = true;
				glBegin(GL_LINES); 
				for(int i=0; i<(int)bp->size(); i++) {
            		if (!bp->isActive(i)) continue;
					Vec3 pos = (*bp)[i].pos; 
					if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
					mMaxVal = std::max( norm(pdi->get(i)), mMaxVal );
					Vec3 val = pdi->get(i) * scale;
					glColor3f(0.5,0.0,0);
					glVertex(pos, dx); 
					pos += val;
					glColor3f(0.5,1.0,0);
					glVertex(pos, dx); 
				}   
				glEnd();
				infoStr << "Pdata '"<<pdi->getName()<<"' #"<<pdNum<<", vec3\n";
				} break;

			default: {
					// skip...
				} break;
			}
		}

		mPdataInfo = infoStr.str(); 
		// enforce refresh upon change
		if(mLastPdata!=pdataId) {
			mLastPdata = pdataId;
			updateText();
		}

		// otherwise draw center
		if(!drewPoints) {
			glPointSize(1.5);
			glBegin(GL_POINTS);

			for(int i=0; i<(int)bp->size(); i++) {
				Vec3 pos = (*bp)[i].pos;
				if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
				
				if(!bp->isActive(i) ) {
					glColor3f(1.0, 0., 0.); // deleted, red
				} else if(bp->getStatus(i) & ParticleBase::PNEW ) {
					glColor3f(0.0, 1.0, 0.); // new, greem
				} else {
					glColor3f(0, 0.0, 1.0); // regular, blue
				}
				glVertex(pos, dx);
				
			}   
			glEnd();
		}
        
		// draw basic part sys done
    }

    glPointSize(1.0);
    glEnable(GL_DEPTH_TEST); 
}

} // namespace

