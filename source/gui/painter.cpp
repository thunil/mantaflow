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

#include "painter.h"
#include <QtOpenGL>
#include <sstream>
#include <iomanip>

using namespace std;

namespace Manta {
    
//******************************************************************************
// Base class
    
void LockedObjPainter::doEvent(int e, int param) {
    // try to obtain valid handle
     if (!mObject)
         nextObject();
        
    // filter update events
    if (e == UpdateFull) {
        // always update
        if (mObject) {
            mObject->lock();
            update();
            mObject->unlock();
            mRequestUpdate = false;
        }
    } else if (e == UpdateRequest) {
        // update if resource is available, otherwise wait until next step
        mRequestUpdate = true;        
        if (mObject) {
            if (mObject->tryLock()) {
                update();
                mRequestUpdate = false;
                mObject->unlock();
            }
        }
    } else if (e == UpdateStep) {
        // update if requested only
        if (mRequestUpdate) {
            if (mObject) {
                mObject->lock();
                update();
                mObject->unlock();
                mRequestUpdate = false;
            }
        }
    } else {    
        // pass on all other events
        processKeyEvent((PainterEvent)e, param);
    }
}

void LockedObjPainter::nextObject() {
    if (PbClass::getNumInstances() == 0) return;
    
    int oldIndex = mObjIndex;
    for(;;) {
        mObjIndex = (mObjIndex + 1) % PbClass::getNumInstances();
        if (oldIndex == mObjIndex) break;
        
        PbClass* obj = PbClass::getInstance(mObjIndex);
        if (obj->canConvertTo(getID()) && !obj->isHidden()) {
            mObject = obj;
            doEvent(UpdateRequest);
            return;
        }
        if (oldIndex < 0) oldIndex = 0; // prevent endless loop on first run
    } 
}

//******************************************************************************
// Grid painter

template<class T>
GridPainter<T>::GridPainter(FlagGrid** flags, QWidget* par) 
    : LockedObjPainter(par), mHide(false), mHideLocal(false), mCentered(true), mLocalGrid(NULL), mFlags(flags), mMaxVal(0), mMax(0)
{
    mDim = 2; // Z plane
    mPlane = 0;
    mInfo = new QLabel();
    
}

template<class T>
GridPainter<T>::~GridPainter() {
    if (mLocalGrid)
        delete mLocalGrid;
}

template<class T>
void GridPainter<T>::attachWidget(QLayout* layout) {
    layout->addWidget(mInfo);
}

template<class T>
void GridPainter<T>::update() {
    Grid<T>* src = (Grid<T>*) mObject;
    
    if (!mLocalGrid) {
        mLocalGrid = new Grid<T>(src->getParent());
        // int grid is base for resolution
        if (src->getType() & GridBase::TypeInt)
            emit setViewport(src->getSize());
    }
    // reallocate if dimensions changed
    if (mLocalGrid->getSize() != src->getSize()) {
        delete mLocalGrid;
        mLocalGrid = new Grid<T>(src->getParent());
        // int grid is base for resolution
        if (src->getType() & GridBase::TypeInt)
            emit setViewport(src->getSize());
    }
    
    *mLocalGrid = *src; // copy grid data and type marker
    mLocalGrid->setName(src->getName());
    mLocalGrid->setParent(src->getParent());    
    mMaxVal = mLocalGrid->getMaxAbsValue();
    
    mPlane = clamp(mPlane, 0, mLocalGrid->getSize()[mDim]-1);
    
    updateText();    
}

template<> string GridPainter<int>::getID() { return "Grid<int>"; }
template<> string GridPainter<Vec3>::getID() { return "Grid<Vec3>"; }
template<> string GridPainter<Real>::getID() { return "Grid<Real>"; }

template<class T>
void GridPainter<T>::processKeyEvent(PainterEvent e, int param)
{
    if (e == EventSetDim) {
        mDim = param;
        if (mLocalGrid->is2D()) mDim = 2;
    } else if (e == EventSetMax) {
        mMax = param;
    } else if (e == EventSetPlane) {
        mPlane = param;
        if (mObject) {
            if (mMax>0)
                mPlane = mPlane * mLocalGrid->getSize()[mDim] / mMax;
            mPlane = clamp(mPlane, 0, mLocalGrid->getSize()[mDim]-1);
        }
    } else if (e == EventToggleGridDisplay)
        mHide = !mHide;
    else
        processSpecificKeyEvent(e, param);
    
    updateText();
}

// get scale value for current grid from map, or create new
template<class T>
Real GridPainter<T>::getScale() {
    if (!mObject) return 0;
    
    if (mValScale.find(mObject) == mValScale.end()) {
        // init new scale value
        Real s = 1.0;
        if (mLocalGrid->getType() & GridBase::TypeVec3)
            s = 0.4;
        else if (mLocalGrid->getType() & GridBase::TypeLevelset)
            s = 1.0; 
        mValScale[mObject] = s;
    }
    return mValScale[mObject];
    
}

//******************************************************************************
// Grid painter class specializations

template<>
void GridPainter<int>::processSpecificKeyEvent(PainterEvent e, int param) {
    if (e == EventNextInt)
        nextObject();
}

template<>
void GridPainter<Real>::processSpecificKeyEvent(PainterEvent e, int param) {
    if (e == EventNextReal)
        nextObject();
    else if (e == EventScaleRealDown && mObject)
        mValScale[mObject] = getScale() * 0.5;
    else if (e == EventScaleRealUp && mObject)
        mValScale[mObject] = getScale() * 2.0;
}

template<>
void GridPainter<Vec3>::processSpecificKeyEvent(PainterEvent e, int param) {
    if (e == EventNextVec)
        nextObject();
    else if (e == EventScaleVecDown && mObject)
        mValScale[mObject] = getScale() * 0.5;
    else if (e == EventScaleVecUp && mObject)
        mValScale[mObject] = getScale() * 2.0;
    else if (e == EventToggleVels)
        mHideLocal = !mHideLocal;    
    else if (e== EventToggleCentered)
        mCentered = !mCentered;
}

template<> void GridPainter<int>::updateText() {
    stringstream s;
    s << "Display Plane " << mPlane << " [" << (char)('X' + mDim) << "]" << endl << endl;
    
    if (mObject) {
        s << "Solver '" << mLocalGrid->getParent()->getName() << "'" << endl;
        s << "Grid resolution [" << mLocalGrid->getSizeX() << ", " << mLocalGrid->getSizeY() << ", " << mLocalGrid->getSizeZ() << "]" << endl;
        if (!mHide)
            s << "Int Grid '" << mLocalGrid->getName() << "'" << endl;
    }    
    mInfo->setText(s.str().c_str());    
}

template<> void GridPainter<Real>::updateText() {
    stringstream s;
    if (mObject && !mHide) {
        s << "Real Grid '" << mLocalGrid->getName() << "'" << endl;
        s << "-> Max " << fixed << setprecision(2) << mMaxVal << "  Scale " << getScale() << endl;
    }
    mInfo->setText(s.str().c_str());    
}

template<> void GridPainter<Vec3>::updateText() {
    stringstream s;
    if (mObject && !mHide && !mHideLocal) {
        s << "Vec Grid '" << mLocalGrid->getName() << "'" << endl;
        s << "-> Max norm " << fixed << setprecision(2) << mMaxVal << "  Scale " << getScale() << endl;
    }
    mInfo->setText(s.str().c_str());
}

// compute line intersection with the display plane
Vec3i getQuad(const Vec3& l0, const Vec3& l1, int dim, int plane, Real dx) {
    Vec3 n(0.); n[dim] = 1;
    Vec3 p0 = n*(plane+0.5);
    Vec3 e = (l1-l0)/dx;
    Vec3 e0 = l0/dx;
    Real dotP = dot(p0-e0,n);
    Real dotE = dot(e,n);
    if (dotE == 0) 
        return Vec3i(-1,-1,-1);
    Vec3 s = e0 + (dotP/dotE)*e;
    return toVec3i(s);
}

template<> string GridPainter<int>::clickLine(const Vec3& p0, const Vec3& p1) { 
    if (!mObject) return "";
    Vec3i s = getQuad(p0,p1,mDim,mPlane,mLocalGrid->getDx());
    if (!mLocalGrid->isInBounds(s)) return "";
    stringstream m;
    m << "Grid [ " << s.x << ", " << s.y << ", " << s.z << " ]" << endl << mLocalGrid->getName() << ": " << mLocalGrid->get(s) << endl;
    return m.str();
}

template<> string GridPainter<Real>::clickLine(const Vec3& p0, const Vec3& p1) { 
    if (!mObject) return "";
    Vec3i s = getQuad(p0,p1,mDim,mPlane,mLocalGrid->getDx());
    if (!mLocalGrid->isInBounds(s)) return "";
    stringstream m;
    m << mLocalGrid->getName() << ": " << setprecision(2) << mLocalGrid->get(s) << endl;
    return m.str();
}

template<> string GridPainter<Vec3>::clickLine(const Vec3& p0, const Vec3& p1) {
    if (!mObject) return "";
    Vec3i s = getQuad(p0,p1,mDim,mPlane,mLocalGrid->getDx());
    if (!mLocalGrid->isInBounds(s)) return "";
    stringstream m;
    m << mLocalGrid->getName() << ": [ " << setprecision(2) << mLocalGrid->get(s).x << ", " <<
                                         mLocalGrid->get(s).y << ", " << mLocalGrid->get(s).z << " ]" << endl;
    return m.str();
}


//******************************************************************************
// Actual painting functions

// GL helper functions

// Macro to iterate through one plane
#define FOR_P_SLICE(__g,__dim,__plane) \
    for(Vec3i __g0(__fRange(Vec3i(0,0,0),__dim,__plane)), __g1(__fRange((__g)->getSize(),__dim,__plane+1)), p(__g0); p.z<__g1.z; p.z++) \
        for(p.y=__g0.y; p.y < __g1.y; p.y++) \
            for(p.x=__g0.x; p.x < __g1.x; p.x++)
inline Vec3i __fRange(Vec3i size, int dim, int plane) { Vec3i p(size); p[dim]=plane; return p; }
              
// coordinate system :
// cell center(i,j,k) -> (i+0.5,j+0.5,k+0.5) / N
// 

void getCellCoordinates(const Vec3i& pos, Vec3 box[4], int dim) {
    int dim2=(dim+1)%3;
    Vec3 p0(pos.x, pos.y, pos.z);
    Vec3 p1(pos.x+1, pos.y+1, pos.z+1);
    p1[dim] = p0[dim] = pos[dim] + 0.5;
    box[0] = p0;
    box[3] = p0; box[3][dim2] = p1[dim2];
    box[1] = p1; box[1][dim2] = p0[dim2];
    box[2] = p1;
}
static inline void glVertex(const Vec3& v, const float dx) {
    glVertex3f(v.x * dx, v.y * dx, v.z * dx);
}
void glBox(const Vec3& p0, const Vec3& p1, const float dx) {
    const int box[24] = {0,1,0,2,0,4,7,6,7,5,7,3,1,3,1,5,2,3,2,6,4,5,4,6};
    for (int i=0;i<24;i++) {
        const int b = box[i];
        glVertex(Vec3( (b&1) ? p1.x : p0.x, (b&2) ? p1.y : p0.y, (b&4) ? p1.z : p0.z), dx);
    }
}

// Paint gridlines
template<> void GridPainter<int>::paint() {
     if (!mObject || mPlane <0 || mPlane >= mLocalGrid->getSize()[mDim])
        return;
    float dx = mLocalGrid->getDx();
    Vec3 box[4];
    glColor3f(0.5,0,0);
    
    bool rbox = true;
    bool lines = mLocalGrid->getSize().max() <= 40; 
lines=true; // debug
    if (lines) {
        //glDepthFunc(GL_LESS);
        glBegin(GL_LINES);
        FOR_P_SLICE(mLocalGrid, mDim, mPlane) {

			int flag = 0;
			flag = mLocalGrid->get(p);

			if (flag & FlagGrid::TypeObstacle) {
    			glColor3f(0.2,0.2,0.2);
			} else if (flag & FlagGrid::TypeOutflow) {
    			glColor3f(0.9,0.2,0);
			} else if (flag & FlagGrid::TypeEmpty) {
    			glColor3f(0.25,0,0);
			} else if (flag & FlagGrid::TypeFluid) {
    			glColor3f(0,0,0.75);
			} else {
    			glColor3f(0.5,0,0); // unknown
			}

            getCellCoordinates(p, box, mDim); 
            for (int n=1;n<=8;n++)
                glVertex(box[(n/2)%4], dx);
        }
        glEnd();
        //glDepthFunc(GL_ALWAYS);        
    }
    
    if (rbox) {
        Vec3 p0(0.0), p1(toVec3(mLocalGrid->getSize())),p(p0);
        glDepthFunc(GL_LESS);
        glBegin(GL_LINES);
        glBox(p0,p1,dx);
        glEnd();
        glDepthFunc(GL_ALWAYS);        
    }
}

// Paint box colors
template<> void GridPainter<Real>::paint() {
    if (!mObject || mHide || mHideLocal || mPlane <0 || mPlane >= mLocalGrid->getSize()[mDim] || !mFlags || !(*mFlags))
        return;
    
    FlagGrid *flags = *mFlags;
    if (flags->getSize() != mLocalGrid->getSize()) flags = 0;
    float dx = mLocalGrid->getDx();
    Vec3 box[4];
    glBegin(GL_QUADS);
    Real scaler = 1.0/getScale();
    bool isLevelset = mLocalGrid->getType() & GridBase::TypeLevelset;
    //glPolygonOffset(1.0,1.0);
    //glDepthFunc(GL_LESS);

	/*FOR_P_SLICE(mLocalGrid, mDim, mPlane) { 
		int flag = FlagGrid::TypeFluid;
		if (flags && (mLocalGrid->getType() & GridBase::TypeLevelset) == 0) flag = flags->get(p);
		if (flag & FlagGrid::TypeObstacle)
			glColor3f(0.15,0.15,0.15);
		else if (flag & FlagGrid::TypeOutflow)
			glColor3f(0.3,0.0,0.0);
		else if (flag & FlagGrid::TypeEmpty)
			glColor3f(0.,0.2,0.);
		else {
			Real v = mLocalGrid->get(p) * scaler;
			
			if (isLevelset) {
				v = max(min(v*0.2, 1.0),-1.0);
				if (v>=0)
					glColor3f(v,0,0.5);
				else
					glColor3f(0.5, 1.0+v, 0.);
			} else {
				if (v>0)
					glColor3f(v,0,0);
				else
					glColor3f(0,0,-v);
			}
		}
		
		if ((flag & FlagGrid::TypeEmpty) == 0) {
			getCellCoordinates(p, box, mDim);
			for (int n=0;n<4;n++) 
				glVertex(box[n], dx);
		}
	}
	glEnd();    */

	// ignore flags, its a bit dangerous :)

	FOR_P_SLICE(mLocalGrid, mDim, mPlane) 
	{ 
		Real v = mLocalGrid->get(p) * scaler; 
		if (isLevelset) {
			v = max(min(v*0.2, 1.0),-1.0);
			if (v>=0)
				glColor3f(v,0,0.5);
			else
				glColor3f(0.5, 1.0+v, 0.);
		} else {
			if (v>0)
				glColor3f(v,0,0);
			else
				glColor3f(0,0,-v);
		}

		getCellCoordinates(p, box, mDim);
		for (int n=0;n<4;n++) 
			glVertex(box[n], dx);
	}
	glEnd();    

    //glDepthFunc(GL_ALWAYS);    
    //glPolygonOffset(0,0);    
}

// Paint velocity vectors
template<> void GridPainter<Vec3>::paint() {
    if (!mObject || mHide || mHideLocal || mPlane <0 || mPlane >= mLocalGrid->getSize()[mDim])
        return;
    
    float dx = mLocalGrid->getDx();
    bool mac = mLocalGrid->getType() & GridBase::TypeMAC;
    Real scaler = getScale();
    glBegin(GL_LINES);
        
    FOR_P_SLICE(mLocalGrid, mDim, mPlane) {        
        Vec3 vel = mLocalGrid->get(p) * scaler;
        Vec3 pos (p.x+0.5, p.y+0.5, p.z+0.5);
        if (mCentered) {
            if (mac) {
                if (p.x < mLocalGrid->getSizeX()-1) 
                    vel.x = 0.5 * (vel.x + scaler * mLocalGrid->get(p.x+1,p.y,p.z).x);
                if (p.y < mLocalGrid->getSizeY()-1) 
                    vel.y = 0.5 * (vel.y + scaler * mLocalGrid->get(p.x,p.y+1,p.z).y);
                if (p.z < mLocalGrid->getSizeZ()-1) 
                    vel.z = 0.5 * (vel.z + scaler * mLocalGrid->get(p.x,p.y,p.z+1).z);
            }
            glColor3f(1,1,1);
            glVertex(pos, dx);
            glColor3f(1,1,1);
            glVertex(pos+vel*1.2, dx);
        } else {
            for (int d=0; d<3; d++) {
                if (fabs(vel[d]) < 1e-2) continue;
                Vec3 p1(pos);
                if (mac)
                    p1[d] -= 0.5f;
                Vec3 color(0.0);
                color[d] = 1;
                glColor3f(color.x, color.y, color.z);
                glVertex(p1, dx);
                glColor3f(1,1,1);
                p1[d] += vel[d];
                glVertex(p1, dx);
            }
        }
    }
    glEnd();    
    
}


// explicit instantiation
template class GridPainter<int>;
template class GridPainter<Real>;
template class GridPainter<Vec3>;
    
} // namespace
