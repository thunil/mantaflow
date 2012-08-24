/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Base class for particle systems
 *
 ******************************************************************************/

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <vector>
#include "grid.h"
#include "vectorbase.h"
#include "integrator.h"
namespace Manta {
// fwd decl
template<int DIM, class T> class Grid;

//! Baseclass for particle systems. Does not implement any data
PYTHON class ParticleBase : public PbClass {
public:
    enum SystemType { BASE, PARTICLE, VELPART, VORTEX, FILAMENT, FLIP, TRACER };
    enum ParticleType {
        PNONE         = 0,
        PDELETE       = (1<<10), // mark as deleted, will be deleted in next compress() step
        PINVALID     = (1<<30), // unused
    };
    
    PYTHON ParticleBase(FluidSolver* parent) : PbClass(parent), mSize(0) {}

    virtual SystemType getType() const { return BASE; }
    inline int size() const { return mSize; }    
    virtual ParticleBase* clone() { return NULL; }
    
protected:
    int mSize;
};


//! Main class for particle systems
/*! Basetype S must at least contain flag, pos fields */
PYTHON template<class S> class ParticleSystem : public ParticleBase {
public:    
    PYTHON ParticleSystem(FluidSolver* parent) : ParticleBase(parent), mDeletes(0), mDeleteChunk(0) {}
    virtual ~ParticleSystem() {}
    
    virtual SystemType getType() const { return S::getType(); };
    
    // accessors
    inline S& operator[](int i) { return mData[i]; }
    inline const S& operator[](int i) const { return mData[i]; }
    
    // adding and deleting
    inline void kill(int i) { mData[i].flag |= PDELETE; if (++mDeletes > mDeleteChunk) compress(); }
    inline bool isActive(int i) { return (mData[i].flag & PDELETE) == 0; }    
    int add(const S& data);
    void clear();
    
    // plugins
    
    //! Advect particle in grid velocity field
    PYTHON void advectInGrid(FlagGrid3& flaggrid, MACGrid3& vel, int integrationMode);
    
    virtual ParticleBase* clone();
    
    protected:  
    void compress();
    
    int mDeletes, mDeleteChunk;    
    std::vector<S> mData;    
};

//! Simplest data class for particle systems
struct BasicParticleData {
    Vec3 pos;
    int flag;
    static ParticleBase::SystemType getType() { return ParticleBase::PARTICLE; }
};






//******************************************************************************
// Implementation
//******************************************************************************

const int DELETE_PART = 20; // chunk size for compression
   
template<class S>
void ParticleSystem<S>::clear() {
    mDeleteChunk = mDeletes = 0;
    mData.clear();    
    mSize = mData.size();    
}

template<class S>
int ParticleSystem<S>::add(const S& data) {
    mData.push_back(data); 
    mDeleteChunk = mData.size() / DELETE_PART;
    mSize = mData.size();
    return mData.size()-1;
}

DefineIntegrator(integrateMeshMAC, MACGrid3, getInterpolated);

KERNEL(pts) template<class S, IntegrationMode mode>
KnAdvectInGrid(ParticleSystem<S>& p, MACGrid3& vel, FlagGrid3& flaggrid, Real dt) {
    if (!p.isActive(i)) return;
    
    // from integrator.h
    p[i].pos += integrateMeshMAC<mode>(p[i].pos, vel, dt);
    
    // TODO: else if(flaggrid.isObstacle(pos)) reproject
    if ((!flaggrid.isInBounds(p[i].pos,1) || flaggrid.isObstacle(p[i].pos)) && p[i].pos.x > 5)
        p[i].flag |= ParticleBase::PDELETE;
}

// advection plugin
template<class S>
void ParticleSystem<S>::advectInGrid(FlagGrid3& flaggrid, MACGrid3& vel, int integrationMode) {
    const Real dt = mParent->getDt();
    switch((IntegrationMode)integrationMode) {
        case EULER: KnAdvectInGrid<S, EULER>(*this, vel, flaggrid, dt); break;
        case RK2: KnAdvectInGrid<S, RK2>(*this, vel, flaggrid, dt); break;
        case RK4: KnAdvectInGrid<S, RK4>(*this, vel, flaggrid, dt); break;
        default: throw Error("invalid integration mode");
    }
}

template<class S>
void ParticleSystem<S>::compress() {
    int nextRead = mData.size();
    for (size_t i=0; i<mData.size(); i++) {
        while ((mData[i].flag & PDELETE) != 0) {
            nextRead--;
            mData[i] = mData[nextRead];
            mData[nextRead].flag = 0;           
        }
    }
    mData.resize(nextRead);
    mDeletes = 0;
    mDeleteChunk = mData.size() / DELETE_PART;
    mSize = mData.size();    
}

template<class S>
ParticleBase* ParticleSystem<S>::clone() {
    ParticleSystem<S>* nm = new ParticleSystem<S>(getParent());
    compress();
    
    nm->mData = mData;
    nm->mSize = mSize;
    nm->setName(getName());
    return nm;
}

} // namespace





#endif