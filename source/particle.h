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

namespace Manta {
// fwd decl
template<class T> class Grid;

//! Baseclass for particle systems. Does not implement any data
PYTHON class ParticleBase : public PbClass {
public:
    enum SystemType { BASE, PARTICLE, VELPART, VORTEX, FILAMENT, FLIP };
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
    PYTHON ParticleSystem(FluidSolver* parent);
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
    PYTHON void advectInGrid(FlagGrid& flaggrid, MACGrid& vel, int integrationMode);
    
    /*
     * for FLIP
    //! Sample velocity from grid. Requires vel field
    void velocitiesToGrid(FlagGrid& flaggrid, MACGrid& grid, bool setFlags);
    //! Reproject velocity to grid. Requires vel field
    void velocitiesFromGrid(MACGrid& grid, bool differential);
    */
    
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

} // namespace

#endif