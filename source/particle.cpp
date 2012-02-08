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

#include "particle.h"
#include "grid.h"
#include "kernel.h"
#include "integrator.h"

#ifdef MESHCODE
    #include "vortexpart.h"
#endif

using namespace std;
namespace Manta {

//******************************************************************************
// ParticleSystem
//******************************************************************************
   
const int DELETE_PART = 20; // chunk size for compression
   
template<class S>
ParticleSystem<S>::ParticleSystem(FluidSolver* parent) 
    : ParticleBase(parent), mDeletes(0), mDeleteChunk(0) 
{ }

template<class S>
void ParticleSystem<S>::clear() {
    mDeleteChunk = mDeletes = 0;
    mData.clear();    
    mSize = mData.size();    
}

template<class S>
int ParticleSystem<S>::add(const S& e) {
    mData.push_back(e); 
    mDeleteChunk = mData.size() / DELETE_PART;
    mSize = mData.size();
    return mData.size()-1;
}

DefineIntegrator(integrateMeshMAC, MACGrid, getInterpolated);

KERNEL(pts) template<class S, IntegrationMode mode>
KnAdvectInGrid(ParticleSystem<S>& p, MACGrid& vel, FlagGrid& flaggrid, Real dt) {
    if (!p.isActive(i)) return;
    
    // from integrator.h
    p[i].pos += integrateMeshMAC<mode>(p[i].pos, vel, dt);
    
    // TODO: else if(flaggrid.isObstacle(pos)) reproject
    if (!flaggrid.isInBounds(p[i].pos,1)) 
        p[i].flag |= ParticleBase::PDELETE;
}

// advection plugin
template<class S>
void ParticleSystem<S>::advectInGrid(FlagGrid& flaggrid, MACGrid& vel, int integrationMode) {
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

/* for FLIP 
 * 
//! Set velocities on a particle system from grid
KERNEL(pts) template<class S>
CopyVelocitiesFromGrid(ParticleSystem<S>& p, MACGrid& vel, bool differential) {
    if (!p.isActive(i)) return;
    
    Vec3 v = vel.getInterpolated(p[i].pos);
    if (differential)
        p[i].vel += v;
    else
        p[i].vel = v;
}

template<class S>
void ParticleSystem<S>::velocitiesFromGrid(MACGrid& grid, bool differential) {
    CopyVelocitiesFromGrid<S>(*this, grid, differential);
}

//! Set velocities on the grid from the particle system
KERNEL(pts) template<class S> 
CopyVelocitiesToGrid(ParticleSystem<S>& p, MACGrid& vel, Grid<Vec3>& tmp, FlagGrid& flaggrid, bool setFlags) {
    if (!p.isActive(i)) return;
    
    vel.setInterpolated(p[i].pos, p[i].vel, &tmp[0]);
    if (setFlags) 
        flaggrid(toVec3i(p[i].pos)) = FlagGrid::TypeFluid;        
}

template<class S>
void ParticleSystem<S>::velocitiesToGrid(FlagGrid& flaggrid, MACGrid& grid, bool setFlags) {
    Grid<Vec3> tmp(mParent);
    CopyVelocitiesToGrid<S>(*this, grid, tmp, flaggrid, setFlags);
    grid.safeDivide(tmp);
}
*/

#ifdef MESHCODE
    // explicit instantiation
    template class ParticleSystem<VortexParticleData>;
#endif

} // namespace
