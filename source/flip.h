/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * FLIP (fluid implicit particles)
 *
 ******************************************************************************/

#ifndef _FLIP_H
#define _FLIP_H

#include "particle.h"
#include "grid.h"
#include "randomstream.h"

namespace Manta {

class LevelsetGrid;
    
struct FlipData {
    FlipData() : pos(_0),vel(_0),flag(0) {}
    FlipData(const Vec3& p, const Vec3& v) : pos(p),vel(v),flag(0) {}
    Vec3 pos, vel;
    int flag;
    static ParticleBase::SystemType getType() { return ParticleBase::FLIP; }
};

//! FLIP particle system
PYTHON class FlipSystem : public ParticleSystem<FlipData> {
public:
    PYTHON FlipSystem(FluidSolver* parent) : ParticleSystem<FlipData>(parent), mOldVel(parent), mRand(1238943) {}
  
    //! Copy velocities from grid with given PIC/FLIP ratio
    PYTHON void velocitiesFromGrid(FlagGrid& flags, MACGrid& vel, Real flipRatio=0.95);
	//! Write back velocities to grid
    PYTHON void velocitiesToGrid(FlagGrid& flags, MACGrid& vel);
	//! Ensure minimum/maximum number of particles per cell
    PYTHON void adjustNumber(MACGrid& vel, FlagGrid& flags, int minParticles=8, int maxParticles=12, LevelsetGrid* phi=NULL );
    //! Mark cells with particles as fluid, otherwise empty
    PYTHON void markFluidCells(FlagGrid& flags);
    //! Initialize grid of particles + jitter. No surface alignment yet
    PYTHON void initialize(FlagGrid& flags, int discretization=2, Real randomness=0.2 );
    
    virtual ParticleBase* clone();
private:
	MACGrid mOldVel;
    RandomStream mRand;
};

} // namespace

#endif
