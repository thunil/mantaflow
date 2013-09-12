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

#include "flip.h"
#include "levelset.h"

using namespace std;
namespace Manta {

//! Set velocities from grid with given PIC/FLIP mixture
KERNEL(pts) 
void CopyVelocitiesFromGrid(FlipSystem& p, FlagGrid& flags, MACGrid& vel, MACGrid& oldVel, Real flipRatio) {
    unusedParameter(flags);
    
    if (!p.isActive(i)) return;
    
    /*if (!flags.isFluid(p[i].pos)) {
        p[i].flag |= ParticleBase::PDELETE;
        return;
    }*/
    
    Vec3 v = vel.getInterpolated(p[i].pos);
    Vec3 delta = v - oldVel.getInterpolated(p[i].pos);
    
    p[i].vel = flipRatio * (p[i].vel + delta) + (1.0f - flipRatio) * v;    
}

//! Set velocities on the grid from the particle system
KERNEL(pts, single) 
void CopyVelocitiesToGrid(FlipSystem& p, FlagGrid& flags, MACGrid& vel, Grid<Vec3>& tmp) {
    unusedParameter(flags);
    
    if (!p.isActive(i)) return;
    
    vel.setInterpolated(p[i].pos, p[i].vel, &tmp[0]);
}

void FlipSystem::velocitiesFromGrid(FlagGrid& flags, MACGrid& vel, Real flipRatio) {
    //assertMsg(vel.is3D(), "Only 3D grids supported so far");
    CopyVelocitiesFromGrid(*this, flags, vel, mOldVel, flipRatio);
}

void FlipSystem::velocitiesToGrid(FlagGrid& flags, MACGrid& vel) {
    //assertMsg(vel.is3D(), "Only 3D grids supported so far");
    
    // interpol -> grid. tmpgrid for particle contribution weights
    Grid<Vec3> tmp(mParent);
    vel.clear();
    CopyVelocitiesToGrid(*this, flags, vel, tmp);
	// NT_DEBUG, note - stomp small values in tmp to zero? 
    vel.safeDivide(tmp);
    
    // store diff
    mOldVel = vel;
}

void FlipSystem::initialize(FlagGrid& flags, int discretization, Real randomness ) {
	bool is3D = flags.is3D();
    Real jlen = randomness / discretization;
    Vec3 disp (1.0 / discretization, 1.0 / discretization, 1.0/discretization);
 
    clear(); 
    FOR_IJK(flags) {
        if (flags.isFluid(i,j,k)) {
            Vec3 pos (i,j,k);
            for (int dk=0; dk<(is3D ? discretization : 1); dk++)
            for (int dj=0; dj<discretization; dj++)
            for (int di=0; di<discretization; di++) {
                Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
                subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
				if(!is3D) subpos[2] = 0.5;
                add(FlipData(subpos, Vec3::Zero));
            }
        }
    }
}


void FlipSystem::adjustNumber( MACGrid& vel, FlagGrid& flags, int minParticles, int maxParticles, LevelsetGrid* phi ) 
{
	const Real SURFACE_LS = -1.5; // which levelset to use as threshold
    Grid<int> tmp(mParent);
    
    // count particles in cells, and delete excess particles
    for (int i=0; i<(int)mData.size(); i++) {
        if (isActive(i)) {
            Vec3i p = toVec3i(mData[i].pos);
            int num = tmp(p);

			bool atSurface = false;
			if( phi && (phi->getInterpolated(mData[i].pos) > SURFACE_LS) ) atSurface = true;
            
            // dont delete particles in non fluid cells here, the particles are "always right"
            if ( num > maxParticles && (!atSurface) ) {
                mData[i].flag |= PDELETE;
			} else
                tmp(p) = num+1;
        }
    }
    
    compress();
    
    // seed new particles
    FOR_IJK(tmp) {
        int cnt = tmp(i,j,k);
		
		// skip surface
		if( phi && ((*phi)(i,j,k) > SURFACE_LS) ) continue;

        if (flags.isFluid(i,j,k) && cnt < minParticles) {
            for (int m=cnt; m < minParticles; m++) { 
                Vec3 rndPos (i + mRand.getReal(), j + mRand.getReal(), k + mRand.getReal());
                add(FlipData(rndPos, vel.getInterpolated(rndPos)));
            }
        }
    }
}

void FlipSystem::markFluidCells(FlagGrid& flags) {
    // remove all fluid cells
    FOR_IJK(flags) {
        if (flags.isFluid(i,j,k)) {
            flags(i,j,k) = (flags(i,j,k) | FlagGrid::TypeEmpty) & ~FlagGrid::TypeFluid;
        }
    }
    
    // mark all particles in flaggrid as fluid
    for(int i=0;i<(int)mData.size();i++) {
        const Vec3i p = toVec3i(mData[i].pos);
        if (flags.isInBounds(p) && flags.isEmpty(p))
            flags(p) = (flags(p) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
    }
}



// compute simple levelset without interpolation (fast, low quality), to be used during simulation
KERNEL(pts, single) 
void ComputeUnionLevelset(FlipSystem& p, LevelsetGrid& phi, Real radius=1.) {
    if (!p.isActive(i)) return;

	const Vec3 pos = p[i].pos - Vec3(0.5); // offset for centered LS value
	//const Vec3i size = phi.getSize();
    const int xi = (int)pos.x, yi = (int)pos.y, zi = (int)pos.z; 

	int radiusInt  = int(2. * radius) + 1;
	int radiusIntZ = phi.is3D() ? radiusInt : 0;
	for(int zj=zi-radiusIntZ; zj<=zi+radiusIntZ; zj++) 
	for(int yj=yi-radiusInt ; yj<=yi+radiusInt ; yj++) 
	for(int xj=xi-radiusInt ; xj<=xi+radiusInt ; xj++) 
	{
		if (! phi.isInBounds(Vec3i(xj,yj,zj)) ) continue;
		phi(xj,yj,zj) = std::min( phi(xj,yj,zj) , fabs( norm(Vec3(xj,yj,zj)-pos) )-radius );
	}
}
PYTHON void unionParticleLevelset (FlipSystem& p, LevelsetGrid& phi, Real radiusFactor=1.) {
	// make sure we cover at least 1 cell by default (1% safety margin)
	Real radius = 0.5 * (phi.is3D() ? sqrt(3.) : sqrt(2.) ) * (radiusFactor+.01);

	// reset
	FOR_IJK(phi) {
		phi(i,j,k) = radius + VECTOR_EPSILON;
	} 
	ComputeUnionLevelset(p, phi, radius);
}

ParticleBase* FlipSystem::clone() {
    FlipSystem* nm = new FlipSystem(getParent());
    compress();
    
    nm->mData = mData;
    nm->setName(getName());
	this->cloneParticleData(nm);
    return nm;
}

// ----
FlipSystem::~FlipSystem() { };


//*****************************************************************************

// NT_DEBUG , warning - duplicate functions for now, clean up at some point!

// Set velocities on the grid from the particle system

KERNEL(pts, single) 
void mapLinearVec3ToMACGrid( BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, Grid<Vec3>& tmp, 
	ParticleDataImpl<Vec3>& pvel ) 
{
    unusedParameter(flags);
    if (!p.isActive(i)) return;
    vel.setInterpolated( p[i].pos, pvel[i], &tmp[0] );
}

PYTHON void mapLinear_PartToMAC( FlagGrid& flags, MACGrid& vel , MACGrid& velOld , 
		BasicParticleSystem& parts , ParticleDataImpl<Vec3>& partVel ) 
{
    // interpol -> grid. tmpgrid for particle contribution weights
    Grid<Vec3> tmp(flags.getParent());
    vel.clear();
    mapLinearVec3ToMACGrid( parts, flags, vel, tmp, partVel );

	// NT_DEBUG, note - stomp small values in tmp to zero? use Real grid?
    vel.safeDivide(tmp);
    
    // store original state
    velOld = vel;
}


KERNEL(pts, single) template<class T>
void knMapLinear( BasicParticleSystem& p, FlagGrid& flags, Grid<T>& gdst, Grid<T>& gtmp, 
	ParticleDataImpl<T>& psource ) 
{
    unusedParameter(flags);
    if (!p.isActive(i)) return;
	//debMsg("p "<<(i) <<": "<<p[i].pos <<" "<< psource[i] <<" ", 1); // debug
    gdst.setInterpolated( p[i].pos, psource[i], gtmp );

    /*FOR_IJK_BND(gdst, 0) {
		if(gdst(i,j,k)!=0.) {
			debMsg("at "<<Vec3i(i,j,k) <<": "<<gdst(i,j,k) <<" "<<gtmp(i,j,k) <<" ", 1); // debug
		}
	}*/
}

PYTHON void mapLinearReal( FlagGrid& flags, Grid<Real>& gdst , 
		BasicParticleSystem& parts , ParticleDataImpl<Real>& partSrc ) 
{
    Grid<Real> tmp(flags.getParent());
    gdst.clear();
    knMapLinear<Real>( parts, flags, gdst, tmp, partSrc ); 

	// debug out
    /*FOR_IJK_BND(gdst, 0) {
		if(gdst(i,j,k)!=0.) {
			debMsg("at "<<Vec3i(i,j,k) <<": "<<gdst(i,j,k) <<" "<<tmp(i,j,k) <<" ", 1); // debug
		}
	}*/

    gdst.safeDivide(tmp);
}
/*
PYTHON template<class T>
void mapLinear( FlagGrid& flags, Grid<T>& gdst , 
		BasicParticleSystem& parts , ParticleDataImpl<T>& partSrc ) 
{
    Grid<T> tmp(gdst);
    gdst.clear();
    knMapLinear<T>( parts, flags, gdst, tmp, partSrc ); 
    gdst.safeDivide(tmp);
}

template void mapLinear<Real>( FlagGrid& flags, Grid<Real>& gdst , BasicParticleSystem& parts , ParticleDataImpl<Real>& partVel );
*/


// Get velocities from grid

KERNEL(pts) 
void mapLinearMACGridToVec3( BasicParticleSystem& p, FlagGrid& flags, MACGrid& vel, ParticleDataImpl<Vec3>& pvel ) 
{
    if (!p.isActive(i)) return;
    pvel[i] = vel.getInterpolated( p[i].pos );
}

PYTHON void mapLinear_MACToPart(FlagGrid& flags, MACGrid& vel , 
		BasicParticleSystem& parts , ParticleDataImpl<Vec3>& partVel ) {
    mapLinearMACGridToVec3( parts, flags, vel, partVel );
}

// init

PYTHON void sampleLevelsetWithParticles( LevelsetGrid& phi, FlagGrid& flags, BasicParticleSystem& parts, 
		int discretization, Real randomness ) 
{
	bool is3D = phi.is3D();
    Real jlen = randomness / discretization;
    Vec3 disp (1.0 / discretization, 1.0 / discretization, 1.0/discretization);
    RandomStream mRand(9832);
 
    //clear(); 

    FOR_IJK_BND(phi, 0) {
        if ( flags.isObstacle(i,j,k) ) continue;
        if ( phi(i,j,k) < 1.733 ) {
            Vec3 pos (i,j,k);
            for (int dk=0; dk<(is3D ? discretization : 1); dk++)
            for (int dj=0; dj<discretization; dj++)
            for (int di=0; di<discretization; di++) {
                Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
                subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
				if(!is3D) subpos[2] = 0.5; 
				if( phi.getInterpolated(subpos) > 0. ) continue; 
                parts.add( BasicParticleData(subpos) );
            }
        }
    }

	debMsg(" NT_DEBUG "<< parts.debugInfoPdata() , 1);
}

PYTHON void markFluidCells(BasicParticleSystem& parts, FlagGrid& flags) {
    // remove all fluid cells
    FOR_IJK(flags) {
        if (flags.isFluid(i,j,k)) {
            flags(i,j,k) = (flags(i,j,k) | FlagGrid::TypeEmpty) & ~FlagGrid::TypeFluid;
        }
    }
    
    // mark all particles in flaggrid as fluid
    for(int i=0;i<parts.size();i++) {
    	if (!parts.isActive(i)) continue;
        const Vec3i p = toVec3i( parts.getPos(i) ); // NT_DEBUG, check -0.5 offset?
        if (flags.isInBounds(p) && flags.isEmpty(p))
            flags(p) = (flags(p) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
    }
}

} // namespace

