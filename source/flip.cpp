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
    
    // interpol -> grid. tmpgrid for counting
    Grid<Vec3> tmp(mParent);
    vel.clear();
    CopyVelocitiesToGrid(*this, flags, vel, tmp);
	// NT_DEBUG, note - stomp small values in tmp to zero? use Real grid?
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
    for (size_t i=0; i<mData.size(); i++) {
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
    for(size_t i=0;i<mData.size();i++) {
        const Vec3i p = toVec3i(mData[i].pos);
        if (flags.isInBounds(p) && flags.isEmpty(p))
            flags(p) = (flags(p) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
    }
}

ParticleBase* FlipSystem::clone() {
    FlipSystem* nm = new FlipSystem(getParent());
    compress();
    
    nm->mData = mData;
    nm->setName(getName());
    return nm;
}


// compute simple levelset without interpolation (fast, low quality), to be used during simulation
KERNEL(pts, single) 
void ComputeUnionLevelset(FlipSystem& p, LevelsetGrid& phi, Real radius=1.) {
    if (!p.isActive(i)) return;

	const Vec3 pos = p[i].pos - Vec3(0.5); // offset for centered LS value
	const Vec3i size = phi.getSize();
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


// ----


PbClass* FlipSystem::create(PbType t, const string& name) {        
    _args.add("nocheck",true);
    if (t.str == "")
        errMsg("Specify particle data type to create");
    
    PbClass* pyObj = PbClass::createPyObject(t.str, name, _args, this->getParent() );

	ParticleDataBase* pdata = dynamic_cast<ParticleDataBase*>(pyObj);
	if(!pdata) {
		errMsg("Unable to get particle data pointer from newly created object. Only create ParticleData type with a flipSystem.creat() call, eg, PdataReal, PdataVec3 etc.");
		delete pyObj;
		return NULL;
	} else {
		pdata->setParticleSys(this);
		mPartData.push_back(pdata);
		debMsg("ok! " << pdata->size() , 1);
	}

	return pyObj;
}

ParticleDataBase::ParticleDataBase(FluidSolver* parent) : 
		PbClass(parent) , mpParticleSys(NULL) {
};

template<class T>
ParticleDataImpl<T>::ParticleDataImpl(FluidSolver* parent) : ParticleDataBase(parent) {
}

template<class T>
ParticleDataImpl<T>::~ParticleDataImpl() {
}

template<class T>
int  ParticleDataImpl<T>::size() {
	return mData.size();
}
template<class T>
void ParticleDataImpl<T>::add() {
}
template<class T>
void ParticleDataImpl<T>::kill(int i) {
}

// explicit instantiation
template class ParticleDataImpl<int>;
template class ParticleDataImpl<Real>;
template class ParticleDataImpl<Vec3>;

} // namespace
