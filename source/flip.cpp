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


void FlipSystem::adjustNumber(MACGrid& vel, FlagGrid& flags, int minParticles, int maxParticles) {
    Grid<int> tmp(mParent);
    
    // count particles in cells, and delete excess particles
    for (size_t i=0; i<mData.size(); i++) {
        if (isActive(i)) {
            Vec3i p = toVec3i(mData[i].pos);
            int num = tmp(p);
            
            // dont delete particles in non fluid cells here, the particles are "always right"
            if ( num > maxParticles)
                mData[i].flag |= PDELETE;
            else
                tmp(p) = num+1;
        }
    }
    
    compress();
    
    // seed new particles
    FOR_IJK(tmp) {
        int cnt = tmp(i,j,k);
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


// add simple extrapolation step
PYTHON void extrapolateMACSimple (FlagGrid& flags, MACGrid& vel, int distance = 4) {
    Grid<int> tmp(flags);
	int dim = (flags.is3D() ? 3:2);
	Vec3i nb[6] = { 
		Vec3i(1 ,0,0), Vec3i(-1,0,0),
		Vec3i(0,1 ,0), Vec3i(0,-1,0),
		Vec3i(0,0,1 ), Vec3i(0,0,-1) };

	for(int c=0; c<dim; ++c) {
		Vec3i dir = 0;
		dir[c] = 1;
		tmp.clear();

		// remove all fluid cells
		FOR_IJK_BND(flags,1) {
			Vec3i p(i,j,k);
			if (flags.isFluid(p) || flags.isFluid(p-dir) ) {
				tmp(p) = 1;
			}
		}

		// debug init! , enable for testing only...
		/*FOR_IJK_BND(flags,1) {
			if (tmp(i,j,k) == 0) continue;
			vel(i,j,k)[c] = (i+j+k+c+1.)*0.1;
		}*/
		
		// extrapolate for distance
		for(int d=1; d<1+distance; ++d) {

			FOR_IJK_BND(flags,1) {
				if (tmp(i,j,k) != 0) continue;

				// copy from initialized neighbors
				Vec3i p(i,j,k);
				int nbs = 0;
				Real avgVel = 0.;
				for (int n=0; n<2*dim; ++n) {
					if (tmp(p+nb[n]) == d) {
						//vel(p)[c] = (c+1.)*0.1;
						avgVel += vel(p+nb[n])[c];
						nbs++;
					}
				}

				if(nbs>0) {
					tmp(p)    = d+1;
					vel(p)[c] = avgVel / nbs;
				}
			}

		} // d

	}
}



} // namespace
