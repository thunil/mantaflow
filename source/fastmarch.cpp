/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Fast marching
 *
 ******************************************************************************/

#include "fastmarch.h"
#include "levelset.h"
#include "kernel.h"
#include <algorithm>

using namespace std;

namespace Manta {
    
template<class COMP, int TDIR>
FastMarch<COMP,TDIR>::FastMarch(FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& levelset, Real maxTime, MACGrid* velTransport)
    : mLevelset(levelset), mFlags(flags), mFmFlags(fmFlags)
{
    if (velTransport)
        mVelTransport.initMarching(velTransport, &flags);
    
    mMaxTime = maxTime * TDIR;
}

// helper for individual components to calculateDistance
template<class COMP, int TDIR> template<int C>
Real FastMarch<COMP,TDIR>::calcWeights(int& okcnt, int& invcnt, Real* v, const Vec3i& idx) {
    Real val = 0.;
    Vec3i idxPlus(idx), idxMinus(idx);
    idxPlus[C]++;
    idxMinus[C]--;
    
    mWeights[C*2] = mWeights[C*2+1] = 0.;
    if (mFmFlags(idxPlus)==FlagInited) {
        // somewhat arbitrary - choose +1 value over -1 ...
        val = mLevelset(idxPlus);
        v[okcnt] = val; okcnt++;
        mWeights[C*2] = 1.;
    } else if (mFmFlags(idxMinus)==FlagInited) {
        val = mLevelset(idxMinus);
        v[okcnt] = val; okcnt++;
        mWeights[C*2+1] = 1.;
    } 
    else {
        invcnt++;
    }
    return val;
}

template<class COMP, int TDIR>
inline Real FastMarch<COMP,TDIR>::calculateDistance(const Vec3i& idx) {
    //int invflag = 0;
    int invcnt = 0;
    Real v[3];
    int okcnt = 0;
    
    Real aVal = calcWeights<0>(okcnt, invcnt, v, idx);
    Real bVal = calcWeights<1>(okcnt, invcnt, v, idx);
    Real cVal = 0.;
	if (mLevelset.is3D())   cVal = calcWeights<2>(okcnt, invcnt, v, idx);
	else					{ invcnt++; mWeights[4] = mWeights[5] = 0.; }

    Real ret = InvalidTime();
    switch(invcnt) {
    case 0: {
        // take all values
        const Real ca=v[0], cb=v[1], cc=v[2];
        const Real csqrt = max(0. , 
                -2.*(ca*ca+cb*cb- cb*cc + cc*cc - ca*(cb+cc)) + 3 );
        // clamp to make sure the sqrt is valid
        ret = 0.333333*( ca+cb+cc+ TDIR*sqrt(csqrt) );

        // weights needed for transport (transpTouch)
        mWeights[0] *= fabs(ret-ca);
        mWeights[1] *= fabs(ret-ca);
        mWeights[2] *= fabs(ret-cb);
        mWeights[3] *= fabs(ret-cb);
        mWeights[4] *= fabs(ret-cc);
        mWeights[5] *= fabs(ret-cc);

        Real norm = 0.0; // try to force normalization
        for(int i=0;i<6;i++) { 
            norm += mWeights[i]; 
        }   
        norm = 1.0/norm;
        for(int i=0;i<6;i++) { mWeights[i] *= norm; } 

        } break; 
    case 1: {
        // take just the 2 ok values
        // t=0.5*( a+b+ (2*g*g-(b-a)*(b-a))^0.5) 
        const Real csqrt = max(0. , 2.-(v[1]-v[0])*(v[1]-v[0]) );
        // clamp to make sure the sqrt is valid
        ret = 0.5*( v[0]+v[1]+ TDIR*sqrt(csqrt) );

        // weights needed for transport (transpTouch)
        mWeights[0] *= fabs(ret-aVal);
        mWeights[1] *= fabs(ret-aVal);
        mWeights[2] *= fabs(ret-bVal);
        mWeights[3] *= fabs(ret-bVal);
        mWeights[4] *= fabs(ret-cVal);
        mWeights[5] *= fabs(ret-cVal);

        Real norm = 0.0; // try to force normalization
        for(int i=0;i<6;i++) { 
            norm += mWeights[i]; 
        }   
        norm = 1.0/norm;
        for(int i=0;i<6;i++) { mWeights[i] *= norm; } 
        // */

        } break; 
    case 2: {
        // just use the one remaining value
        ret = v[0]+ (Real)(TDIR) ; // direction = +- 1
        } break; 
    default:
        throw Error("FastMarch :: Invalid invcnt");
        break;
    }
    return ret;
}

template<class COMP, int TDIR>
void FastMarch<COMP,TDIR>::addToList(const Vec3i& p, const Vec3i& src) {
    if (!mFlags.isInBounds(p,1)) return;
    const int idx = mFlags.index(p);
    
    // already known value?
    // value alreay set to valid value?
    if(mFmFlags[idx] == FlagInited) return;

    // discard by source time now , TODO do instead before calling all addtolists?
    Real srct = mLevelset(src);
    if(COMP::compare(srct, mMaxTime)) return;

    Real ttime = calculateDistance(p);
    
    // remove old entry if larger
    bool found=false;

    Real oldt = mLevelset[idx];
    if (mFmFlags[idx] == FlagIsOnHeap) {
        found = true;
        // is old time better?
        if(COMP::compare(ttime,oldt)) return;        
    }

    // update field
    mFmFlags[idx] = FlagIsOnHeap;
    mLevelset[idx] = ttime;
    
    if (mVelTransport.isInitialized())
        mVelTransport.transpTouch(p.x, p.y, p.z, mWeights, ttime);

    if(!found) {
        // add list entry with source value
        COMP entry;
        entry.p = p;
        entry.time  = &mLevelset[idx];

        mHeap.push( entry );
    }
}

//! Enforce delta_phi = 0 on boundaries
KERNEL(single)
void SetLevelsetBoundaries (LevelsetGrid& phi) {
    if (i==0)      phi(i,j,k) = phi(1,j,k);
    if (i==maxX-1) phi(i,j,k) = phi(i-1,j,k);

    if (j==0)      phi(i,j,k) = phi(i,1,k);
    if (j==maxY-1) phi(i,j,k) = phi(i,j-1,k);

	if(phi.is3D()) {
    	if (k==0)      phi(i,j,k) = phi(i,j,1);
    	if (k==maxZ-1) phi(i,j,k) = phi(i,j,k-1);
	}
}

/*****************************************************************************/
//! Walk...
template<class COMP, int TDIR>
void FastMarch<COMP,TDIR>::performMarching() {
    mReheapVal = 0.0;
    while(mHeap.size() > 0) {
        
        const COMP& ce = mHeap.top();
        Vec3i p = ce.p; 
        mFmFlags(p) = FlagInited;
        mHeap.pop();
        
        addToList(Vec3i(p.x-1,p.y,p.z), p);
        addToList(Vec3i(p.x+1,p.y,p.z), p);
        addToList(Vec3i(p.x,p.y-1,p.z), p);
        addToList(Vec3i(p.x,p.y+1,p.z), p);
		if(mLevelset.is3D()) {
        	addToList(Vec3i(p.x,p.y,p.z-1), p);
        	addToList(Vec3i(p.x,p.y,p.z+1), p);        
		}
    }
    
    // set boundary for plain array
    SetLevelsetBoundaries setls(mLevelset);
}

// explicit instantiation
template class FastMarch<FmHeapEntryIn, -1>;
template class FastMarch<FmHeapEntryOut, +1>;


// a simple extrapolation step , used for cases where there's on levelset
// (note, less accurate than fast marching extrapolation!)
PYTHON void extrapolateMACSimple (FlagGrid& flags, MACGrid& vel, int distance = 4) {
    Grid<int> tmp( flags.getParent() );
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

		// debug init! , enable for testing only - set varying velocities inside
		//FOR_IJK_BND(flags,1) { if (tmp(i,j,k) == 0) continue; vel(i,j,k)[c] = (i+j+k+c+1.)*0.1; }
		
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
