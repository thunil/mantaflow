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
    if(mFmFlags(idxPlus)==FlagInited) {
        // somewhat arbitrary - choose +1 value over -1 ...
        val = mLevelset(idxPlus);
        v[okcnt] = val; okcnt++;
        mWeights[C*2] = 1.;
    } else if(mFmFlags(idxMinus)==FlagInited) {
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
    Real cVal = calcWeights<2>(okcnt, invcnt, v, idx);

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
        //debMsg("RET","a="<<a<<" b="<<b<<" ret="<<ret );

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
    if(mHeapComp.compare(srct, mMaxTime)) return;

    Real ttime = calculateDistance(p);
    
    // remove old entry if larger
    bool found=false;

    Real oldt = mLevelset[idx];
    if (mFmFlags[idx] == FlagIsOnHeap) {
        found = true;
        // is old time better?
        if(mHeapComp.compare(ttime,oldt)) return;
    }

    // update field
    mFmFlags[idx] = FlagIsOnHeap;
    mLevelset[idx] = ttime;

    if (mVelTransport.isInitialized())
        mVelTransport.transpTouch(p.x, p.y, p.z, mWeights, ttime);

    if(!found) {
        // add list entry with source value
        FmHeapEntry entry;
        entry.p = p;
        entry.time  = &mLevelset[idx];

        mHeap.push_back( entry );

        push_heap(mHeap.begin(), mHeap.end(), mHeapComp );        
    }
}

//! Enforce delta_phi = 0 on boundaries
KERNEL 
void SetLevelsetBoundaries (LevelsetGrid& phi) {
    if (i==0) phi(i,j,k) = phi(1,j,k);
    if (j==0) phi(i,j,k) = phi(i,1,k);
    if (k==0) phi(i,j,k) = phi(i,j,1);
    if (i==maxX-1) phi(i,j,k) = phi(i-1,j,k);
    if (j==maxY-1) phi(i,j,k) = phi(i,j-1,k);
    if (k==maxZ-1) phi(i,j,k) = phi(i,j,k-1);
}

/*****************************************************************************/
//! Walk...
template<class COMP, int TDIR>
void FastMarch<COMP,TDIR>::performMarching() {
    
    make_heap(mHeap.begin(), mHeap.end(), mHeapComp ); 
    
    mReheapVal = 0.0;
    while(mHeap.size() > 0) {
        if(mReheapVal<0.0) {
            make_heap(mHeap.begin(), mHeap.end(), mHeapComp ); 
            mReheapVal = 0.0;            
        }

        pop_heap( mHeap.begin(), mHeap.end() , mHeapComp );
        
        Vec3i p = mHeap[mHeap.size()-1].p; 
        mFmFlags(p) = FlagInited;
        mHeap.pop_back();
    
        addToList(Vec3i(p.x-1,p.y,p.z), p);
        addToList(Vec3i(p.x+1,p.y,p.z), p);
        addToList(Vec3i(p.x,p.y-1,p.z), p);
        addToList(Vec3i(p.x,p.y+1,p.z), p);
        addToList(Vec3i(p.x,p.y,p.z-1), p);
        addToList(Vec3i(p.x,p.y,p.z+1), p);        
    }

    // set boundary for plain array
    SetLevelsetBoundaries setls(mLevelset);
}

// explicit instantiation
template class FastMarch<FmHeapComparatorIn, -1>;
template class FastMarch<FmHeapComparatorOut, +1>;

} // namespace
