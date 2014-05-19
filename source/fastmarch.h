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

#ifndef _FASTMARCH_H
#define _FASTMARCH_H

#include <queue>
#include "levelset.h"

namespace Manta {
    
//! Fast marching. Transport certain values
template<class T>
class FmValueTransport {
public:
    FmValueTransport() : mpVel(0), mpFlags(0) { };
    ~FmValueTransport() { };

    void initMarching(MACGrid* vel, FlagGrid* flags) {
        mpVel = vel;
        mpFlags = flags;
    }
    
    inline bool isInitialized() {
        return mpVel != 0;
    }

    //! cell is touched by marching from source cell
    inline void transpTouch(int x,int y,int z, Real *weights, Real time) {
        if(!mpVel || !mpFlags->isEmpty(x,y,z)) return;
        
        T val = T(0.); 
        if(weights[0]>0.0) val += mpVel->get(x+1, y+0, z+0) * weights[0];
        if(weights[1]>0.0) val += mpVel->get(x-1, y+0, z+0) * weights[1];
        if(weights[2]>0.0) val += mpVel->get(x+0, y+1, z+0) * weights[2];
        if(weights[3]>0.0) val += mpVel->get(x+0, y-1, z+0) * weights[3];
		if(mpVel->is3D()) {
        	if(weights[4]>0.0) val += mpVel->get(x+0, y+0, z+1) * weights[4];
        	if(weights[5]>0.0) val += mpVel->get(x+0, y+0, z-1) * weights[5];
		}
        
        // set velocity components if adjacent is empty
        if (mpFlags->isEmpty(x-1,y,z)) (*mpVel)(x,y,z).x = val.x;
        if (mpFlags->isEmpty(x,y-1,z)) (*mpVel)(x,y,z).y = val.y;
		if(mpVel->is3D()) {
        	if (mpFlags->isEmpty(x,y,z-1)) (*mpVel)(x,y,z).z = val.z;            
		}
    }; 

protected:
    MACGrid* mpVel;
    FlagGrid* mpFlags;
};

class FmHeapEntryOut {
public:
    Vec3i p;
    // quick time access for sorting
    Real time;
    static inline bool compare(const Real x, const Real y) { 
        return x > y;
    }

    inline bool operator< (const FmHeapEntryOut& o) const {
        const Real d = fabs((time) - ((o.time)));
        if (d > 0.) return (time) > ((o.time)); 
        if (p.z != o.p.z) return p.z > o.p.z;
        if (p.y != o.p.y) return p.y > o.p.y;
        return p.x > o.p.x;
    };

};

class FmHeapEntryIn {
public:
    Vec3i p;
    // quick time access for sorting
    Real time;
    static inline bool compare(const Real x, const Real y) { 
        return x < y;
    }

    inline bool operator< (const FmHeapEntryIn& o) const {
        const Real d = fabs((time) - ((o.time)));
        if (d > 0.) return (time) < ((o.time)); 
        if (p.z != o.p.z) return p.z < o.p.z;
        if (p.y != o.p.y) return p.y < o.p.y;
        return p.x < o.p.x;
    };
};


//! fast marching algorithm wrapper class
template<class T, int TDIR>
class FastMarch {

public:
    // MSVC doesn't allow static const variables in template classes
    static inline Real InvalidTime() { return -1000; }
    static inline Real InvtOffset() { return 500; }

    enum SpecialValues { FlagInited = 1, FlagIsOnHeap = 2};

    FastMarch(FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& levelset, Real maxTime, MACGrid* velTransport = NULL); 
    ~FastMarch() {}
    
    //! advect level set function with given velocity */
    void performMarching();

    //! test value for invalidity
    inline bool isInvalid(Real v) const { return (v <= InvalidTime()); }

    void addToList(const Vec3i& p, const Vec3i& src);

    //! convert phi to time value
    inline Real phi2time(Real phival) { return (phival-InvalidTime()+ InvtOffset()) * -1.0; }
    
    //! ... and back
    inline Real time2phi(Real tval) { return (InvalidTime() - InvtOffset() - tval); }

    inline Real _phi(int i, int j, int k) { return mLevelset(i,j,k); }
protected:   
    LevelsetGrid& mLevelset;
    Grid<int>& mFmFlags;
    FlagGrid& mFlags;
    
    FmValueTransport<Vec3> mVelTransport;
    
    //! maximal time to march for
    Real mMaxTime;

    //! fast marching list
    std::priority_queue<T, std::vector<T>, std::less<T> > mHeap;
    Real mReheapVal;

    //! weights for touching points
    Real mWeights[6];

    template<int C> inline Real calcWeights(int& okCnt, int& invcnt, Real* v, const Vec3i& idx);
    
    inline Real calculateDistance(const Vec3i& pos);
};

} // namespace
#endif

