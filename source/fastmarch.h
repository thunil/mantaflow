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

#include <vector>
#include "levelset.h"

namespace Manta {
    
//! Fast marching. Transport certain values
template<int DIM, class T>
class FmValueTransport {
public:
    FmValueTransport() : mpFlags(0), mpVel(0) { };
    ~FmValueTransport() { };

    void initMarching(MACGrid<DIM>* vel, FlagGrid<DIM>* flags) {
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
        if(DIM==3 && weights[4]>0.0) val += mpVel->get(x+0, y+0, z+1) * weights[4];
        if(DIM==3 && weights[5]>0.0) val += mpVel->get(x+0, y+0, z-1) * weights[5];
        
        // set velocity components if adjacent is empty
        if (mpFlags->isEmpty(x-1,y,z)) (*mpVel)(x,y,z).x = val.x;
        if (mpFlags->isEmpty(x,y-1,z)) (*mpVel)(x,y,z).y = val.y;
        if (DIM==3 && mpFlags->isEmpty(x,y,z-1)) (*mpVel)(x,y,z).z = val.z;            
    }; 

protected:
    MACGrid<DIM>* mpVel;
    FlagGrid<DIM>* mpFlags;
};

class FmHeapEntry {
public:
    Vec3i p;
    // quick time access for sorting
    Real *time;
};

//! heap comparison object for outwards marching
class FmHeapComparatorOut {
public:
    static inline bool compare(const Real x, const Real y) { 
        return x > y;
    }

    inline bool operator() (const FmHeapEntry& x, const FmHeapEntry& y) const {
        return (*(x.time) > *(y.time));
    };
};

//! heap comparison object for inwards marching
class FmHeapComparatorIn {
public:
    static inline bool compare(const Real x, const Real y) {
        return x < y;
    }

    inline bool operator() (const FmHeapEntry& x, const FmHeapEntry& y) const {
        return (*(x.time) < *(y.time));
    };
};


// MSVC doesn't support static const 
#define InvalidTime -1000.0
#define InvtOffset    500.0

//! fast marching algorithm wrapper class
template<int DIM, class COMP, int TDIR>
class FastMarch {

public:

    enum SpecialValues { FlagInited = 1, FlagIsOnHeap = 2};

    FastMarch(FlagGrid<DIM>& flags, Grid<DIM,int>& fmFlags, LevelsetGrid<DIM>& levelset, Real maxTime, MACGrid<DIM>* velTransport = NULL); 
    ~FastMarch() {}
    
    //! advect level set function with given velocity */
    void performMarching();

    //! test value for invalidity
    inline bool isInvalid(Real v) const { return (v <= InvalidTime); }

    void addToList(const Vec3i& p, const Vec3i& src);

    //! convert phi to time value
    inline Real phi2time(Real phival) { return (phival-InvalidTime+ InvtOffset) * -1.0; }
    
    //! ... and back
    inline Real time2phi(Real tval) { return (InvalidTime - InvtOffset - tval); }

protected:   

    LevelsetGrid<DIM>& mLevelset;
    Grid<DIM,int>& mFmFlags;
    FlagGrid<DIM>& mFlags;
    
    FmValueTransport<DIM, Vec3> mVelTransport;
    
    //! maximal time to march for
    Real mMaxTime;

    //! fast marching list
    std::vector<FmHeapEntry> mHeap;
    COMP mHeapComp;
    Real mReheapVal;

    //! weights for touching points
    Real mWeights[6];

    template<int C> inline Real calcWeights(int& okCnt, int& invcnt, Real* v, const Vec3i& idx);
    
    inline Real calculateDistance(const Vec3i& pos);
};

} // namespace
#endif

