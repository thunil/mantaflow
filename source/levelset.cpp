/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Levelset
 *
 ******************************************************************************/

#include "levelset.h"
#include "fastmarch.h"
#include "kernel.h"

using namespace std;
namespace Manta {

//************************************************************************
// Helper functions and kernels for marching

static const int FlagInited = FastMarch<FmHeapComparatorOut, +1>::FlagInited;

// neighbor lookup vectors
static const Vec3i neighbors[6] = { Vec3i(-1,0,0), Vec3i(1,0,0), Vec3i(0,-1,0), Vec3i(0,1,0), Vec3i(0,0,-1), Vec3i(0,0,1) };
    
KERNEL(bnd=1) InitFmIn (FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& phi, bool ignoreWalls) {
    const int idx = flags.index(i,j,k);
    const Real v = phi[idx];
    if (v>=0 && (!ignoreWalls || !flags.isObstacle(idx)))
        fmFlags[idx] = FlagInited;
    else
        fmFlags[idx] = 0;
}

KERNEL(bnd=1) InitFmOut (FlagGrid& flags, Grid<int>& fmFlags, LevelsetGrid& phi, bool ignoreWalls) {
    const int idx = flags.index(i,j,k);
    const Real v = phi[idx];
    if (ignoreWalls) {
        fmFlags[idx] = (v<0) ? FlagInited : 0;
        if (flags.isObstacle(idx)) {
            fmFlags[idx] = 0;
            phi[idx] = 0;
        }
    }
    else
        fmFlags[idx] = (v<0) ? FlagInited : 0;
}

KERNEL(bnd=1) SetUninitialized (Grid<int>& fmFlags, LevelsetGrid& phi, const Real val) {
    if (fmFlags(i,j,k) != FlagInited)
        phi(i,j,k) = val;
}

template<bool inward>
inline bool isAtInterface(Grid<int>& fmFlags, LevelsetGrid& phi, const Vec3i& p) {
    // check for interface
    for (int nb=0; nb<6; nb++) {
        const Vec3i pn(p + neighbors[nb]);
        if (!fmFlags.isInBounds(pn)) continue;
        
        if (fmFlags(pn) != FlagInited) continue;
        if ((inward && phi(pn) >= 0) || 
            (!inward && phi(pn) < 0)) return true;
    }
    return false;
}

//************************************************************************
// Levelset class def

LevelsetGrid::LevelsetGrid(FluidSolver* parent, bool show) 
    : Grid<Real>(parent, show) 
{ 
    mType = (GridType)(TypeLevelset | TypeReal);    
}    

extern void updateQtGui(bool full, int frame); // HACK

Real LevelsetGrid::invalidTimeValue() {
    return FastMarch<FmHeapComparatorOut, 1>::InvalidTime;
}

void LevelsetGrid::reinitMarching(FlagGrid& flags, Real maxTime, MACGrid* velTransport, bool ignoreWalls, bool correctOuterLayer)
{
    Grid<int> fmFlags(mParent);
    LevelsetGrid& phi = *this;
    
    FastMarch<FmHeapComparatorOut, +1> marchOut(flags, fmFlags, phi, maxTime, velTransport);
    FastMarch<FmHeapComparatorIn, -1> marchIn(flags, fmFlags, phi, maxTime, NULL);
    
    // march inside
    InitFmIn (flags, fmFlags, phi, ignoreWalls);
    
    FOR_IJK_BND(flags, 1) {
        if (fmFlags(i,j,k) == FlagInited) continue;
        if (flags.isObstacle(i,j,k)) continue;
        const Vec3i p(i,j,k);
                
        if(isAtInterface<true>(fmFlags, phi, p)) {
            // set value
            fmFlags(p) = FlagInited;
            
            // add neighbors that are not at the interface
            for (int nb=0; nb<6; nb++) {
                const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                if (flags.isObstacle(pn)) continue;
                
                // check neighbors of neighbor
                if (phi(pn) < 0 && !isAtInterface<true>(fmFlags, phi, pn)) {
                    marchIn.addToList(pn, p);
                }
            }            
        }
    }
    marchIn.performMarching(); 
        
    // set un initialized regions
    SetUninitialized (fmFlags, phi, -maxTime); 
    
    // done with inwards marching, now march out...    
    InitFmOut (flags, fmFlags, phi, ignoreWalls);
    
    // by default, correctOuterLayer is on
    if (correctOuterLayer) {    
        // normal version, inwards march is done, now add all outside values (0..2] to list
        // note, this might move the interface a bit! but keeps a nice signed distance field...        
        FOR_IJK_BND(flags, 1) {
            if (flags.isObstacle(i,j,k)) continue;
            const Vec3i p(i,j,k);
            
            // check nbs
            for (int nb=0; nb<6; nb++) {
                const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                
                if (fmFlags(pn) != FlagInited) continue;
                if (flags.isObstacle(pn)) continue;
                
                const Real nbPhi = phi(pn);
                
                // only add nodes near interface, not e.g. outer boundary vs. invalid region                
                if (nbPhi < 0 && nbPhi >= -2)
                    marchOut.addToList(p, pn);
            }
        }         
    } else {
        // alternative version, keep interface, do not distort outer cells
        // add all ouside values, but not those at the IF layer
        FOR_IJK_BND(flags, 1) {
            if (flags.isObstacle(i,j,k)) continue;
            
            // only look at ouside values
            const Vec3i p(i,j,k);
            if (phi(p) < 0) continue;
            
            if (isAtInterface<false>(fmFlags, phi, p)) {
                // now add all non, interface neighbors
                fmFlags(p) = FlagInited;
                
                // add neighbors that are not at the interface
                for (int nb=0; nb<6; nb++) {
                    const Vec3i pn(p + neighbors[nb]); // index always valid due to bnd=1                
                    if (flags.isObstacle(pn)) continue;
                
                    // check neighbors of neighbor
                    if (phi(pn) > 0 && !isAtInterface<false>(fmFlags, phi, pn))
                        marchOut.addToList(pn, p);
                }            
            }
        }
    }    
    marchOut.performMarching();
    
    // set un initialized regions
    SetUninitialized (fmFlags, phi, +maxTime);    
}

} //namespace