/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for pressure correction:
 * - solve_pressure
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"

using namespace std;

namespace Manta { 

static inline bool isNotFluid(FlagGrid& flags, int i, int j, int k)
{
    if ( flags.isFluid(i,j,k)   ) return false;
    if ( flags.isFluid(i-1,j,k) ) return false; 
    if ( flags.isFluid(i,j-1,k) ) return false; 
	if ( flags.is3D() ) {
    	if ( flags.isFluid(i,j,k-1) ) return false;
	}
	return true;
}

//! Semi-Lagrange interpolation kernel
KERNEL(bnd=1) template<class T> 
void SemiLagrange (FlagGrid& flags, MACGrid& vel, Grid<T>& dst, Grid<T>& src, Real dt, bool isLevelset) 
{
    if (flags.isObstacle(i,j,k)) {
        dst(i,j,k) = 0;
        return;
    }
    if (!isLevelset && isNotFluid(flags,i,j,k) ) {
        dst(i,j,k) = src(i,j,k);
        return;
    }
    
    // SL traceback
    Vec3 pos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getCentered(i,j,k) * dt;
    dst(i,j,k) = src.getInterpolated(pos);
}

static inline bool isNotFluidMAC(FlagGrid& flags, int i, int j, int k)
{
    if ( flags.isFluid(i,j,k)   ) return false;
	return true;
}

//! Semi-Lagrange interpolation kernel for MAC grids
KERNEL(bnd=1)
void SemiLagrangeMAC(FlagGrid& flags, MACGrid& vel, MACGrid& dst, MACGrid& src, Real dt) 
{
    if (flags.isObstacle(i,j,k)) {
        dst(i,j,k) = 0;
        return;
    }
    if ( isNotFluidMAC(flags,i,j,k) ) {
        dst(i,j,k) = src(i,j,k);
        return;
    }
    
    // get currect velocity at MAC position
    // no need to shift xpos etc. as lookup field is also shifted
    Vec3 xpos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACX(i,j,k) * dt;
    Real vx = src.getInterpolatedComponent<0>(xpos);
    Vec3 ypos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACY(i,j,k) * dt;
    Real vy = src.getInterpolatedComponent<1>(ypos);
    Vec3 zpos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACZ(i,j,k) * dt;
    Real vz = src.getInterpolatedComponent<2>(zpos);
    
    dst(i,j,k) = Vec3(vx,vy,vz);
}

//! Kernel: Correct based on forward and backward SL steps (for both centered & mac grids)
KERNEL(idx) template<class T> 
void MacCormackCorrect(FlagGrid& flags, Grid<T>& dst, Grid<T>& old, Grid<T>& fwd,  Grid<T>& bwd, 
                       Real strength, bool isLevelSet, bool isMAC=false )
{
	// note, replacement for isNotFluidMAC and isNotFluid
	bool skip = false;
    if (!flags.isFluid(idx                   ) )     skip = true;
	if(!isMAC) {
    if (!flags.isFluid(idx-flags.getStrideX()) )     skip = true; 
    if (!flags.isFluid(idx-flags.getStrideY()) )     skip = true; 
	if ( flags.is3D() ) {
    	if (!flags.isFluid(idx-flags.getStrideY()) ) skip = true;
	} }
    if ( skip ) {
        dst[idx] = isLevelSet ? fwd[idx] : (T)0.0;
        return;
    }
    
    // note, strenth of correction can be modified here
	dst[idx] = fwd[idx] + strength * 0.5 * (old[idx] - bwd[idx]);
}

// Helper to collect min/max in a template
template<class T> inline void getMinMax(T& minv, T& maxv, const T& val) {
    if (val < minv) minv = val;
    if (val > maxv) maxv = val;
}
template<> inline void getMinMax<Vec3>(Vec3& minv, Vec3& maxv, const Vec3& val) {
    getMinMax(minv.x, maxv.x, val.x);
    getMinMax(minv.y, maxv.y, val.y);
    getMinMax(minv.z, maxv.z, val.z);
}

    
//! Helper function for clamping non-mac grids
template<class T>
inline T doClampComponent(const Vec3i& upperClamp, Grid<T>& orig, T dst, const Vec3i& posFwd) {
    // clamp forward lookup to grid
    const int i0 = clamp(posFwd.x, 0, upperClamp.x-1);
    const int j0 = clamp(posFwd.y, 0, upperClamp.y-1);
    const int k0 = clamp(posFwd.z, 0, (orig.is3D() ? (upperClamp.z-1) : 1) );
    const int i1 = i0+1, j1 = j0+1, k1= (orig.is3D() ? (k0+1) : k0);
    
    if (!orig.isInBounds(Vec3i(i0,j0,k0),1)) {
	return dst;
    }

        // find min/max around fwd pos
        T minv = orig(i0,j0,k0), maxv = minv;
        getMinMax(minv, maxv, orig(i1,j0,k0));
        getMinMax(minv, maxv, orig(i0,j1,k0));
        getMinMax(minv, maxv, orig(i1,j1,k0));
        getMinMax(minv, maxv, orig(i0,j0,k1));
        getMinMax(minv, maxv, orig(i1,j0,k1));
        getMinMax(minv, maxv, orig(i0,j1,k1));
        getMinMax(minv, maxv, orig(i1,j1,k1)); 
        
        // write clamped value
        return clamp(dst, minv, maxv);
}
    
//! Helper function for clamping MAC grids
template<int c> 
inline Real doClampComponentMAC(const Vec3i& upperClamp, MACGrid& orig, Real dst, const Vec3i& posFwd) {
    // clamp forward lookup to grid
    const int i0 = clamp(posFwd.x, 0, upperClamp.x-1);
    const int j0 = clamp(posFwd.y, 0, upperClamp.y-1);
    const int k0 = clamp(posFwd.z, 0, (orig.is3D() ? (upperClamp.z-1) : 1) );
    const int i1 = i0+1, j1 = j0+1, k1= (orig.is3D() ? (k0+1) : k0);
    if (!orig.isInBounds(Vec3i(i0,j0,k0),1)) 
        return dst;
    
    // find min/max around fwd pos
    Real minv = orig(i0,j0,k0)[c], maxv = minv;
    getMinMax(minv, maxv, orig(i1,j0,k0)[c]);
    getMinMax(minv, maxv, orig(i0,j1,k0)[c]);
    getMinMax(minv, maxv, orig(i1,j1,k0)[c]);
    getMinMax(minv, maxv, orig(i0,j0,k1)[c]);
    getMinMax(minv, maxv, orig(i1,j0,k1)[c]);
    getMinMax(minv, maxv, orig(i0,j1,k1)[c]);
    getMinMax(minv, maxv, orig(i1,j1,k1)[c]);
    
    return clamp(dst, minv, maxv);    
}

//! Kernel: Clamp obtained value to min/max in source area, and reset values that point out of grid or into boundaries
//          (note - MAC grids are handled below)
KERNEL(bnd=1) template<class T>
void MacCormackClamp(FlagGrid& flags, MACGrid& vel, Grid<T>& dst, Grid<T>& orig, Grid<T>& fwd, Real dt)
{
    if (flags.isObstacle(i,j,k))
        return;
    if ( isNotFluid(flags,i,j,k) ) {
        dst(i,j,k) = fwd(i,j,k);
        return;
    }

    T     dval       = dst(i,j,k);
    Vec3i upperClamp = flags.getSize() - 1;
    
    // lookup forward/backward
    Vec3i posFwd = toVec3i( Vec3(i,j,k) - vel.getCentered(i,j,k) * dt );
    Vec3i posBwd = toVec3i( Vec3(i,j,k) + vel.getCentered(i,j,k) * dt );
    
    dval = doClampComponent<T>(upperClamp, orig, dval, posFwd );
    
    // test if lookups point out of grid or into obstacle
    if (posFwd.x < 0 || posFwd.y < 0 || posFwd.z < 0 ||
        posBwd.x < 0 || posBwd.y < 0 || posBwd.z < 0 ||
        posFwd.x > upperClamp.x || posFwd.y > upperClamp.y || ((posFwd.z > upperClamp.z)&&flags.is3D()) ||
        posBwd.x > upperClamp.x || posBwd.y > upperClamp.y || ((posBwd.z > upperClamp.z)&&flags.is3D()) ||
        flags.isObstacle(posFwd) || flags.isObstacle(posBwd) ) 
    {        
        dval = fwd(i,j,k);
    }
    dst(i,j,k) = dval;
}

//! Kernel: same as MacCormackClamp above, but specialized version for MAC grids
KERNEL(bnd=1) 
void MacCormackClampMAC (FlagGrid& flags, MACGrid& vel, MACGrid& dst, MACGrid& orig, MACGrid& fwd, Real dt)
{
    if (flags.isObstacle(i,j,k))
        return;
    if ( isNotFluidMAC(flags,i,j,k) ) {
        dst(i,j,k) = fwd(i,j,k);
        return;
    }
    
    Vec3  pos(i,j,k);
    Vec3  dval       = dst(i,j,k);
    Vec3i upperClamp = flags.getSize() - 1;
    
    // get total fwd lookup
    Vec3i posFwd = toVec3i( Vec3(i,j,k) - vel.getCentered(i,j,k) * dt );
    Vec3i posBwd = toVec3i( Vec3(i,j,k) + vel.getCentered(i,j,k) * dt );
    
    // clamp individual components
    dval.x = doClampComponentMAC<0>(upperClamp, orig, dval.x, toVec3i( pos - vel.getAtMACX(i,j,k) * dt) );
    dval.y = doClampComponentMAC<1>(upperClamp, orig, dval.y, toVec3i( pos - vel.getAtMACY(i,j,k) * dt) );
    dval.z = doClampComponentMAC<2>(upperClamp, orig, dval.z, toVec3i( pos - vel.getAtMACZ(i,j,k) * dt) );
    
    // test if lookups point out of grid or into obstacle
    if (posFwd.x < 0 || posFwd.y < 0 || posFwd.z < 0 ||
        posBwd.x < 0 || posBwd.y < 0 || posBwd.z < 0 ||
        posFwd.x > upperClamp.x || posFwd.y > upperClamp.y || ((posFwd.z > upperClamp.z)&&flags.is3D()) ||
        posBwd.x > upperClamp.x || posBwd.y > upperClamp.y || ((posBwd.z > upperClamp.z)&&flags.is3D()) 
        //|| flags.isObstacle(posFwd) || flags.isObstacle(posBwd)  // note - this unfortunately introduces asymmetry... TODO update
		) 
    {        
        dval = fwd(i,j,k);
    }
 
    // writeback
    dst(i,j,k) = dval;
}

//! template function for performing SL advection
template<class GridType> 
void fnAdvectSemiLagrange(FluidSolver* parent, FlagGrid& flags, MACGrid& vel, GridType& orig, int order, Real strength) {
    typedef typename GridType::BASETYPE T;
    
    Real dt = parent->getDt();
    bool levelset = orig.getType() & GridBase::TypeLevelset;
    
    // forward step
    GridType fwd(parent);
    SemiLagrange<T> (flags, vel, fwd, orig, dt, levelset);
    
    if (order == 1) {
        orig.swap(fwd);
    }
    else if (order == 2) { // MacCormack
        GridType bwd(parent);
        GridType newGrid(parent);
    
        // bwd <- backwards step
        SemiLagrange<T> (flags, vel, bwd, fwd, -dt, levelset);
        
        // newGrid <- compute correction
        MacCormackCorrect<T> (flags, newGrid, orig, fwd, bwd, strength, levelset);
        
        // clamp values
        MacCormackClamp<T> (flags, vel, newGrid, orig, fwd, dt);
        
        orig.swap(newGrid);
    }
}

//! template function for performing SL advection: specialized version for MAC grids
template<> 
void fnAdvectSemiLagrange<MACGrid>(FluidSolver* parent, FlagGrid& flags, MACGrid& vel, MACGrid& orig, int order, Real strength) {
    Real dt = parent->getDt();
    
    // forward step
    MACGrid fwd(parent);    
    SemiLagrangeMAC (flags, vel, fwd, orig, dt);
    
    if (order == 1) {
        orig.swap(fwd);
    }
    else if (order == 2) { // MacCormack
        MACGrid bwd(parent);
        MACGrid newGrid(parent);
        
        // bwd <- backwards step
        SemiLagrangeMAC (flags, vel, bwd, fwd, -dt);
        
        // newGrid <- compute correction
        MacCormackCorrect<Vec3> (flags, newGrid, orig, fwd, bwd, strength, false, true);
        
        // clamp values
        MacCormackClampMAC (flags, vel, newGrid, orig, fwd, dt);
        
        orig.swap(newGrid);
    }
}

//! Perform semi-lagrangian advection of target Real- or Vec3 grid
PYTHON void advectSemiLagrange (FlagGrid* flags, MACGrid* vel, GridBase* grid, 
                           int order = 1, Real strength = 1.0)
{    
    assertMsg(order==1 || order==2, "AdvectSemiLagrange: Only order 1 (regular SL) and 2 (MacCormack) supported");
    
    // determine type of grid    
    if (grid->getType() & GridBase::TypeReal) {
        fnAdvectSemiLagrange< Grid<Real> >(parent, *flags, *vel, *((Grid<Real>*) grid), order, strength);
    }
    else if (grid->getType() & GridBase::TypeMAC) {    
        fnAdvectSemiLagrange< MACGrid >(parent, *flags, *vel, *((MACGrid*) grid), order, strength);
    }
    else if (grid->getType() & GridBase::TypeVec3) {    
        fnAdvectSemiLagrange< Grid<Vec3> >(parent, *flags, *vel, *((Grid<Vec3>*) grid), order, strength);
    }
    else
        errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");    
}

} // end namespace DDF 

