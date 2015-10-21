/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2015 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for advection
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"
#include <limits>

using namespace std;

namespace Manta { 


//! Semi-Lagrange interpolation kernel
KERNEL(bnd=1) template<class T> 
void SemiLagrange (FlagGrid& flags, MACGrid& vel, Grid<T>& dst, Grid<T>& src, Real dt, bool isLevelset, int orderSpace) 
{
	// traceback position
	Vec3 pos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getCentered(i,j,k) * dt;
	dst(i,j,k) = src.getInterpolatedHi(pos, orderSpace);
}

//! Semi-Lagrange interpolation kernel for MAC grids
KERNEL(bnd=1)
void SemiLagrangeMAC(FlagGrid& flags, MACGrid& vel, MACGrid& dst, MACGrid& src, Real dt, int orderSpace) 
{
	// get currect velocity at MAC position
	// no need to shift xpos etc. as lookup field is also shifted
	Vec3 xpos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACX(i,j,k) * dt;
	Real vx = src.getInterpolatedComponentHi<0>(xpos, orderSpace);
	Vec3 ypos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACY(i,j,k) * dt;
	Real vy = src.getInterpolatedComponentHi<1>(ypos, orderSpace);
	Vec3 zpos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getAtMACZ(i,j,k) * dt;
	Real vz = src.getInterpolatedComponentHi<2>(zpos, orderSpace);
	
	dst(i,j,k) = Vec3(vx,vy,vz);
}


//! Kernel: Correct based on forward and backward SL steps (for both centered & mac grids)
KERNEL(idx) template<class T> 
void MacCormackCorrect(FlagGrid& flags, Grid<T>& dst, Grid<T>& old, Grid<T>& fwd,  Grid<T>& bwd, 
					   Real strength, bool isLevelSet, bool isMAC=false )
{
	dst[idx] = fwd[idx];

	if (flags.isFluid(idx)) {
		// only correct inside fluid region; note, strenth of correction can be modified here
		dst[idx] += strength * 0.5 * (old[idx] - bwd[idx]);
	}
}

//! Kernel: Correct based on forward and backward SL steps (for both centered & mac grids)
KERNEL() template<class T> 
void MacCormackCorrectMAC(FlagGrid& flags, Grid<T>& dst, Grid<T>& old, Grid<T>& fwd,  Grid<T>& bwd, 
					   Real strength, bool isLevelSet, bool isMAC=false )
{
	bool skip[3] = { false, false, false };

	if (!flags.isFluid(i,j,k)) skip[0] = skip[1] = skip[2] = true;
	if(isMAC) {
		if( (i>0) && (!flags.isFluid(i-1,j,k) )) skip[0] = true; 
		if( (j>0) && (!flags.isFluid(i,j-1,k) )) skip[1] = true; 
		if( (k>0) && (!flags.isFluid(i,j,k-1) )) skip[2] = true; 
	}

	for(int c=0; c<3; ++c ) {
		if ( skip[c] ) {
			dst(i,j,k)[c] = fwd(i,j,k)[c];
		} else { 
			// perform actual correction with given strength
			dst(i,j,k)[c] = fwd(i,j,k)[c] + strength * 0.5 * (old(i,j,k)[c] - bwd(i,j,k)[c]);
		}
	}
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
inline T doClampComponent(const Vec3i& gridSize, T dst, Grid<T>& orig, T fwd, const Vec3& pos, const Vec3& vel ) 
{
	T minv( std::numeric_limits<Real>::max()), maxv( -std::numeric_limits<Real>::max());

	// forward (and optionally) backward
	Vec3i positions[2];
	positions[0] = toVec3i(pos - vel);
	positions[1] = toVec3i(pos + vel);

	for(int l=0; l<2; ++l) {
		Vec3i& currPos = positions[l];

		// clamp forward lookup to grid
		const int i0 = clamp(currPos.x, 0, gridSize.x-1);
		const int j0 = clamp(currPos.y, 0, gridSize.y-1);
		const int k0 = clamp(currPos.z, 0, (orig.is3D() ? (gridSize.z-1) : 1) );
		const int i1 = i0+1, j1 = j0+1, k1= (orig.is3D() ? (k0+1) : k0);
		if( (!orig.isInBounds(Vec3i(i0,j0,k0),0)) || (!orig.isInBounds(Vec3i(i1,j1,k1),0)) )
			return fwd; 

		// find min/max around source pos
		getMinMax(minv, maxv, orig(i0,j0,k0));
		getMinMax(minv, maxv, orig(i1,j0,k0));
		getMinMax(minv, maxv, orig(i0,j1,k0));
		getMinMax(minv, maxv, orig(i1,j1,k0));

		if(orig.is3D()) {
		getMinMax(minv, maxv, orig(i0,j0,k1));
		getMinMax(minv, maxv, orig(i1,j0,k1));
		getMinMax(minv, maxv, orig(i0,j1,k1));
		getMinMax(minv, maxv, orig(i1,j1,k1)); } 
	}

	dst = clamp(dst, minv, maxv);
	return dst;
}
	
//! Helper function for clamping MAC grids
template<int c> 
inline Real doClampComponentMAC(const Vec3i& gridSize, Real dst, MACGrid& orig, Real fwd, const Vec3& pos, const Vec3& vel ) 
{
	Real minv = std::numeric_limits<Real>::max(), maxv = -std::numeric_limits<Real>::max();

	// forward (and optionally) backward
	Vec3i positions[2];
	positions[0] = toVec3i(pos - vel);
	positions[1] = toVec3i(pos + vel);

	for(int l=0; l<2; ++l) {
		Vec3i& currPos = positions[l];

		// clamp forward lookup to grid
		const int i0 = clamp(currPos.x, 0, gridSize.x-1);
		const int j0 = clamp(currPos.y, 0, gridSize.y-1);
		const int k0 = clamp(currPos.z, 0, (orig.is3D() ? (gridSize.z-1) : 1) );
		const int i1 = i0+1, j1 = j0+1, k1= (orig.is3D() ? (k0+1) : k0);
		if( (!orig.isInBounds(Vec3i(i0,j0,k0),0)) || (!orig.isInBounds(Vec3i(i1,j1,k1),0)) )
			return fwd; 

		// find min/max around source pos
		getMinMax(minv, maxv, orig(i0,j0,k0)[c]);
		getMinMax(minv, maxv, orig(i1,j0,k0)[c]);
		getMinMax(minv, maxv, orig(i0,j1,k0)[c]);
		getMinMax(minv, maxv, orig(i1,j1,k0)[c]);

		if(orig.is3D()) {
		getMinMax(minv, maxv, orig(i0,j0,k1)[c]);
		getMinMax(minv, maxv, orig(i1,j0,k1)[c]);
		getMinMax(minv, maxv, orig(i0,j1,k1)[c]);
		getMinMax(minv, maxv, orig(i1,j1,k1)[c]); } 
	}

	dst = clamp(dst, minv, maxv);
	return dst;
}

//! Kernel: Clamp obtained value to min/max in source area, and reset values that point out of grid or into boundaries
//          (note - MAC grids are handled below)
KERNEL(bnd=1) template<class T>
void MacCormackClamp(FlagGrid& flags, MACGrid& vel, Grid<T>& dst, Grid<T>& orig, Grid<T>& fwd, Real dt)
{
	T     dval       = dst(i,j,k);
	Vec3i gridUpper  = flags.getSize() - 1;
	
	dval = doClampComponent<T>(gridUpper, dval, orig, fwd(i,j,k), Vec3(i,j,k), vel.getCentered(i,j,k) * dt );

	// lookup forward/backward , round to closest NB
	Vec3i posFwd = toVec3i( Vec3(i,j,k) + Vec3(0.5,0.5,0.5) - vel.getCentered(i,j,k) * dt );
	Vec3i posBwd = toVec3i( Vec3(i,j,k) + Vec3(0.5,0.5,0.5) + vel.getCentered(i,j,k) * dt );
	
	// test if lookups point out of grid or into obstacle (note doClampComponent already checks sides, below is needed for valid flags access)
	if (posFwd.x < 0 || posFwd.y < 0 || posFwd.z < 0 ||
		posBwd.x < 0 || posBwd.y < 0 || posBwd.z < 0 ||
		posFwd.x > gridUpper.x || posFwd.y > gridUpper.y || ((posFwd.z > gridUpper.z)&&flags.is3D()) ||
		posBwd.x > gridUpper.x || posBwd.y > gridUpper.y || ((posBwd.z > gridUpper.z)&&flags.is3D()) ||
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
	Vec3  pos(i,j,k);
	Vec3  dval       = dst(i,j,k);
	Vec3  dfwd       = fwd(i,j,k);
	Vec3i gridUpper  = flags.getSize() - 1;
	
	dval.x = doClampComponentMAC<0>(gridUpper, dval.x, orig, dfwd.x, pos, vel.getAtMACX(i,j,k) * dt);
	dval.y = doClampComponentMAC<1>(gridUpper, dval.y, orig, dfwd.y, pos, vel.getAtMACY(i,j,k) * dt);
	dval.z = doClampComponentMAC<2>(gridUpper, dval.z, orig, dfwd.z, pos, vel.getAtMACZ(i,j,k) * dt);

	// note - the MAC version currently does not check whether source points were inside an obstacle! (unlike centered version)
	// this would need to be done for each face separately to stay symmetric...

	dst(i,j,k) = dval;
}


//! template function for performing SL advection
//! (Note boundary width only needed for specialization for MAC grids below)
template<class GridType> 
void fnAdvectSemiLagrange(FluidSolver* parent, FlagGrid& flags, MACGrid& vel, GridType& orig, int order, Real strength, int orderSpace, bool openBounds, int bWidth) {
	typedef typename GridType::BASETYPE T;
	
	Real dt = parent->getDt();
	bool levelset = orig.getType() & GridBase::TypeLevelset;
	
	// forward step
	GridType fwd(parent);
	SemiLagrange<T> (flags, vel, fwd, orig, dt, levelset, orderSpace);
	
	if (order == 1) {
		orig.swap(fwd);
	}
	else if (order == 2) { // MacCormack
		GridType bwd(parent);
		GridType newGrid(parent);
	
		// bwd <- backwards step
		SemiLagrange<T> (flags, vel, bwd, fwd, -dt, levelset, orderSpace);
		
		// newGrid <- compute correction
		MacCormackCorrect<T> (flags, newGrid, orig, fwd, bwd, strength, levelset);
		
		// clamp values
		MacCormackClamp<T> (flags, vel, newGrid, orig, fwd, dt);
		
		orig.swap(newGrid);
	}
}

// outflow functions

//! calculate local propagation velocity for cell (i,j,k)
Vec3 getBulkVel(FlagGrid& flags, MACGrid& vel, int i, int j, int k){
	Vec3 avg = Vec3(0,0,0);
	int count = 0;
	int size=1; // stencil size
	int nmax = (flags.is3D() ? size : 0);
	// average the neighboring fluid / outflow cell's velocity
	for (int n = -nmax; n<=nmax;n++){
		for (int m = -size; m<=size; m++){
			for (int l = -size; l<=size; l++){
				if (flags.isInBounds(Vec3i(i+l,j+m,k+n),0) && (flags.isFluid(i+l,j+m,k+n)||flags.isOutflow(i+l,j+m,k+n))){
					avg += vel(i+l,j+m,k+n);
					count++;
				}
			}
		}
	}
	return count>0 ? avg/count : avg;
}

//! extrapolate normal velocity components into outflow cell
KERNEL() void extrapolateVelConvectiveBC(FlagGrid& flags, MACGrid& vel, MACGrid& velDst, MACGrid& velPrev, Real timeStep, int bWidth){
	if (flags.isOutflow(i,j,k)){
		Vec3 bulkVel = getBulkVel(flags,vel,i,j,k);
		bool done=false;
		int dim = flags.is3D() ? 3 : 2;
		Vec3i cur, low, up, flLow, flUp;
		cur = low = up = flLow = flUp = Vec3i(i,j,k);
		// iterate over each velocity component x, y, z
		for (int c = 0; c<dim; c++){
			Real factor = timeStep*max((Real)1.0,bulkVel[c]); // prevent the extrapolated velocity from exploding when bulk velocity below 1
			low[c] = flLow[c] = cur[c]-1;
			up[c] = flUp[c] = cur[c]+1;
			// iterate over bWidth to allow for extrapolation into more distant outflow cells
			for (int d = 0; d<bWidth+1; d++){
				if (cur[c]>d && flags.isFluid(flLow)) {	
					velDst(i,j,k)[c] = ((vel(i,j,k)[c] - velPrev(i,j,k)[c]) / factor) + vel(low)[c];
					done=true;
				}
				if (cur[c]<flags.getSize()[c]-d-1 && flags.isFluid(flUp)) {
					// check for cells equally far away from two fluid cells -> average value between both sides
					if (done) velDst(i,j,k)[c] = 0.5*(velDst(i,j,k)[c] + ((vel(i,j,k)[c] - velPrev(i,j,k)[c]) / factor) + vel(up)[c]);
					else velDst(i,j,k)[c] = ((vel(i,j,k)[c] - velPrev(i,j,k)[c]) / factor) + vel(up)[c];
					done=true;
				}
				flLow[c]=flLow[c]-1;
				flUp[c]=flUp[c]+1;
				if (done) break;
			}
			low = up = flLow = flUp = cur;
			done=false;
		}
	}
}

//! copy extrapolated velocity components
KERNEL() void copyChangedVels(FlagGrid& flags, MACGrid& velDst, MACGrid& vel){ if (flags.isOutflow(i,j,k)) vel(i, j, k) = velDst(i, j, k); }

//! extrapolate normal velocity components into open boundary cells (marked as outflow cells)
void applyOutflowBC(FlagGrid& flags, MACGrid& vel, MACGrid& velPrev, double timeStep, int bWidth=1) {
	MACGrid velDst(vel.getParent()); // do not overwrite vel while it is read
	extrapolateVelConvectiveBC(flags, vel, velDst, velPrev, max(1.0,timeStep*4), bWidth);
	copyChangedVels(flags,velDst,vel);
}

// advection helpers

//! prevent parts of the surface getting "stuck" in obstacle regions
KERNEL() void knResetPhiInObs (FlagGrid& flags, Grid<Real>& sdf) {
	if( flags.isObstacle(i,j,k) && (sdf(i,j,k)<0.) ) {
		sdf(i,j,k) = 0.1;
	}
}
PYTHON() void resetPhiInObs (FlagGrid& flags, Grid<Real>& sdf) { knResetPhiInObs(flags, sdf); }

// advection main calls

//! template function for performing SL advection: specialized version for MAC grids
template<> 
void fnAdvectSemiLagrange<MACGrid>(FluidSolver* parent, FlagGrid& flags, MACGrid& vel, MACGrid& orig, int order, Real strength, int orderSpace, bool openBounds, int bWidth) {
	Real dt = parent->getDt();
	
	// forward step
	MACGrid fwd(parent);    
	SemiLagrangeMAC (flags, vel, fwd, orig, dt, orderSpace);
	
	if (orderSpace != 1) { debMsg("Warning higher order for MAC grids not yet implemented...",1); }

	if (order == 1) {
		if (openBounds) applyOutflowBC(flags, fwd, orig, dt, bWidth);
		orig.swap(fwd);
	}
	else if (order == 2) { // MacCormack 
		MACGrid bwd(parent);
		MACGrid newGrid(parent);
		
		// bwd <- backwards step
		SemiLagrangeMAC (flags, vel, bwd, fwd, -dt, orderSpace);
		
		// newGrid <- compute correction
		MacCormackCorrectMAC<Vec3> (flags, newGrid, orig, fwd, bwd, strength, false, true);
		
		// clamp values
		MacCormackClampMAC (flags, vel, newGrid, orig, fwd, dt); 
		
		if (openBounds) applyOutflowBC(flags, newGrid, orig, dt, bWidth);
		orig.swap(newGrid);
	}
}


//! Perform semi-lagrangian advection of target Real- or Vec3 grid
//! Open boundary handling needs information about width of border
PYTHON() void advectSemiLagrange (FlagGrid* flags, MACGrid* vel, GridBase* grid,
						   int order = 1, Real strength = 1.0, int orderSpace = 1, bool openBounds = false, int boundaryWidth = 1)
{    
	assertMsg(order==1 || order==2, "AdvectSemiLagrange: Only order 1 (regular SL) and 2 (MacCormack) supported");
	
	// determine type of grid    
	if (grid->getType() & GridBase::TypeReal) {
		fnAdvectSemiLagrange< Grid<Real> >(flags->getParent(), *flags, *vel, *((Grid<Real>*) grid), order, strength, orderSpace, openBounds, boundaryWidth);
	}
	else if (grid->getType() & GridBase::TypeMAC) {    
		fnAdvectSemiLagrange< MACGrid >(flags->getParent(), *flags, *vel, *((MACGrid*) grid), order, strength, orderSpace, openBounds, boundaryWidth);
	}
	else if (grid->getType() & GridBase::TypeVec3) {    
		fnAdvectSemiLagrange< Grid<Vec3> >(flags->getParent(), *flags, *vel, *((Grid<Vec3>*) grid), order, strength, orderSpace, openBounds, boundaryWidth);
	}
	else
		errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");    
}

} // end namespace DDF 

