/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugins for pressure correction: solve_pressure, and ghost fluid helpers
 *
 ******************************************************************************/
#include "vectorbase.h"
#include "kernel.h"
#include "conjugategrad.h"
#include "multigrid.h"

using namespace std;
namespace Manta {

enum Preconditioner { PcNone = 0, PcMIC = 1, PcMG = 2 };

//! Kernel: Construct the right-hand side of the poisson equation
KERNEL(bnd=1, reduce=+) returns(int cnt=0) returns(double sum=0)
void MakeRhs (FlagGrid& flags, Grid<Real>& rhs, MACGrid& vel, 
			  Grid<Real>* perCellCorr, MACGrid* fractions)
{
	if (!flags.isFluid(i,j,k)) {
		rhs(i,j,k) = 0;
		return;
	}

	// compute divergence 
	// no flag checks: assumes vel at obstacle interfaces is set to zero
	Real set(0);
	if(!fractions) {
		set =               vel(i,j,k).x - vel(i+1,j,k).x + 
				 			vel(i,j,k).y - vel(i,j+1,k).y; 
		if(vel.is3D()) set+=vel(i,j,k).z - vel(i,j,k+1).z;
	}else{
		set =               (*fractions)(i,j,k).x * vel(i,j,k).x - (*fractions)(i+1,j,k).x * vel(i+1,j,k).x + 
							(*fractions)(i,j,k).y * vel(i,j,k).y - (*fractions)(i,j+1,k).y * vel(i,j+1,k).y; 
		if(vel.is3D()) set+=(*fractions)(i,j,k).z * vel(i,j,k).z - (*fractions)(i,j,k+1).z * vel(i,j,k+1).z;
	}
	
	// per cell divergence correction (optional)
	if(perCellCorr) 
		set += perCellCorr->get(i,j,k);
	
	// obtain sum, cell count
	sum += set;
	cnt++;
	
	rhs(i,j,k) = set;
}

//! Kernel: Apply velocity update from poisson equation
KERNEL(bnd = 1)
void CorrectVelocity(FlagGrid& flags, MACGrid& vel, Grid<Real>& pressure) 
{
	IndexInt idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if (flags.isFluid(i-1,j,k)) vel[idx].x -= (pressure[idx] - pressure(i-1,j,k));
		if (flags.isFluid(i,j-1,k)) vel[idx].y -= (pressure[idx] - pressure(i,j-1,k));
		if (flags.is3D() && flags.isFluid(i,j,k-1)) vel[idx].z -= (pressure[idx] - pressure(i,j,k-1));
 
		if (flags.isEmpty(i-1,j,k)) vel[idx].x -= pressure[idx];
		if (flags.isEmpty(i,j-1,k)) vel[idx].y -= pressure[idx];
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx].z -= pressure[idx];
	}
	else if (flags.isEmpty(idx)&&!flags.isOutflow(idx)) // don't change velocities in outflow cells
	{
		if (flags.isFluid(i-1,j,k)) vel[idx].x += pressure(i-1,j,k);
		else                        vel[idx].x  = 0.f;
		if (flags.isFluid(i,j-1,k)) vel[idx].y += pressure(i,j-1,k);
		else                        vel[idx].y  = 0.f;
		if (flags.is3D() ) {
		if (flags.isFluid(i,j,k-1)) vel[idx].z += pressure(i,j,k-1);
		else                        vel[idx].z  = 0.f;
		}
	}
}

// *****************************************************************************
// Ghost fluid helpers

// calculate fraction filled with liquid (note, assumes inside value is < outside!)
inline static Real thetaHelper(Real inside, Real outside)
{
	Real denom = inside-outside;
	if (denom > -1e-04) return 0.5; // should always be neg, and large enough...
	return std::max(Real(0), std::min(Real(1), inside/denom));
}

// calculate ghost fluid factor, cell at idx should be a fluid cell
inline static Real ghostFluidHelper(IndexInt idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return alpha = gfClamp;
	return (1-(1/alpha)); 
}

//! Kernel: Adapt A0 for ghost fluid
KERNEL(bnd=1) 
void ApplyGhostFluidDiagonal(Grid<Real> &A0, const FlagGrid &flags, const Grid<Real> &phi, Real gfClamp)
{
	const int X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	IndexInt idx = flags.index(i,j,k);
	if (!flags.isFluid(idx)) return;

	if (flags.isEmpty(i-1,j,k)) A0[idx] -= ghostFluidHelper(idx, -X, phi, gfClamp);
	if (flags.isEmpty(i+1,j,k)) A0[idx] -= ghostFluidHelper(idx, +X, phi, gfClamp);
	if (flags.isEmpty(i,j-1,k)) A0[idx] -= ghostFluidHelper(idx, -Y, phi, gfClamp);
	if (flags.isEmpty(i,j+1,k)) A0[idx] -= ghostFluidHelper(idx, +Y, phi, gfClamp);
	if (flags.is3D()) {
		if (flags.isEmpty(i,j,k-1)) A0[idx] -= ghostFluidHelper(idx, -Z, phi, gfClamp);
		if (flags.isEmpty(i,j,k+1)) A0[idx] -= ghostFluidHelper(idx, +Z, phi, gfClamp);
	}
}

//! Kernel: Apply velocity update: ghost fluid contribution
KERNEL(bnd=1)
void CorrectVelocityGhostFluid(MACGrid &vel, const FlagGrid &flags, const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp)
{
	const IndexInt X = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	const IndexInt idx = flags.index(i,j,k);
	if (flags.isFluid(idx))
	{
		if (flags.isEmpty(i-1,j,k)) vel[idx][0] += pressure[idx] * ghostFluidHelper(idx, -X, phi, gfClamp);
		if (flags.isEmpty(i,j-1,k)) vel[idx][1] += pressure[idx] * ghostFluidHelper(idx, -Y, phi, gfClamp);
		if (flags.is3D() && flags.isEmpty(i,j,k-1)) vel[idx][2] += pressure[idx] * ghostFluidHelper(idx, -Z, phi, gfClamp);
	}
	else if (flags.isEmpty(idx)&&!flags.isOutflow(idx)) // do not change velocities in outflow cells
	{
		if (flags.isFluid(i-1,j,k)) vel[idx][0] -= pressure(i-1,j,k) * ghostFluidHelper(idx-X, +X, phi, gfClamp);
		else                        vel[idx].x  = 0.f;
		if (flags.isFluid(i,j-1,k)) vel[idx][1] -= pressure(i,j-1,k) * ghostFluidHelper(idx-Y, +Y, phi, gfClamp);
		else                        vel[idx].y  = 0.f;
		if (flags.is3D() ) {
		if (flags.isFluid(i,j,k-1)) vel[idx][2] -= pressure(i,j,k-1) * ghostFluidHelper(idx-Z, +Z, phi, gfClamp);
		else                        vel[idx].z  = 0.f;
		}
	}
}


// improve behavior of clamping for large time steps:

inline static Real ghostFluidWasClamped(IndexInt idx, int offset, const Grid<Real> &phi, Real gfClamp)
{
	Real alpha = thetaHelper(phi[idx], phi[idx+offset]);
	if (alpha < gfClamp) return true;
	return false;
}

KERNEL(bnd=1)
void ReplaceClampedGhostFluidVels(MACGrid &vel, FlagGrid &flags, 
		const Grid<Real> &pressure, const Grid<Real> &phi, Real gfClamp )
{
	const IndexInt idx = flags.index(i,j,k);
	const IndexInt X   = flags.getStrideX(), Y = flags.getStrideY(), Z = flags.getStrideZ();
	if (!flags.isEmpty(idx)) return;

	if( (flags.isFluid(i-1,j,k)) && ( ghostFluidWasClamped(idx-X, +X, phi, gfClamp) ) )
		vel[idx][0] = vel[idx-X][0];
	if( (flags.isFluid(i,j-1,k)) && ( ghostFluidWasClamped(idx-Y, +Y, phi, gfClamp) ) )
		vel[idx][1] = vel[idx-Y][1];
	if( flags.is3D() &&
	  ( (flags.isFluid(i,j,k-1)) && ( ghostFluidWasClamped(idx-Z, +Z, phi, gfClamp) ) ))
		vel[idx][2] = vel[idx-Z][2];

	if( (flags.isFluid(i+1,j,k)) && ( ghostFluidWasClamped(idx+X, -X, phi, gfClamp)) )
		vel[idx][0] = vel[idx+X][0];
	if( (flags.isFluid(i,j+1,k)) && ( ghostFluidWasClamped(idx+Y, -Y, phi, gfClamp)) )
		vel[idx][1] = vel[idx+Y][1];
	if( flags.is3D() && 
	   (flags.isFluid(i,j,k+1))  && ( ghostFluidWasClamped(idx+Z, -Z, phi, gfClamp)) )
		vel[idx][2] = vel[idx+Z][2];
}

//! Kernel: Compute min value of Real grid
KERNEL(idx, reduce=+) returns(int numEmpty=0)
int CountEmptyCells(FlagGrid& flags) {
	if (flags.isEmpty(idx) ) numEmpty++;
}

// *****************************************************************************
// Main pressure solve

IndexInt solvePressureBase(MACGrid& vel, Grid<Real>& pressure, FlagGrid& flags, Grid<Real>& rhs, Real cgAccuracy = 1e-3,
	Grid<Real>* phi = 0,
	Grid<Real>* perCellCorr = 0,
	MACGrid* fractions = 0,
	Real gfClamp = 1e-04,
	Real cgMaxIterFac = 1.5,
	int preconditioner = PcMIC,
	bool enforceCompatibility = false,
	bool useL2Norm = false,
	bool zeroPressureFixing = false,
	Grid<Real>* retRhs = NULL)
{
	// reserve temp grids
	FluidSolver* parent = flags.getParent();
	Grid<Real> residual(parent);
	Grid<Real> search(parent);
	Grid<Real> A0(parent);
	Grid<Real> Ai(parent);
	Grid<Real> Aj(parent);
	Grid<Real> Ak(parent);
	Grid<Real> tmp(parent);
	Grid<Real> pca0(parent);
	Grid<Real> pca1(parent);
	Grid<Real> pca2(parent);
	Grid<Real> pca3(parent);
		
	// setup matrix and boundaries 
	MakeLaplaceMatrix (flags, A0, Ai, Aj, Ak, fractions);

	if (phi) {
		ApplyGhostFluidDiagonal(A0, flags, *phi, gfClamp);
	}
	
	// compute divergence and init right hand side
	MakeRhs kernMakeRhs (flags, rhs, vel, perCellCorr, fractions);
	
	if (enforceCompatibility)
		rhs += (Real)(-kernMakeRhs.sum / (Real)kernMakeRhs.cnt);
	
	// check whether we need to fix some pressure value...
	// (manually enable, or automatically for high accuracy, can cause asymmetries otherwise)
	IndexInt fixPidx = -1;
	if(zeroPressureFixing || cgAccuracy<1e-07) 
	{
		if(FLOATINGPOINT_PRECISION==1) debMsg("Warning - high CG accuracy with single-precision floating point accuracy might not converge...", 2);

		int numEmpty = CountEmptyCells(flags);
		if(numEmpty==0) {
			FOR_IJK_BND(flags,1) {
				if(flags.isFluid(i,j,k)) {
					fixPidx = flags.index(i,j,k);
					break;
				}
			}
			//debMsg("No empty cells! Fixing pressure of cell "<<fixPidx<<" to zero",1);
		}
		if(fixPidx>=0) {
			// adjustment for approx. symmetric zeroPressureFixing cell (top center)
			fixPidx = flags.index(flags.getSizeX() / 2, flags.getSizeY() - 2, flags.is3D() ? flags.getSizeZ() / 2 : 0);

			flags[fixPidx] |= FlagGrid::TypeZeroPressure;
			rhs[fixPidx] = 0.; 
			debMsg("Pinning pressure of cell "<<fixPidx<<" to zero", 2);
		}
	}


	if (preconditioner == PcNone || preconditioner == PcMIC)
	{
		// CG setup
		// note: the last factor increases the max iterations for 2d, which right now can't use a preconditioner 
		const int maxIter = (int)(cgMaxIterFac * flags.getSize().max()) * (flags.is3D() ? 1 : 4);
		GridCgInterface *gcg;
		if (vel.is3D())
			gcg = new GridCg<ApplyMatrix>  (pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
		else
			gcg = new GridCg<ApplyMatrix2D>(pressure, rhs, residual, search, flags, tmp, &A0, &Ai, &Aj, &Ak );
	
		gcg->setAccuracy( cgAccuracy ); 
		gcg->setUseL2Norm( useL2Norm );

		// optional preconditioning
		gcg->setPreconditioner( preconditioner ? GridCgInterface::PC_mICP : GridCgInterface::PC_None, &pca0, &pca1, &pca2, &pca3);

		for (int iter=0; iter<maxIter; iter++) {
			if (!gcg->iterate()) iter=maxIter;
		} 
		//debMsg("FluidSolver::solvePressureBase iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
		delete gcg;
	}
	else if (preconditioner == PcMG)
	{
		const int maxIter = 10;
		GridMg* gmg = new GridMg(pressure);

		gmg->setA(flags, &A0, &Ai, &Aj, &Ak);
		gmg->setRhs(rhs);

		gmg->setAccuracy(cgAccuracy);

		for (int iter = 0; iter<maxIter; iter++) {
			if (!gmg->doVCycle(pressure)) iter=maxIter;
		}
		//debMsg("FluidSolver::solvePressureBase iterations:"<<gcg->getIterations()<<", res:"<<gcg->getSigma(), 1);
		delete gmg;
	}
	
	return fixPidx;
}

//! Perform pressure projection of the velocity grid
PYTHON() void solvePressure(MACGrid& vel, Grid<Real>& pressure, FlagGrid& flags, Real cgAccuracy = 1e-3,
    Grid<Real>* phi = 0, 
    Grid<Real>* perCellCorr = 0, 
    MACGrid* fractions = 0,
    Real gfClamp = 1e-04,
    Real cgMaxIterFac = 1.5,
	int preconditioner = PcMIC,
	bool enforceCompatibility = false,
    bool useL2Norm = false, 
	bool zeroPressureFixing = false,
	Grid<Real>* retRhs = NULL )
{
	Grid<Real> rhs(vel.getParent()); 
	IndexInt fixPidx = solvePressureBase(vel, pressure, flags, rhs, cgAccuracy, phi, perCellCorr, fractions, gfClamp, cgMaxIterFac, preconditioner, enforceCompatibility, useL2Norm, zeroPressureFixing, retRhs);
	
	CorrectVelocity(flags, vel, pressure ); 
	if (phi) {
		CorrectVelocityGhostFluid (vel, flags, pressure, *phi, gfClamp);
		// improve behavior of clamping for large time steps:
		ReplaceClampedGhostFluidVels (vel, flags, pressure, *phi, gfClamp);
	}

	if(fixPidx>=0)
		flags[fixPidx] &= ~FlagGrid::TypeZeroPressure;

	// optionally , return RHS
	if(retRhs) {
		retRhs->copyFrom( rhs );
	}
}

} // end namespace

