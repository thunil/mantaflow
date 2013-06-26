/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Tools to setup fields and inflows
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "shapes.h"
#include "commonkernels.h"
#include "noisefield.h"

using namespace std;

namespace Manta {
    
//! Apply noise to grid
KERNEL 
void KnApplyNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma) 
{
    if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
    Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
    
    Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
    if (density(i,j,k) < target)
        density(i,j,k) = target;
}

//! Init noise-modulated density inside shape
PYTHON void densityInflow(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Shape* shape, Real scale=1.0, Real sigma=0)
{
    Grid<Real> sdf = shape->computeLevelset();    
    KnApplyNoise(flags, density, noise, sdf, scale, sigma);
}


    
//! Apply vector noise to grid
KERNEL 
void KnApplyNoiseVec(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma) 
{
    if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
    Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
    
    Vec3 noiseVec3 = noise.evaluateVec(Vec3(i,j,k)) * scale * factor;

	target(i,j,k) += noiseVec3;
    //if (dens(i,j,k) < target) dens(i,j,k) = target;
}

//! Apply vector-based wavelet noise to target grid
PYTHON void applyNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, Shape* shape, Real scale=1.0, Real sigma=0)
{
    Grid<Real> sdf = shape->computeLevelset();    
    KnApplyNoiseVec(flags, target, noise, sdf, scale, sigma);
}


//!interpolate grid from one size to another size
KERNEL 
void KnInterpolateGrid(Grid<Real>& source, Grid<Real>& target, const Vec3& sourceFactor)
{
	Vec3 pos = Vec3(i,j,k) * sourceFactor;
	target(i,j,k) = source.getInterpolated(pos);
}

PYTHON void interpolateGrid(Grid<Real>& source, Grid<Real>& target)
{
	Vec3 sourceFactor = Vec3( 
		Real(source.getSizeX())/target.getSizeX(), 
		Real(source.getSizeY())/target.getSizeY(), 
		Real(source.getSizeZ())/target.getSizeZ() );

	KnInterpolateGrid(source, target, sourceFactor);
}


//!interpolate a mac velocity grid from one size to another size
KERNEL 
void KnInterpolateMACGrid(MACGrid& target, MACGrid& source, const Vec3& sourceFactor)
{
	Vec3 pos = Vec3(i,j,k) * sourceFactor;

	Real vx = source.getInterpolated(pos - Vec3(0.5,0,0))[0];
	Real vy = source.getInterpolated(pos - Vec3(0,0.5,0))[1];
	Real vz = 0.f;
	if(source.is3D()) vz = source.getInterpolated(pos - Vec3(0,0,0.5))[2];

	target(i,j,k) = Vec3(vx,vy,vz);
}

PYTHON void interpolateMACGrid(MACGrid& target, MACGrid& source)
{
	Vec3 sourceFactor = Vec3( 
		Real(source.getSizeX())/target.getSizeX(), 
		Real(source.getSizeY())/target.getSizeY(), 
		Real(source.getSizeZ())/target.getSizeZ() );

	KnInterpolateMACGrid(target, source, sourceFactor);
}


    
} // namespace
