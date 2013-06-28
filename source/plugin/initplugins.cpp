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
void KnApplyNoiseVec(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, 
					  Real scale, Grid<Real>* weight ) 
{
	if ( !flags.isFluid(i,j,k) ) return;

	Real factor = 1;
	if(weight) factor = (*weight)(i,j,k);

	Vec3 noiseVec3 = noise.evaluateVec( Vec3(i,j,k) ) * scale * factor;

	target(i,j,k) += noiseVec3;
}

//! Apply vector-based wavelet noise to target grid
PYTHON void applyNoiseVec3(FlagGrid& flags, Grid<Vec3>& target, WaveletNoiseField& noise, 
							Real scale=1.0 , Grid<Real>* weight=NULL )
{
    // note - passing a MAC grid here is slightly inaccurate, we should
    // evaluate each component separately
	KnApplyNoiseVec(flags, target, noise, scale , weight );
}



//! Compute energy of a staggered velocity field (at cell center)
KERNEL 
void KnApplyComputeEnergy( FlagGrid& flags, MACGrid& vel, Grid<Real>& energy ) 
{
    Real e = 0.f;
    if ( flags.isFluid(i,j,k) ) {
        Vec3 v = vel.getCentered(i,j,k);
        e = 0.5 * v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    }
    energy(i,j,k) = e;
}

PYTHON void computeEnergy( FlagGrid& flags, MACGrid& vel, Grid<Real>& energy )
{
    KnApplyComputeEnergy( flags, vel, energy );
}



//!interpolate grid from one size to another size
KERNEL 
void KnInterpolateGrid(Grid<Real>& target, Grid<Real>& source, const Vec3& sourceFactor)
{
	Vec3 pos = Vec3(i,j,k) * sourceFactor;
    if(!source.is3D()) pos[2] = 0;
	target(i,j,k) = source.getInterpolated(pos);
}

PYTHON void interpolateGrid( Grid<Real>& target, Grid<Real>& source )
{
	Vec3 sourceFactor = Vec3( 
		Real(source.getSizeX())/target.getSizeX(), 
		Real(source.getSizeY())/target.getSizeY(), 
		Real(source.getSizeZ())/target.getSizeZ() );

	// a brief note on a mantaflow specialty: the target grid has to be the first argument here!
	// the parent fluidsolver object is taken from the first grid, and it determines the size of the
	// loop for the kernel call. as we're writing into target, it's important to loop exactly over
	// all cells of the target grid... (note, when calling the plugin in python, it doesnt matter anymore).

	KnInterpolateGrid(target, source, sourceFactor);
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

	// see interpolateGrid for why the target grid needs to come first in the parameters!

	KnInterpolateMACGrid(target, source, sourceFactor);
}

PYTHON void computeWaveletCoeffs(Grid<Real>& input)
{
    Grid<Real> temp1(input.getParent()), temp2(input.getParent());
    WaveletNoiseField::computeCoefficients(input, temp1, temp2);
}

} // namespace
