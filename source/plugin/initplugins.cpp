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
#include "particle.h"
#include "noisefield.h"
#include "simpleimage.h"
#include "mesh.h"

using namespace std;

namespace Manta {
	
//! Apply noise to grid
KERNEL() 
void KnApplyNoiseInfl(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma) 
{
	if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	Real factor = clamp(1.0-0.5/sigma * (sdf(i,j,k)+sigma), 0.0, 1.0);
	
	Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
	if (density(i,j,k) < target)
		density(i,j,k) = target;
}

//! Init noise-modulated density inside shape
PYTHON() void densityInflow(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Shape* shape, Real scale=1.0, Real sigma=0)
{
	Grid<Real> sdf = shape->computeLevelset();
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
}
//! Apply noise to real grid based on an SDF
KERNEL() void KnAddNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf, Real scale) {
	if (!flags.isFluid(i,j,k) || (sdf && (*sdf)(i,j,k) > 0.) ) return;
	density(i,j,k) += noise.evaluate(Vec3(i,j,k)) * scale;
}
PYTHON() void addNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Grid<Real>* sdf=NULL, Real scale=1.0 ) {
	KnAddNoise(flags, density, noise, sdf, scale );
}

//! sample noise field and set pdata with its values (for convenience, scale the noise values)
KERNEL(pts) template<class T>
void knSetPdataNoise(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) {
	pdata[idx] = noise.evaluate( parts.getPos(idx) ) * scale;
}
KERNEL(pts) template<class T>
void knSetPdataNoiseVec(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) {
	pdata[idx] = noise.evaluateVec( parts.getPos(idx) ) * scale;
}
PYTHON() void setNoisePdata    (BasicParticleSystem& parts, ParticleDataImpl<Real>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<Real>(parts, pd,noise,scale); }
PYTHON() void setNoisePdataVec3(BasicParticleSystem& parts, ParticleDataImpl<Vec3>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoiseVec<Vec3>(parts, pd,noise,scale); }
PYTHON() void setNoisePdataInt (BasicParticleSystem& parts, ParticleDataImpl<int >& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<int> (parts, pd,noise,scale); }

//! SDF gradient from obstacle flags
PYTHON() Grid<Vec3> obstacleGradient(FlagGrid& flags) {
	LevelsetGrid levelset(flags.getParent(),false);
	Grid<Vec3> gradient(flags.getParent());
	
	// rebuild obstacle levelset
	FOR_IDX(levelset) {
		levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
	}
	levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);
	
	// build levelset gradient
	GradientOp(gradient, levelset);
	
	FOR_IDX(levelset) {
		Vec3 grad = gradient[idx];
		Real s = normalize(grad);
		if (s <= 0.1 || levelset[idx] >= 0) 
			grad=Vec3(0.);        
		gradient[idx] = grad * levelset[idx];
	}
	
	return gradient;
}

PYTHON() LevelsetGrid obstacleLevelset(FlagGrid& flags) {
   LevelsetGrid levelset(flags.getParent(),false);
	Grid<Vec3> gradient(flags.getParent());

	// rebuild obstacle levelset
	FOR_IDX(levelset) {
		levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
	}
	levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);

	return levelset;
}    


//*****************************************************************************
// blender init functions 

KERNEL() 
void KnApplyEmission(FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute) 
{
	if (!flags.isFluid(i,j,k) || emission(i,j,k) == 0.) return;
	if (isAbsolute)
		density(i,j,k) = emission(i,j,k);
	else
		density(i,j,k) += emission(i,j,k);
}

//! Add emission values
//isAbsolute: whether to add emission values to existing, or replace
PYTHON() void applyEmission(FlagGrid& flags, Grid<Real>& density, Grid<Real>& emission, bool isAbsolute) {
	KnApplyEmission(flags, density, emission, isAbsolute);
}

// blender init functions for meshes

KERNEL() 
void KnApplyDensity(FlagGrid& flags, Grid<Real>& density, Grid<Real>& sdf, Real value, Real sigma) 
{
	if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
	density(i,j,k) = value;
}
//! Init noise-modulated density inside mesh
PYTHON() void densityInflowMeshNoise(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Mesh* mesh, Real scale=1.0, Real sigma=0)
{
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 1.);
	KnApplyNoiseInfl(flags, density, noise, sdf, scale, sigma);
}

//! Init constant density inside mesh
PYTHON() void densityInflowMesh(FlagGrid& flags, Grid<Real>& density, Mesh* mesh, Real value=1., Real cutoff = 7, Real sigma=0)
{
	LevelsetGrid sdf(density.getParent(), false);
	mesh->computeLevelset(sdf, 2., cutoff);
	KnApplyDensity(flags, density, sdf, value, sigma);
}


//*****************************************************************************

//! check for symmetry , optionally enfore by copying
PYTHON() void checkSymmetry( Grid<Real>& a, Grid<Real>* err=NULL, bool symmetrize=false, int axis=0, int bound=0)
{
	const int c  = axis; 
	const int s = a.getSize()[c];
	FOR_IJK(a) { 
		Vec3i idx(i,j,k), mdx(i,j,k);
		mdx[c] = s-1-idx[c];
		if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

		if(err) (*err)(idx) = fabs( (double)(a(idx) - a(mdx) ) ); 
		if(symmetrize && (idx[c]<s/2)) {
			a(idx) = a(mdx);
		}
	}
}
//! check for symmetry , mac grid version
PYTHON() void checkSymmetryVec3( Grid<Vec3>& a, Grid<Real>* err=NULL, bool symmetrize=false , int axis=0, 
								int bound=0, int disable=0)
{
	if(err) err->setConst(0.);

	// each dimension is measured separately for flexibility (could be combined)
	const int c  = axis;
	const int o1 = (c+1)%3;
	const int o2 = (c+2)%3;

	// x
	if(! (disable&1) ) {
		const int s = a.getSize()[c]+1; 
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if(mdx[c] >= a.getSize()[c]) continue; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			// special case: center "line" of values , should be zero!
			if(mdx[c] == idx[c] ) {
				if(err) (*err)(idx) += fabs( (double)( a(idx)[c] ) ); 
				if(symmetrize) a(idx)[c] = 0.;
				continue; 
			}

			// note - the a(mdx) component needs to be inverted here!
			if(err) (*err)(idx) += fabs( (double)( a(idx)[c]- (a(mdx)[c]*-1.) ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[c] = a(mdx)[c] * -1.;
			}
		}
	}

	// y
	if(! (disable&2) ) {
		const int s = a.getSize()[c];
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			if(err) (*err)(idx) += fabs( (double)( a(idx)[o1]-a(mdx)[o1] ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[o1] = a(mdx)[o1];
			}
		}
	} 

	// z
	if(! (disable&4) ) {
		const int s = a.getSize()[c];
		FOR_IJK(a) { 
			Vec3i idx(i,j,k), mdx(i,j,k);
			mdx[c] = s-1-idx[c]; 
			if( bound>0 && ((!a.isInBounds(idx,bound)) || (!a.isInBounds(mdx,bound))) ) continue;

			if(err) (*err)(idx) += fabs( (double)( a(idx)[o2]-a(mdx)[o2] ) ); 
			if(symmetrize && (idx[c]<s/2)) {
				a(idx)[o2] = a(mdx)[o2];
			}
		}
	} 

}


// from simpleimage.cpp
void projectImg( SimpleImage& img, Grid<Real>& val, int shadeMode=0, Real scale=1.);

//! output shaded (all 3 axes at once for 3D)
//! shading modes: 0 smoke, 1 surfaces
PYTHON() void projectPpmFull( Grid<Real>& val, string name, int shadeMode=0, Real scale=1.)
{
	SimpleImage img;
	projectImg( img, val, shadeMode, scale );
	img.writePpm( name );
}

// helper functions for pdata operator tests

//! init some test particles at the origin
PYTHON() void addTestParts( BasicParticleSystem& parts, int num)
{
	for(int i=0; i<num; ++i)
		parts.addBuffered( Vec3(0,0,0) );

	parts.doCompress();
	parts.insertBufferedParticles();
}

// calculate the difference between two pdata fields (note - slow!, not parallelized)
PYTHON() Real pdataMaxDiff ( ParticleDataBase* a, ParticleDataBase* b )
{    
	double maxVal = 0.;
	//debMsg(" PD "<< a->getType()<<"  as"<<a->getSizeSlow()<<"  bs"<<b->getSizeSlow() , 1);
	assertMsg(a->getType()     == b->getType()    , "pdataMaxDiff problem - different pdata types!");
	assertMsg(a->getSizeSlow() == b->getSizeSlow(), "pdataMaxDiff problem - different pdata sizes!");
	
	if (a->getType() & ParticleDataBase::TypeReal) 
	{
		ParticleDataImpl<Real>& av = *dynamic_cast<ParticleDataImpl<Real>*>(a);
		ParticleDataImpl<Real>& bv = *dynamic_cast<ParticleDataImpl<Real>*>(b);
		FOR_PARTS(av) {
			maxVal = std::max(maxVal, (double)fabs( av[idx]-bv[idx] ));
		}
	} else if (a->getType() & ParticleDataBase::TypeInt) 
	{
		ParticleDataImpl<int>& av = *dynamic_cast<ParticleDataImpl<int>*>(a);
		ParticleDataImpl<int>& bv = *dynamic_cast<ParticleDataImpl<int>*>(b);
		FOR_PARTS(av) {
			maxVal = std::max(maxVal, (double)fabs( (double)av[idx]-bv[idx] ));
		}
	} else if (a->getType() & ParticleDataBase::TypeVec3) {
		ParticleDataImpl<Vec3>& av = *dynamic_cast<ParticleDataImpl<Vec3>*>(a);
		ParticleDataImpl<Vec3>& bv = *dynamic_cast<ParticleDataImpl<Vec3>*>(b);
		FOR_PARTS(av) {
			double d = 0.;
			for(int c=0; c<3; ++c) { 
				d += fabs( (double)av[idx][c] - (double)bv[idx][c] );
			}
			maxVal = std::max(maxVal, d );
		}
	} else {
		errMsg("pdataMaxDiff: Grid Type is not supported (only Real, Vec3, int)");    
	}

	return maxVal;
}

//*****************************************************************************
// helper functions for volume fractions (which are needed for second order obstacle boundaries)


KERNEL() void kninitVortexVelocity(Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius) {
	
	if(phiObs(i,j,k) >= -1.) {

		Real dx = i - center.x; if(dx>=0) dx -= .5; else dx += .5;
		Real dy = j - center.y;
		Real r = std::sqrt(dx*dx+dy*dy);
		Real alpha = atan2(dy,dx);

		vel(i,j,k).x = -std::sin(alpha)*(r/radius);

		dx = i - center.x;
		dy = j - center.y; if(dy>=0) dy -= .5; else dy += .5;
		r = std::sqrt(dx*dx+dy*dy);
		alpha = atan2(dy,dx);

		vel(i,j,k).y = std::cos(alpha)*(r/radius);

	}

}

PYTHON() void initVortexVelocity(Grid<Real> &phiObs, MACGrid& vel, const Vec3 &center, const Real &radius) {
	kninitVortexVelocity(phiObs,  vel, center, radius);
}

inline static Real calcFraction(Real phi1, Real phi2)
{
	if(phi1>0. && phi2>0.) return 1.;
	if(phi1<0. && phi2<0.) return 0.;

	// make sure phi1 < phi2
	if (phi2<phi1) { Real t = phi1; phi1= phi2; phi2 = t; }
	Real denom = phi1-phi2;
	if (denom > -1e-04) return 0.5; 

	Real frac = 1. - phi1/denom;
	if(frac<0.01) frac = 0.; // skip , dont mark as fluid
	return std::max(Real(0), std::min(Real(1), frac ));
}

KERNEL (bnd=1) 
void KnUpdateFractions(FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth) {

	// walls at domain bounds and inner objects
	fractions(i,j,k).x = calcFraction( phiObs(i,j,k) , phiObs(i-1,j,k));
	fractions(i,j,k).y = calcFraction( phiObs(i,j,k) , phiObs(i,j-1,k));
    if(phiObs.is3D()) {
	fractions(i,j,k).z = calcFraction( phiObs(i,j,k) , phiObs(i,j,k-1));
	}

	// remaining BCs at the domain boundaries 
	const int w = boundaryWidth;
	// only set if not in obstacle
 	if(phiObs(i,j,k)<0.) return;

	// x-direction boundaries
	if(i <= w+1) {                     //min x
		if( (flags.isInflow(i-1,j,k)) ||
			(flags.isOutflow(i-1,j,k)) ||
			(flags.isOpen(i-1,j,k)) ) {
				fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
	}
	if(i >= flags.getSizeX()-w-2) {    //max x
		if(	(flags.isInflow(i+1,j,k)) ||
			(flags.isOutflow(i+1,j,k)) ||
			(flags.isOpen(i+1,j,k)) ) {
			fractions(i+1,j,k).x = fractions(i+1,j,k).y = 1.; if(flags.is3D()) fractions(i+1,j,k).z = 1.;
		}
	}
	// y-direction boundaries
 	if(j <= w+1) {                     //min y
		if(	(flags.isInflow(i,j-1,k)) ||
			(flags.isOutflow(i,j-1,k)) ||
			(flags.isOpen(i,j-1,k)) ) {
			fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
 	}
 	if(j >= flags.getSizeY()-w-2) {      //max y
		if(	(flags.isInflow(i,j+1,k)) ||
			(flags.isOutflow(i,j+1,k)) ||
			(flags.isOpen(i,j+1,k)) ) {
			fractions(i,j+1,k).x = fractions(i,j+1,k).y = 1.; if(flags.is3D()) fractions(i,j+1,k).z = 1.;
		}
 	}
	// z-direction boundaries
	if(flags.is3D()) {
	if(k <= w+1) {                 //min z
		if(	(flags.isInflow(i,j,k-1)) ||
			(flags.isOutflow(i,j,k-1)) ||
			(flags.isOpen(i,j,k-1)) ) {
			fractions(i,j,k).x = fractions(i,j,k).y = 1.; if(flags.is3D()) fractions(i,j,k).z = 1.;
		}
	}
	if(j >= flags.getSizeZ()-w-2) { //max z
		if(	(flags.isInflow(i,j,k+1)) ||
			(flags.isOutflow(i,j,k+1)) ||
			(flags.isOpen(i,j,k+1)) ) {
			fractions(i,j,k+1).x = fractions(i,j,k+1).y = 1.; if(flags.is3D()) fractions(i,j,k+1).z = 1.;
		}
	}
	}

}

PYTHON() void updateFractions(FlagGrid& flags, Grid<Real>& phiObs, MACGrid& fractions, const int &boundaryWidth=0) {
	fractions.setConst( Vec3(0.) );
	KnUpdateFractions(flags, phiObs, fractions, boundaryWidth);
}

KERNEL (bnd=1) 
void KnUpdateFlags(FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs) {

	Real test = 0.;
	test += fractions.get(i  ,j,k).x;
	test += fractions.get(i+1,j,k).x;
	test += fractions.get(i,j  ,k).y;
	test += fractions.get(i,j+1,k).y;
	if (flags.is3D()) {
	test += fractions.get(i,j,k  ).z;
	test += fractions.get(i,j,k+1).z; }

	if(test==0. && phiObs(i,j,k) < 0.) flags(i,j,k) = FlagGrid::TypeObstacle; 
	else flags(i,j,k) = FlagGrid::TypeEmpty; 
}

PYTHON() void setObstacleFlags(FlagGrid& flags, MACGrid& fractions, Grid<Real>& phiObs) {
	KnUpdateFlags(flags,fractions, phiObs);
}

} // namespace

