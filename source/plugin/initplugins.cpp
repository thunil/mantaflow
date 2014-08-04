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

//! sample noise field and set pdata with its values (for convenience, scale the noise values)
KERNEL(pts) template<class T>
void knSetPdataNoise(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) {
	pdata[idx] = noise.evaluate( parts.getPos(idx) ) * scale;
}
KERNEL(pts) template<class T>
void knSetPdataNoiseVec(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) {
	pdata[idx] = noise.evaluateVec( parts.getPos(idx) ) * scale;
}
PYTHON void setNoisePdata    (BasicParticleSystem& parts, ParticleDataImpl<Real>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<Real>(parts, pd,noise,scale); }
PYTHON void setNoisePdataVec3(BasicParticleSystem& parts, ParticleDataImpl<Vec3>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoiseVec<Vec3>(parts, pd,noise,scale); }
PYTHON void setNoisePdataInt (BasicParticleSystem& parts, ParticleDataImpl<int >& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<int> (parts, pd,noise,scale); }

//! SDF gradient from obstacle flags
PYTHON Grid<Vec3> obstacleGradient(FlagGrid& flags) {
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

PYTHON LevelsetGrid obstacleLevelset(FlagGrid& flags) {
   LevelsetGrid levelset(flags.getParent(),false);
	Grid<Vec3> gradient(flags.getParent());

	// rebuild obstacle levelset
	FOR_IDX(levelset) {
		levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
	}
	levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);

	return levelset;
}    

//! check for symmetry , optionally enfore by copying
PYTHON void checkSymmetry( Grid<Real>& a, Grid<Real>* err=NULL, bool symmetrize=false, int axis=0, int bound=0)
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
PYTHON void checkSymmetryVec3( Grid<Vec3>& a, Grid<Real>* err=NULL, bool symmetrize=false , int axis=0, 
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


//! project data onto a plane, and write ppm
PYTHON void projectPpmOut( Grid<Real>& val, string name, int axis=2, Real scale=1.)
{
	const int c  = axis; 
	const int o1 = (c+1)%3;
	const int o2 = (c+2)%3;

	SimpleImage img;
	img.init( val.getSize()[o1],  val.getSize()[o2] );
	Real s = 1. / (Real)val.getSize()[c];

	FOR_IJK(val) { 
		Vec3i idx(i,j,k); 
		//if(idx[c]==val.getSize()[c]/2) img( idx[o1], idx[o2] ) = val(idx); 
		img( idx[o1], idx[o2] ) += s * val(idx);
	}
	img.mapRange( 1./scale );
	img.writePpm( name );
}

KERNEL 
void KnPrecompLight(Grid<Real>& density, Grid<Real>& L, Vec3 light = Vec3(1,1,1) , Grid<Vec3>* deb=NULL)
{
	Vec3 n = getGradient( density, i,j,k ) * -1.; 
	normalize(n);

	Real d = dot( light, n );
	L(i,j,k) = d;
	if(deb) (*deb)(i,j,k) = Vec3(d,d,d);
}

// simple shading with pre-computed gradient
static inline void shadeCell(Vec3& dst, int shadeMode, Real src, Real light, int depthPos, Real depthInv) 
{	
	switch(shadeMode) {

	case 1: {
		// surfaces
		Vec3 ambient = Vec3(0.1,0.1,0.1);
		Vec3 diffuse = Vec3(0.9,0.9,0.9); 
		Real alpha = src; 

		// different color for depth?
		diffuse[0] *= ((Real)depthPos * depthInv) * 0.7 + 0.3;
		diffuse[1] *= ((Real)depthPos * depthInv) * 0.7 + 0.3;

		Vec3 col = ambient + diffuse * light; 

		//img( 0+i, j ) = (1.-alpha) * img( 0+i, j ) + alpha * col;
		dst = (1.-alpha) * dst + alpha * col;
		} break;

	default: {
		// volumetrics / smoke
		dst += depthInv * Vec3(src,src,src);
		} break;

	}
}

//! output shaded (all 3 axes at once for 3D)
//! shading modes: 0 smoke, 1 surfaces
PYTHON void projectPpmFull( Grid<Real>& val, string name, int shadeMode=0, Real scale=1.,
		Grid<Vec3>* deb=NULL) // NT_DEBUG
{
	Vec3i s  = val.getSize();
	Vec3  si = Vec3( 1. / (Real)s[0], 1. / (Real)s[1], 1. / (Real)s[2] );

	SimpleImage img;
	int imgSx = s[0];
	if(val.is3D()) imgSx += s[2]+s[0]; // mult views in 3D
	img.init( imgSx, std::max(s[0], std::max( s[1],s[2])) );

	// precompute lighting
	Grid<Real> L(val);
	KnPrecompLight( val, L , Vec3(1,1,1), deb);

	FOR_IJK(val) { 
		Vec3i idx(i,j,k);
		// img( 0+i, j ) += si[2] * val(idx); // averaging
		shadeCell( img( 0+i, j ) , shadeMode, val(idx), L(idx), k, si[2]);
	}

	if( val.is3D() ) {

	FOR_IJK(val) { 
		Vec3i idx(i,j,k);
		//img( s[0]+k, j ) += si[0] * val(idx);
		shadeCell( img( s[0]+k, j ) , shadeMode, val(idx), L(idx), i, si[0]);
	}

	FOR_IJK(val) { 
		Vec3i idx(i,j,k);
		//img( s[0]+s[2]+i, k ) += si[1] * val(idx);
		shadeCell( img( s[0]+s[2]+i, k ) , shadeMode, val(idx), L(idx), j, si[1]);
	}

	} // 3d

	img.mapRange( 1./scale );
	img.writePpm( name );
}

// helper functions for pdata operator tests

//! init some test particles at the origin
PYTHON void addTestParts( BasicParticleSystem& parts, int num)
{
	for(int i=0; i<num; ++i)
		parts.addBuffered( Vec3(0,0,0) );

	parts.doCompress();
	parts.insertBufferedParticles();
}

// calculate the difference between two pdata fields (note - slow!, not parallelized)
PYTHON Real pdataMaxDiff ( ParticleDataBase* a, ParticleDataBase* b )
{    
	double maxVal = 0.;
	//debMsg(" PD "<< a->getType()<<"  as"<<a->getSizeSlow()<<"  bs"<<b->getSizeSlow() , 1);
	assertMsg(a->getType()     == b->getType()    , "pdataMaxDiff problem - different pdata types!");
	assertMsg(a->getSizeSlow() == b->getSizeSlow(), "pdataMaxDiff  problem -different pdata sizes!");
	
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

} // namespace

