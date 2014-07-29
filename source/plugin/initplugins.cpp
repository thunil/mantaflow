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
PYTHON void checkSymmetry( Grid<Real>& a, Grid<Real>* err=NULL, bool symmetrize=false )
{
	// only along x for now!
	const int s = a.getSize()[0];
	FOR_IJK(a) { 
		if(err) (*err)(i,j,k) = fabs( (double)(a(i,j,k) - a(s-1-i,j,k) ) ); 
		if(symmetrize && (i<s/2)) {
			a(i,j,k) = a(s-1-i,j,k);
		}
	}
}
PYTHON void checkSymmetryVec3( Grid<Vec3>& a, Grid<Real>* err=NULL, bool symmetrize=false, int disable=0 )
{
	if(err) err->setConst(0.);

	//FOR_IJK(a) { a(i,j,k)[0] = 0.; }

	// only along x for now!
	// each dimension is measured separately for flexibility (could be combined)

	// x
	if(! (disable&1) ) {
		const int s = a.getSize()[0]+1; 
		FOR_IJK(a) { 
			if(s-1-i >= a.getSize()[0]) continue; 
			if(s-1-i == i             ) continue; 
			Vec3 v1 = a(i    ,j,k);
			Vec3 v2 = a(s-1-i,j,k);
			v1[1] = v2[1] = 0.;
			v1[2] = v2[2] = 0.;
			v1[0] *=  1.;
			v2[0] *= -1.;
			if(err) (*err)(i,j,k) += fabs( (double)( norm(v1 - v2) ) ); 
			if(symmetrize && (i<s/2)) {
				a(i,j,k)[0] = -a(s-1-i,j,k)[0] + 0.;
			}
		}
	}

	// y
	if(! (disable&2) ) {
		const int s = a.getSize()[0];
		FOR_IJK(a) { 
			Vec3 v1=a(i    ,j,k);
			Vec3 v2=a(s-1-i,j,k);
			v1[0] = v2[0] = 0.;
			v1[2] = v2[2] = 0.;
			v1[1] *=  1.;
			v2[1] *=  1.;
			if(err) (*err)(i,j,k) += fabs( (double)( norm(v1 - v2) ) ); 
			if(symmetrize && (i<s/2)) {
				a(i,j,k)[1] = a(s-1-i,j,k)[1];
			}
		}
	} 

	// z
	if(! (disable&4) ) {
		const int s = a.getSize()[0];
		FOR_IJK(a) { 
			Vec3 v1=a(i    ,j,k);
			Vec3 v2=a(s-1-i,j,k);
			v1[0] = v2[0] = 0.;
			v1[1] = v2[1] = 0.;
			v1[2] *=  1.;
			v2[2] *=  1.;
			if(err) (*err)(i,j,k) += fabs( (double)( norm(v1 - v2) ) ); 
			if(symmetrize && (i<s/2)) {
				a(i,j,k)[2] = a(s-1-i,j,k)[2];
			}
		}
	} 

}

// only along x for now!
/*
PYTHON void symmetrizeReal( Grid<Real>& a )
{
	const int s = a.getSizeX();
	FOR_IJK(a) { 
		if(i>s/2) continue;
		a(i,j,k) = a(s-1-i,j,k); 
	}
}
PYTHON void symmetrizeVec3( Grid<Vec3>& a )
{
	// only along x for now!
	if(1) {
		const int s = a.getSizeX();
		FOR_IJK(a) { 
			if(i>s/2) continue;
			Vec3& v1 = a(i    , j, k);
			Vec3& v2 = a(s-1-i, j, k);
			v2[1] = v1[1];
			v2[2] = v1[2];
		}
	} 
	if(1){
		const int s = a.getSizeX()+1; 
		FOR_IJK(a) { 
			if(i>s/2) continue;
			//if(i>s-1) { err(i,j,k) = 0.; }
			//else      
			{ 
				Vec3& v1 = a(i    ,j,k);
				Vec3& v2 = a(s-1-i,j,k);
				v2[0] = -v1[0];
			}
		}
	}
}
*/


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

