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
	pdata[i] = noise.evaluate( parts.getPos(i) ) * scale;
}
KERNEL(pts) template<class T>
void knSetPdataNoiseVec(BasicParticleSystem& parts, ParticleDataImpl<T>& pdata, WaveletNoiseField& noise, Real scale) {
	pdata[i] = noise.evaluateVec( parts.getPos(i) ) * scale;
}
PYTHON void setNoisePdata    (BasicParticleSystem& parts, ParticleDataImpl<Real>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<Real>(parts, pd,noise,scale); }
PYTHON void setNoisePdataVec3(BasicParticleSystem& parts, ParticleDataImpl<Vec3>& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoiseVec<Vec3>(parts, pd,noise,scale); }
PYTHON void setNoisePdataInt (BasicParticleSystem& parts, ParticleDataImpl<int >& pd, WaveletNoiseField& noise, Real scale=1.) { knSetPdataNoise<int> (parts, pd,noise,scale); }

//! SDF gradient from obstacle flags
PYTHON Grid<Vec3> obstacleGradient(FlagGrid& flags) {
    LevelsetGrid levelset(parent,false);
    Grid<Vec3> gradient(parent);
    
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
   LevelsetGrid levelset(parent,false);
    Grid<Vec3> gradient(parent);

    // rebuild obstacle levelset
    FOR_IDX(levelset) {
        levelset[idx] = flags.isObstacle(idx) ? -0.5 : 0.5;
    }
    levelset.reinitMarching(flags, 6.0, 0, true, false, FlagGrid::TypeReserved);

    return levelset;
}    


} // namespace
