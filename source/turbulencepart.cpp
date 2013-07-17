/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Turbulence particles
 *
 ******************************************************************************/

#include "turbulencepart.h"
#include "shapes.h"
#include "randomstream.h"

using namespace std;
namespace Manta {
    
TurbulenceParticleSystem::TurbulenceParticleSystem(FluidSolver* parent, WaveletNoiseField& noise) :
    ParticleSystem<TurbulenceParticleData>(parent), noise(noise)
{ 
}

ParticleBase* TurbulenceParticleSystem::clone() {
    TurbulenceParticleSystem* nm = new TurbulenceParticleSystem(getParent(), noise);
    compress();
    
    nm->mData = mData;
    nm->setName(getName());
    return nm;
}

void TurbulenceParticleSystem::seed(Shape* shape, int num) {
    static RandomStream rand(34894231);
    Vec3 sz = shape->getExtent(), p0 = shape->getCenter() - sz*0.5;
    for (int i=0; i<num; i++) {
        Vec3 p;
        do {
            p = rand.getVec3() * sz + p0;            
        } while(!shape->isInside(p));
        add(TurbulenceParticleData(p));        
    }
}

void TurbulenceParticleSystem::resetTexCoords(int num, const Vec3& inflow ) {
    if (num==0) {
        for (int i=0; i<size(); i++) mData[i].tex0 = mData[i].pos - inflow;
    } else {
        for (int i=0; i<size(); i++) mData[i].tex1 = mData[i].pos - inflow;
    }
}


KERNEL(pts)
void KnSynthesizeTurbulence(TurbulenceParticleSystem& p, FlagGrid& flags, WaveletNoiseField& noise, Grid<Real>& kGrid, 
                            Real alpha, Real dt, int octaves, Real scale, Real invL0) {
    const Real PERSISTENCE = 0.56123f;
    
    const Vec3 pos(p[i].pos);
    if (flags.isInBounds(pos)) { // && !flags.isObstacle(pos)) {
        Real k2 = kGrid.getInterpolated(pos);
        Real ks = k2<0 ? 0.0 : sqrt(k2);
        
        // Wavelet noise lookup
        Real amplitude = scale * ks;
        Real multiplier = invL0;
        Vec3 vel(0.);        
        for (int o=0; o<octaves; o++) {
            //Vec3 ns = noise.evaluateCurl(p[i].pos * multiplier) * amplitude;
            Vec3 n0 = noise.evaluateCurl(p[i].tex0 * multiplier) * amplitude;
            Vec3 n1 = noise.evaluateCurl(p[i].tex1 * multiplier) * amplitude;
            vel += alpha * n0 + (1.0f-alpha) * n1;
            
            // next scale 
            amplitude *= PERSISTENCE;
            multiplier *= 2.0f;
        }
        
        // advection
        Vec3 dx = vel*dt;
        p[i].pos += dx;
        p[i].tex0 += dx;
        p[i].tex1 += dx;
    }
}
    
void TurbulenceParticleSystem::synthesize(FlagGrid& flags, Grid<Real>& k, int octaves, Real switchLength, Real L0, Real scale, Vec3 inflowBias) {
    static Real ctime = 0;
    static Vec3 inflow(0.);
    Real dt = getParent()->getDt();
    
    // collect inflow bias
    inflow += inflowBias * dt;
    
    // alpha: hat function over time
    Real oldAlpha = 2.0f*nmod(ctime/switchLength, 1.0f);        
    ctime += dt;
    Real alpha = 2.0f*nmod(ctime/switchLength, 1.0f);
    
    if (oldAlpha < 1.0f && alpha >= 1.0f) resetTexCoords(0, inflow);
    if (oldAlpha > alpha) resetTexCoords(1, inflow);
    if (alpha>1.0f) alpha=2.0f-alpha;
    alpha = 1.0;
        
    KnSynthesizeTurbulence(*this, flags, noise, k, alpha, dt, octaves, scale, 1.0f/L0);
}
   

} // namespace
