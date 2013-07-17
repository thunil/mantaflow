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

#ifndef _TURBULENCEPART_H_
#define _TURBULENCEPART_H_

#include "particle.h"
#include "noisefield.h"

namespace Manta {
class Shape;
    
struct TurbulenceParticleData {
    TurbulenceParticleData() : pos(_0),tex0(_0),tex1(_0),flag(0) {}
    TurbulenceParticleData(const Vec3& p) : pos(p),tex0(p),tex1(p),flag(0) {}
    Vec3 pos;
    Vec3 tex0, tex1;
    int flag;
    static ParticleBase::SystemType getType() { return ParticleBase::TURBULENCE; }
};

//! Turbulence particles
PYTHON class TurbulenceParticleSystem : public ParticleSystem<TurbulenceParticleData> {
public:
    PYTHON TurbulenceParticleSystem(FluidSolver* parent, WaveletNoiseField& noise);
  
    PYTHON void resetTexCoords(int num);    
    PYTHON void seed(Shape* source, int num);
    PYTHON void synthesize(FlagGrid& flags, Grid<Real>& k, int octaves=2, Real switchLength=10.0, Real L0=0.1, Real scale=1.0);
        
    virtual ParticleBase* clone();
    
private:
    WaveletNoiseField& noise;
};

} // namespace


#endif