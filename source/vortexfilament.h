/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex filament
 *
 ******************************************************************************/

#ifndef _VORTEXFIL_H
#define _VORTEXFIL_H

#include "particle.h"

namespace Manta {
class Mesh;
    
struct VortexFilamentData {
    VortexFilamentData() : pos(0.0f),vorticity(0.0f),sigma(0),flag(0) {}
    VortexFilamentData(const Vec3& p, const Vec3& v, Real sig) : pos(p),vorticity(v),sigma(sig),flag(0) {}
    Vec3 pos, vorticity;
    Real sigma;
    int flag;    
    static ParticleBase::SystemType getType() { return ParticleBase::VORTEX; }
};

//! Vortex particles
PYTHON class VortexFilamentSystem : public ParticleSystem<VortexFilamentData> {
public:
    PYTHON VortexFilamentSystem(FluidSolver* parent);
  
    PYTHON void advectSelf(Real scale=1.0, int integrationMode=RK4);
    PYTHON void applyToMesh(Mesh& mesh, Real scale=1.0, int integrationMode=RK4);
};

} // namespace


#endif