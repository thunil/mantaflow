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
    VortexFilamentData() : idx0(-1),idx1(-1),circulation(0),flag(0) {}
    VortexFilamentData(int i0, int i1, Real c) : idx0(i0),idx1(i1),circulation(c),flag(0) {}
    void renumber(int* _renumber) { idx0 = _renumber[idx0]; idx1 = _renumber[idx1]; }
    
    int idx0, idx1;
    Real circulation;
    int flag;
};


//! Vortex filaments
PYTHON class VortexFilamentSystem : public ConnectedParticleSystem<BasicParticleData, VortexFilamentData> {
public:
    virtual SystemType getType() const { return ParticleBase::FILAMENT; };
        
    PYTHON VortexFilamentSystem(FluidSolver* parent);
  
    PYTHON void advectSelf(Real scale=1.0, Real regularization=0.1, Real cutoff=1e-6, int integrationMode=RK4);
    PYTHON void applyToMesh(Mesh& mesh, Real scale=1.0, Real regularization=0.1, Real cutoff=1e-6, int integrationMode=RK4);
    
    PYTHON void addRing(const Vec3& position, Real circulation, Real radius, const Vec3& normal, int number);
    PYTHON void addLine(const Vec3& p0, const Vec3& p1, Real circulation);
    
protected:
    std::vector<VortexFilamentData> mSegments;
};

} // namespace


#endif