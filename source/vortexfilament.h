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
    
struct VortexRing {
    VortexRing() : circulation(0.),flag(0),isClosed(false) {}
    VortexRing(Real c, bool closed=false) : circulation(c),flag(0),isClosed(closed) {}
    void renumber(int* _renumber);
    inline int size() const { return indices.size(); }
    inline int idx(int i) const { return indices[(i+indices.size()) % indices.size()]; }
    inline int idx0(int i) const { return indices[i]; }
    inline int idx1(int i) const { return indices[ (i+1) % indices.size() ]; }
    
    bool isClosed;
    int flag;
    Real circulation;
    std::vector<int> indices;
};

//! Vortex filaments
PYTHON class VortexFilamentSystem : public ConnectedParticleSystem<BasicParticleData, VortexRing> {
public:
    virtual SystemType getType() const { return ParticleBase::FILAMENT; };
        
    PYTHON VortexFilamentSystem(FluidSolver* parent);
  
    //! self-advect the filament system
    PYTHON void advectSelf(Real scale=1.0, Real regularization=0.1, int integrationMode=IntRK4);
    //! advect a particle system 
    PYTHON void advectParticles(TracerParticleSystem& sys, Real scale=1.0, Real regularization=0.1, int integrationMode=IntRK2);
    //! advect triangle mesh using filaments
    PYTHON void advectMesh(Mesh& mesh, Real scale=1.0, Real regularization=0.1, int integrationMode=IntRK4);
    //! perform doubly-discrete smoke ring flow update
    //! as in [Weissmann,Pinkall 2009]
    PYTHON void doublyDiscreteUpdate(Real regularization=0.1);
    //! remesh long or strongly-curved segments
    PYTHON void remesh(Real maxLen=3.0, Real minLen=1.0);
    
    //! add a filament ring to the system
    PYTHON void addRing(const Vec3& position, Real circulation, Real radius, Vec3 normal, int number);
    //! add a line filament to the system
    PYTHON void addLine(const Vec3& p0, const Vec3& p1, Real circulation);
    
        
    virtual ParticleBase* clone();
protected:
    
    //! Biot-Savart line integration
    void integrate(const std::vector<Vec3>& nodesOld, std::vector<Vec3>& nodesNew, Real scale, Real reg, int integrationMode);
};

} // namespace


#endif