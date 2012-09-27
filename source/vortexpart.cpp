/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Vortex particles
 *
 ******************************************************************************/

#include "vortexpart.h"
#include "integrator.h"
#include "mesh.h"

using namespace std;
namespace Manta {

// vortex particle effect: (cyl coord around wp)
// u = -|wp|*rho*exp( (-rho^2-z^2)/(2sigma^2) ) e_phi
struct VortexKernel {
    VortexKernel(vector<VortexParticleData>& _vp, Real _scale) : vp(_vp), scale(_scale) {}
    
    inline Vec3 eval(const Vec3& p, VortexParticleData& orig, vector<Vec3>& pos) const {
        if (orig.flag & ParticleBase::PDELETE) return Vec3::Zero;
        return integrate(p,pos);
    }
    
    inline Vec3 eval(const Vec3& p, Node& orig, vector<Vec3>& pos) const {
        if (orig.flags & Mesh::NfFixed) return Vec3::Zero;
        return integrate(p,pos);
    }
    
    inline Vec3 integrate(const Vec3& p, vector<Vec3>& pos) const {
        Vec3 u(_0);
        for (size_t i=0; i<pos.size(); i++) {
            if (vp[i].flag & ParticleBase::PDELETE) continue;
            
            // cutoff radius
            const Vec3 r = p - pos[i];
            const Real rlen2 = normSquare(r);   
            const Real sigma2 = square(vp[i].sigma);
            if (rlen2 > 6.0 * sigma2 || rlen2 < 1e-8) continue;
            
            // split vortex strength
            Vec3 vortNorm = vp[i].vorticity;
            Real strength = normalize(vortNorm) * scale;
        
            // transform in cylinder coordinate system
            const Real rlen = sqrt(rlen2);
            const Real z = dot(r, vortNorm);
            const Vec3 ePhi = cross(r, vortNorm) / rlen;
            const Real rho2 = rlen2 - z*z;
        
            Real vortex = 0;
            if (rho2 > 1e-10) {
                // evaluate Kernel      
                vortex = strength * sqrt(rho2) * exp (rlen2 * -0.5/sigma2);  
            }
            u += vortex * ePhi;
        }
        return u;
    }
    
    const vector<VortexParticleData>& vp;
    const Real scale;
};
    
VortexParticleSystem::VortexParticleSystem(FluidSolver* parent) :
    ParticleSystem<VortexParticleData>(parent)
{ 
}

void VortexParticleSystem::advectSelf(Real scale, int integrationMode) {
    VortexKernel kernel(mData, scale);
    integratePointSet(mData, mData, kernel, integrationMode);    
}

void VortexParticleSystem::applyToMesh(Mesh& mesh, Real scale, int integrationMode) {
    VortexKernel kernel(mData, scale);
    integratePointSet(mData, mesh.getNodeData(), kernel, integrationMode);    
}

ParticleBase* VortexParticleSystem::clone() {
    VortexParticleSystem* nm = new VortexParticleSystem(getParent());
    compress();
    
    nm->mData = mData;
    nm->setName(getName());
    return nm;
}

    

} // namespace
