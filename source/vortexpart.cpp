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
    VortexKernel() {}
    VortexKernel(VortexParticleData& d, Real scale) : pos(d.pos) {
        isigma = -0.5/square(d.sigma);
        vortNorm = d.vorticity;
        strength = normalize(vortNorm) * scale;
        cutoff2 = 6.0 * square(d.sigma);
    }
    Vec3 pos;
    Vec3 vortNorm;
    Real strength;
    Real isigma;
    Real cutoff2;
    
    inline bool isValid(const Vec3& p) const {
        const Real rlen2 = normSquare(p-pos);        
        return rlen2 < cutoff2 && rlen2 > 1e-8;
    }
    
    inline Vec3 evaluate(const Vec3& p) const {
        // transform in cylinder coordinate system
        const Vec3 r = p-pos;
        const Real rlen2 = normSquare(r);        
        const Real rlen = sqrt(rlen2);
        const Real z = dot(r, vortNorm);
        const Vec3 ePhi = cross(r, vortNorm) / rlen;
        const Real rho2 = rlen2 - z*z;
        
        Real vortex;
        if (rho2 > 1e-10) {
            // evaluate Kernel      
            vortex = strength * sqrt(rho2) * exp (rlen2 * isigma);  
        } else {
            vortex = 0;
        }
        return vortex * ePhi;
    }
};
    
VortexParticleSystem::VortexParticleSystem(FluidSolver* parent) :
    ParticleSystem<VortexParticleData>(parent)
{ 
}

DefineIntegrator(integrateVortexKernel, VortexKernel, evaluate);

KERNEL(particle) template<IntegrationMode mode> 
advectNodes(vector<Vec3>& nodesNew, const vector<Vec3>& nodesOld, const VortexKernel& kernel, Real dt) {
    const Vec3 p = nodesOld[i];
    if (kernel.isValid(p))
        nodesNew[i] += integrateVortexKernel<mode>(p, kernel, dt);    
}

KERNEL(particle) template<IntegrationMode mode>
advectVortices(vector<VortexParticleData>& nodesNew, const vector<VortexKernel>& nodesOld, const VortexKernel& kernel, Real dt) {
    const Vec3 p = nodesOld[i].pos;
    if (kernel.isValid(p))
        nodesNew[i].pos += integrateVortexKernel<mode>(p, kernel, dt);    
}

void VortexParticleSystem::applyToMesh(Mesh& mesh, Real scale, int integrationMode) {
    const Real dt = getParent()->getDt();
    
    // copy node array
    const int nodes = mesh.numNodes();
    vector<Vec3> nodesOld(nodes), nodesNew(nodes);
    vector<VortexKernel> kernels(nodes);
    for (int i=0; i<nodes; i++) {
        nodesOld[i] = nodesNew[i] = mesh.nodes(i).pos;
    }
    for (size_t i=0; i<mData.size(); i++) {
        kernels[i] = VortexKernel(mData[i], scale);
    }
    
    // loop over the vortices
    for (size_t i=0; i<mData.size(); i++) {
        if (!isActive(i)) continue;
        if (integrationMode==EULER) {
            advectNodes<EULER>(nodesNew, nodesOld, kernels[i], dt);
            advectVortices<EULER>(mData, kernels, kernels[i], dt);
        } else if (integrationMode==RK2) {
            advectNodes<RK2>(nodesNew, nodesOld, kernels[i], dt);
            advectVortices<RK2>(mData, kernels, kernels[i], dt);
        } else if (integrationMode==RK4) {
            advectNodes<RK4>(nodesNew, nodesOld, kernels[i], dt);
            advectVortices<RK4>(mData, kernels, kernels[i], dt);
        } else throw Error("unknown integration type");
    }
    
    // copy back
    for (int i=0; i<nodes; i++) {
        if (!mesh.isNodeFixed(i))
            mesh.nodes(i).pos = nodesNew[i];
    }    
}

    

} // namespace
