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

#include "vortexpart.h"
#include "integrator.h"
#include "mesh.h"

using namespace std;
namespace Manta {




// ****************************************************************************
//  Vortex Filaments System
// ****************************************************************************

template<ParticleSystem::IntegrationMode mode>
void FilamentIntegration(vector<SurfaceNode>& nodes, vector<Vec3>& vortexNodes, const vector<Vec3>& vertexNormals, const Vec3& gamma0, const Vec3& gamma1, const Real a2, const Real cutoff2, const Real idx, const Real strength, const Real dt) {
	// integrate over mesh nodes
	const Vec3 gammaMid = (gamma1 + gamma0) * 0.5;
	Vec3 gamma = gamma1 - gamma0;
	normalize(gamma);
	
#if DDF_OPENMP==1
#pragma omp for schedule(static)
#endif
	for (size_t i=0; i<nodes.size(); i++) {
        if ((nodes[i].flags & SurfaceNode::FIXED) != 0) continue;
        const Vec3 p = Vec3(nodes[i].lastPos[0] * idx, nodes[i].lastPos[1] * idx, nodes[i].lastPos[2] * idx);
		const Vec3 r0 = p - gamma0, r1 = p - gamma1;
		if (normNoSqrt(p - gammaMid) >= cutoff2) continue;
		
		Vec3 du = Integrator<mode>::filament(r0, r1, vertexNormals[i], gamma, a2, strength, dt);
		nodes[i].pos += SmVector3(du.x, du.y, du.z) / idx;
	}
	
	// integrate over vortex nodes
#if DDF_OPENMP==1
#pragma omp for schedule(static)
#endif
	for (size_t i=0; i<vortexNodes.size(); i++) {
		const Vec3 r0 = vortexNodes[i] - gamma0, r1 = vortexNodes[i] - gamma1;
		if (normNoSqrt(vortexNodes[i] - gammaMid) >= cutoff2 || normNoSqrt(r0) < 1e-8 || normNoSqrt(r1) < 1e-8) continue;
		
		vortexNodes[i] += Integrator<mode>::filament(r0, r1, Vec3::ZERO, gamma, a2, strength, dt);
		/*if (vortices[i].x <= 0 || vortices[i].y <= 0 || vortices[i].z <= 0 || 
			vortices[i].x >= 60 || vortices[i].y >= 60 || vortices[i].z >= 60) {
			printf("p [%g,%g,%g] isigma %g strength %g\n",p.x,p.y,p.z, isigma, strength);
			Vec3 rdu = IntegrateEuler(p, isigma, vortNorm, strength, dt);exit(1);
		}*/
	}
}

void ParticleSystemFilaments::integrateNodes(vector<SurfaceNode>& nodes, const vector<Vec3>& normals,const Real scale, const double dt, const Real dx, IntegrationMode mode, FlagGrid* flaggrid) {
	// keep last pos	
	for(std::vector<SurfaceNode>::iterator sni=nodes.begin(); sni!=nodes.end(); sni++)
		sni->lastPos = sni->pos;
	vector<Vec3> lastPos (mPos);
	
	const Real a = 1, cutoff = 10;
	
	const Real idx = 1.0/dx, a2 = a*a, cutoff2 = cutoff*cutoff;
	
	// iterate filaments, mesh nodes
	for (size_t i=0; i< mFil.size(); i++) {
		if ((mFlags[i] & PDELETE) != 0) continue;
		
		// precalculate values for kernel evaluation
		
		const Vec3& gamma0 = lastPos[mFil[i].vertex1];
		const Vec3& gamma1 = lastPos[mFil[i].vertex2];
		const Real& strength = mFil[i].circulation * scale;
		
		// integrate mesh nodes and vortex filaments
		switch(mode) {
			case EULER: FilamentIntegration<EULER>(nodes, mPos, normals, gamma0, gamma1, a2, cutoff2, idx, strength, dt); break;
			case RK2: FilamentIntegration<RK2>(nodes, mPos, normals, gamma0, gamma1, a2, cutoff2, idx, strength, dt); break;
			case RK4: FilamentIntegration<RK4>(nodes, mPos, normals, gamma0, gamma1, a2, cutoff2, idx, strength, dt); break;
		}
	}
	
	// kill stray particles
	for (size_t i=0; i< mPos.size(); i++) {
		if (!flaggrid->checkIndexValidWithBounds(mPos[i].x,mPos[i].y,mPos[i].z,1))
		{}//mFlags[i] |= PDELETE;
	}
}

void ParticleSystemFilaments::addRing(const Vec3& position, Real circulation, const Vec3& radius, int number) {
	Vec3 normal(radius);
	Real rad = normalize(normal);
	
	Vec3 worldup (0,1,0);
	if (norm(normal - worldup) < 1e-5) worldup = Vec3(1,0,0);
	
	Vec3 u = cross(normal, worldup); normalize(u);
	Vec3 v = cross(normal, u); normalize(v);
	
	int firstNum = mPos.size();
	for (int i=0; i<number; i++) {
		Real phi = (Real)i/(Real)number * M_PI * 2.0;
		Vec3 p = position + rad * (u*cos(phi) + v*sin(phi));
		
		int num = ParticleSystem::add(p, PFILAMENT);
		mFil.push_back(FilamentType(num, num+1, circulation));
	}
	mFil.rbegin()->vertex2 = firstNum;
}

void ParticleSystemFilaments::addLine(const Vec3& p0, const Vec3& p1, Real circulation) {
	FilamentType fil(add(p0,PFILAMENT),add(p1,PFILAMENT),circulation);
	mFil.push_back(fil);
}


// vortex filament effect: 
inline Vec3 filamentKernel(const Vec3& r0, const Vec3& r1, const Vec3& n, const Vec3& gamma, const Real a2, const Real strength) {
	// vortex line integral
	const Vec3 ep = cross(gamma, r0);
	const Vec3 r0n = r0 / (a2 + normNoSqrt(r0));
	const Vec3 r1n = r1 / (a2 + normNoSqrt(r1));	
	const Real mag = dot(r1n - r0n, gamma) / (a2 + normNoSqrt(ep)) * strength;
	Vec3 fil = mag * ep;
	
	return fil;
}

// vortex filament effect
struct FilamentKernel {
    FilamentKernel() {}
    FilamentKernel(VortexFilamentData& d, Real scale) : pos(d.pos) {
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
void advectNodes(vector<Vec3>& nodesNew, const vector<Vec3>& nodesOld, const VortexKernel& kernel, Real dt) {
    const Vec3 p = nodesOld[i];
    if (kernel.isValid(p))
        nodesNew[i] += integrateVortexKernel<mode>(p, kernel, dt);    
}

KERNEL(particle) template<IntegrationMode mode>
void advectVortices(vector<VortexParticleData>& nodesNew, const vector<VortexKernel>& nodesOld, const VortexKernel& kernel, Real dt) {
    const Vec3 p = nodesOld[i].pos;
    if (kernel.isValid(p))
        nodesNew[i].pos += integrateVortexKernel<mode>(p, kernel, dt);    
}

void VortexParticleSystem::advectSelf(Real scale, int integrationMode) {
    const Real dt = getParent()->getDt();
    
    // copy kernel array
    vector<VortexKernel> kernels(size());
    for (size_t i=0; i<mData.size(); i++)
        kernels[i] = VortexKernel(mData[i], scale);
    
    // loop over the vortices
    for (size_t i=0; i<mData.size(); i++) {
        if (!isActive(i)) continue;
        if (integrationMode==EULER) advectVortices<EULER>(mData, kernels, kernels[i], dt);
        else if (integrationMode==RK2) advectVortices<RK2>(mData, kernels, kernels[i], dt);
        else if (integrationMode==RK4) advectVortices<RK4>(mData, kernels, kernels[i], dt);
        else errMsg("unknown integration type");
    }
}

void VortexParticleSystem::applyToMesh(Mesh& mesh, Real scale, int integrationMode) {
    const Real dt = getParent()->getDt();
    
    // copy node array
    const int nodes = mesh.numNodes();
    vector<Vec3> nodesOld(nodes), nodesNew(nodes);
    for (int i=0; i<nodes; i++)
        nodesOld[i] = nodesNew[i] = mesh.nodes(i).pos;
    
    // loop over the vortices
    for (size_t i=0; i<mData.size(); i++) {
        if (!isActive(i)) continue;
        VortexKernel kernel(mData[i], scale);
        if (integrationMode==EULER) advectNodes<EULER>(nodesNew, nodesOld, kernel, dt);
        else if (integrationMode==RK2) advectNodes<RK2>(nodesNew, nodesOld, kernel, dt);
        else if (integrationMode==RK4) advectNodes<RK4>(nodesNew, nodesOld, kernel, dt);
        else errMsg("unknown integration type");
    }
    
    // copy back
    for (int i=0; i<nodes; i++) {
        if (!mesh.isNodeFixed(i))
            mesh.nodes(i).pos = nodesNew[i];
    }    
}

    

} // namespace
