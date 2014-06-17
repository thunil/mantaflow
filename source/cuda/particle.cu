/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * GPU SDF creation from triangle mesh
 * 
 ******************************************************************************/

#include "cudatools.h"
#include "particle.h"
#include "grid.h"

/*using namespace std;
namespace Manta {
using namespace std;
using DDF::nVec3i;
using DDF::Vec3;

// vortex particle effect: (cyl coord around wp)
// u = -|wp|*rho*exp( (-rho^2-z^2)/(2sigma^2) ) e_phi
__device__ inline float3 velocityKernel(const float3& r, const float3& vortNorm, const float strength, const float isigma) {
	// transform in cylinder coordinate system
	const float rlen2 = normSqr(r);
	const float rlendiv = rsqrtf(rlen2);
	const float z = dot(r, vortNorm);
	const float3 ePhi = cross(r, vortNorm) * rlendiv;
	const float rho2 = rlen2 - z*z;
	
	float vortex;
	if (rho2 > 1e-10) {
		// evaluate Kernel      
		vortex = strength * sqrtf(rho2) * expf (rlen2 * isigma);  
	} else {
		vortex = 0;
	}
	return vortex * ePhi;
}

__device__ inline float3 rk4(const float3& p, const float3& vortNorm, const float strength, const float isigma, const float dt) {
	float3 k1 = velocityKernel(p, vortNorm, strength, isigma);
	float3 p2 = p + k1 * (0.5*dt);
	float3 k2 = velocityKernel(p2, vortNorm, strength, isigma);
	float3 p3 = p + k2 * (0.5*dt);
	float3 k3 = velocityKernel(p3, vortNorm, strength, isigma);
	float3 p4 = p + dt*k3;
	float3 k4 = velocityKernel(p4, vortNorm, strength, isigma);   
	return (k1 + (k2+k3)*2.0 + k4) * (dt/6.0);
}

const int VortKernelBlockSize = 256;

// apply vortex kernel to nodes
__global__ void integrateNodesRK4 (AlignedVec3Array nodes, const AlignedVec3Array sources, const AlignedVec3Array vortex, const AlignedFloatArray sigma2,
								   const float dt)
{
	const int globalIdx = threadIdx.x + VortKernelBlockSize * blockIdx.x;
	const int localIdx = threadIdx.x;
	__shared__ float3 vpos[VortKernelBlockSize];
	__shared__ float3 vort[VortKernelBlockSize];
	__shared__ float vstr[VortKernelBlockSize];
	__shared__ float vsig[VortKernelBlockSize];

	// load current position
	float3 pos = nodes.get(globalIdx);
	float3 u = pos;
		
	// divide sources into blocks for shared memeory usage
	for (int i=0; i<sources.blocks; i++) {
		const int blockIdx = VortKernelBlockSize * i;
		const int srcGlobal = localIdx + blockIdx;
		
		// load shared data
		vpos[localIdx] = sources.get(srcGlobal);
		vsig[localIdx] = sigma2.get(srcGlobal);
		float3 v = vortex.get(srcGlobal);
		float vnorm = normalize(v);
		vstr[localIdx] = vnorm;
		vort[localIdx] = v;
		__syncthreads();

		if (globalIdx < nodes.len) {            
			for (int j=0; j<VortKernelBlockSize; j++) {
				if (j+blockIdx >= sources.len) continue;
				
				// apply actual vorticity kernel
				float3 r = pos - vpos[j];
				float sig2 = vsig[j];
				float r2 = normSqr(r);
				float cutoff2 = 6.0f*sig2;
				if (r2 > cutoff2 || r2 < 1e-8) continue;
				// RK4 integration
				u += rk4(r, vort[j], vstr[j], -0.5/sig2, dt);                
			}
		}
		__syncthreads();
	}
	
	// writeback
	if (globalIdx < nodes.len)
		nodes.set(globalIdx, u);    
}

void cudaIntegrateVortexParticles(vector<SurfaceNode>& nodes, DDF::ParticleSystemVortex& sys, const float scale, const float dt, const float dx) 
{    
	// count valid source, create host array
	int numSources = 0;
	for(int i=0; i<sys.size(); i++)
		if (sys.isActive(i)) numSources++;
	
	vector<Vec3> srcPos(numSources), srcStr(numSources);
	vector<float> srcSig(numSources); 
	for (int i=0,idx=0; i<sys.size(); i++) {
		if (sys.isActive(i)) {
			srcPos[idx] = sys.pos(i);
			srcStr[idx] = sys.vortex(i).vorticity * scale;
			float sig = sys.vortex(i).sigma;
			srcSig[idx] = sig*sig;
			idx++;
		}
	}
	if (numSources == 0) return;
	
	// count valid dest, create host arrays
	int numDest = 0;
	for(size_t i=0; i<nodes.size(); i++)
		if ((nodes[i].flags & SurfaceNode::FIXED) == 0) numDest++;            
	
	vector<Vec3> dstPos(numDest);
	for(int i=0,idx=0; i<(int)nodes.size(); i++) {
		if ((nodes[i].flags & SurfaceNode::FIXED) == 0) {
			SmVector3 p = nodes[i].pos;
			dstPos[idx] = Vec3(p[0] / dx, p[1] / dx, p[2] / dx);
			idx++;
		}
	}
	
	// upload sources to GPU
	DeviceVec3 dSrcPos(VortKernelBlockSize), dSrcStr(VortKernelBlockSize), dDstPos(VortKernelBlockSize), dVDstPos(VortKernelBlockSize);
	DeviceFloat dSrcSig(VortKernelBlockSize);
	dSrcPos.upload(srcPos);
	dSrcSig.upload(srcSig);
	dSrcStr.upload(srcStr);
		
	// apply to mesh
	dDstPos.upload(dstPos);
	integrateNodesRK4<<<dDstPos.blocks(), VortKernelBlockSize>>>(dDstPos.device, dSrcPos.device, dSrcStr.device, dSrcSig.device, dt);
	dDstPos.errorTest();
	dDstPos.download(dstPos);
	
	// apply to vortex particles
	dVDstPos.upload(srcPos);
	integrateNodesRK4<<<dVDstPos.blocks(), VortKernelBlockSize>>>(dVDstPos.device, dSrcPos.device, dSrcStr.device, dSrcSig.device, dt);
	dVDstPos.errorTest();
	dVDstPos.download(srcPos);
	
	// back to arrays...
	for (int i=0,idx=0; i<sys.size(); i++)
		if (sys.isActive(i)) sys.setPos(i, srcPos[idx++]);

	for (int i=0,idx=0; i<(int)nodes.size(); i++) {
		if ((nodes[i].flags & SurfaceNode::FIXED) == 0) {
			nodes[i].pos = SmVector3(dstPos[idx].x * dx, dstPos[idx].y * dx, dstPos[idx].z * dx);
			idx++;            
		}
	}    
}
*/