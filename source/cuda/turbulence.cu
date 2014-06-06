/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * GPU turbulence synthesis
 * 
 ******************************************************************************/

#include "curlnoise.h"
#include "vortexsheet.h"
#include "grid.h"
#include "commonkernels.h"

using namespace std;
namespace Manta {

CompressedTex gNoiseTex;
	
//! synthesize uniform K41 curl noise onto node list
__global__ void KnSynthesizeK41(CVec3Ptr nodes, CVec3Ptr tex1, CVec3Ptr tex2, Real alpha, float* ke,
								CudaNoiseDev noise, const int octaves, const float iL0, const float str, const int len, const float kmin) 
{
	const int nodeIdx = threadIdx.x + blockDim.x * blockIdx.x;
	if (nodeIdx >= len) return;
	
	// load from array
	float3 p = nodes.get(nodeIdx);
	float3 tc1 = tex1.get(nodeIdx);
	float3 tc2 = tex2.get(nodeIdx);
	float k = ke ? ke[nodeIdx] : 1.0f;
	//tc1=p; // static
	k = k - kmin;
	k = (k<0.0f) ? 0.0f : sqrtf(k);
	
	// interpolate texcoords 
	/*float3 tc = alpha * tc1 + (1.0f-alpha) * tc2;
	float3 v = noise.synthesizeK41(gNoiseTex, tc, octaves, L0);*/
	
	// interpolate velocities
	float3 v1 = noise.synthesizeK41(gNoiseTex, tc1, octaves, iL0);
	float3 v2 = noise.synthesizeK41(gNoiseTex, tc2, octaves, iL0);
	float3 v = alpha * v1 + (1.0f-alpha) * v2;
	
	// apply
	float3 update = (str * k) * v;
	p += update;
	tc1 += update;
	tc2 += update;

	// writeback
	nodes.set(nodeIdx, p);
	tex1.set(nodeIdx, tc1);
	tex2.set(nodeIdx, tc2);
}

//! synthesize K41 curl noise onto mesh
PYTHON void synthesizeK41(VortexSheetMesh& mesh, Grid<Real>* k = NULL,
					 Real scale = 1.0, Real L0 = 0.1, int octaves=3, Real switchLength = 10.0, bool hardReset=false, Real minIntensity=0.1) 
{
	const int blockSize = 256;
	const int blocks = (mesh.numNodes()-1)/blockSize+1;
	const float dt = parent->getDt();
	const float str = dt * scale;
	const float kmin = 1.5 * square(minIntensity);
	
	// hat function over time
	static float ctime = 0;
	float oldAlpha = 2.0f*nmod(ctime/switchLength, 1.0f);        
	ctime += parent->getDt();
	float alpha = 2.0f*nmod(ctime/switchLength, 1.0f);
	if (hardReset) {
		if (oldAlpha > alpha) mesh.resetTex1();
		alpha = 1.0f;        
	} else {
		if (oldAlpha < 1.0f && alpha >= 1.0f) mesh.resetTex2();
		if (oldAlpha > alpha) mesh.resetTex1();
		if (alpha>1.0f) alpha=2.0f-alpha;
	}
	
	// create noise tex on first call
	static CudaNoiseTexture noise;
		
	// upload data
	CVec3Array nodes(mesh.numNodes());
	CVec3Array tc1(mesh.numNodes()), tc2(mesh.numNodes());
	CArray<float> ke(mesh.numNodes());
	for (int i=0; i<mesh.numNodes(); i++) {
		nodes.set(i, mesh.nodes(i).pos);
		tc1.set(i, mesh.tex1(i));
		tc2.set(i, mesh.tex2(i));
	}
	nodes.upload();
	tc1.upload();
	tc2.upload();
	if (k) {
		for (int i=0; i<mesh.numNodes(); i++)
			ke.set(i, k->getInterpolated(mesh.nodes(i).pos));
		ke.upload();
	}
	KnSynthesizeK41<<<blocks, blockSize>>> (nodes.data(), tc1.data(), tc2.data(), alpha, k ? (ke.data()) : 0, 
											noise.bind(gNoiseTex), octaves, 1.0f/L0, str, mesh.numNodes(), kmin);
	
	// download data
	nodes.download();
	tc1.download();
	tc2.download();
	for (int i=0; i<mesh.numNodes(); i++) {
		if (!mesh.isNodeFixed(i)) {
			mesh.nodes(i).pos = nodes[i];
			mesh.tex1(i) = tc1[i];
			mesh.tex2(i) = tc2[i];
		}
	}
}

//! Kernel: synthesize uniform K41 curl noise onto grid
__global__ void KnSynthesizeGridK41(CVec3Ptr nodes, CudaNoiseDev noise, const int octaves, const float iL0, const float str, int stridey, int stridez) 
{
	int3 cell = make_int3(threadIdx.x + blockDim.x*blockIdx.x, threadIdx.y + blockDim.y*blockIdx.y, threadIdx.z + blockDim.z*blockIdx.z);
	float3 p = make_float3(cell.x , cell.y , cell.z );
	
	// interpolate velocities
	float3 v = noise.synthesizeK41(gNoiseTex, p, octaves, iL0);
	
	// writeback
	nodes.set(cell.x + stridey* cell.y + stridez* cell.z, str*v);
}

//! synthesize K41 curl noise onto grid
PYTHON void synthesizeK41Grid(MACGrid& vel, MACGrid& dst, Real scale = 1.0, Real L0 = 0.1, int octaves=3) 
{    
	const float dt = parent->getDt();
	const float str = dt * scale;
	
	// create noise tex on first call
	static CudaNoiseTexture noise;
	
	// upload data
	const Vec3i gridRes = parent->getGridSize();
	const int num = vel.getSizeX()*vel.getSizeY()*vel.getSizeZ();
	CVec3Array nodes(num);
	nodes.upload();

	dim3 blocksize(8, 8, 8);
	dim3 blocks((gridRes.x-1)/8+1, (gridRes.y-1)/8+1, (gridRes.z-1)/8+1);    
	KnSynthesizeGridK41<<<blocks, blocksize>>> (nodes.data(), noise.bind(gNoiseTex), octaves, 1.0f/L0, str, vel.getStrideY(), vel.getStrideZ());
	
	nodes.download();    
	Grid<Vec3> center(parent);
	MACGrid mac(parent);
	for (int i=0; i<num; i++)
		center[i] = nodes[i];
	
	GetMAC(mac, center);
	dst = vel;
	dst += (Grid<Vec3>&)mac;
}
	
} // namespace