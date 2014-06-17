/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * CUDA functions for meshes
 *
 ******************************************************************************/
 
#include <iostream>
#include <map>
#include "cudatools.h"
#include "vortexsheet.h"

using namespace std;

namespace Manta {


const int VortKernelBlockSize = 512;

//! for each Node, integrate all triangles
//! [for simplicity, expects array with multiples of VortKernelBlockSize]
__global__ void VorticityKernel(CVec3Ptr nodes, const CVec3Ptr vortexPos, const CVec3Ptr strength, float reg2, float cutoff2, 
								const int len, const int triBlocks, const int offset) {
	const int nodeIdx = threadIdx.x + VortKernelBlockSize * blockIdx.x + offset;
	
	__shared__ float3 vpos[VortKernelBlockSize];
	__shared__ float3 vstr[VortKernelBlockSize];

	// load current position
	const int triLocal = threadIdx.x;
	float3 pos = nodes.get(nodeIdx);
	float3 u = pos;
		
	// divide triangles into blocks for shared memeory usage
	for (int i=0; i<triBlocks; i++) {
		const int triGlobal = threadIdx.x + VortKernelBlockSize * i;
		
		// load shared data
		vpos[triLocal] = vortexPos.get(triGlobal);
		vstr[triLocal] = strength.get(triGlobal);
		__syncthreads();
				
		//if (nodeIdx < len) {
			for (int j=0; j<VortKernelBlockSize; j++) {
				// actual vorticity kernel
				float3 r = pos - vpos[j];
				float r2 = normSqr(r);
				if (r2 < cutoff2) {
					float l = r2+reg2, l2=l*l;
					float div = rsqrtf(l2*l);
					u += cross (vstr[j], r) * div;
				}
			}
		//}
		__syncthreads();
	}
	
	// writeback
	nodes.set(nodeIdx, u);    
}

PYTHON void meshApplyBuoyancyTotalCuda(VortexSheetMesh& mesh, 
								  Real scale=1e-3, Real regularization=1, Real cutoffCells = 1e10, bool useDiff = false)
{
	Real dt = parent->getDt();
	// align len to blocksize
	int nodelen = ((mesh.numNodes() - 1) / VortKernelBlockSize + 1) * VortKernelBlockSize;
	int trilen = ((mesh.numTris() - 1) / VortKernelBlockSize + 1) * VortKernelBlockSize;
		
	// construct device arrays
	CVec3Array nodes(nodelen);
	CVec3Array tripos(trilen);
	CVec3Array strength(trilen);
	
	// copy values from mesh
	for(size_t i=0; i<mesh.numTris(); i++) {
		Vec3 center = mesh.getFaceCenter(i);        
		Vec3 vort = mesh.sheet(i).vorticity;
		if (useDiff)
			vort -= mesh.sheet(i).vorticitySmoothed;
		Vec3 str = vort * (mesh.getFaceArea(i) * scale * dt);
		
		tripos.set(i, center);
		strength.set(i, str);
	}
	for(size_t i=0; i<mesh.numNodes(); i++)
		nodes.set(i, mesh.nodes(i).pos);
	
	// fill aligned parts
	for(int i=mesh.numTris(); i<trilen; i++) { tripos.set(i, Vec3::Zero); strength.set(i, Vec3::Zero); }
	for(int i=mesh.numNodes(); i<nodelen; i++) { nodes.set(i, Vec3::Zero); }

	nodes.upload();
	tripos.upload();
	strength.upload();
	
	// construct parameters
	float reg2=regularization*regularization;
	float cutoff2 = cutoffCells*cutoffCells;
	
	// to avoid timeout, divide launches
	const int diff = 15000000 / tripos.size() * VortKernelBlockSize;
	cout << "Mesh buoyancy start" << endl;
	for (int offset=0; offset<nodelen; offset+= diff) {
		int blocks = (min(nodelen-offset, diff) - 1) / VortKernelBlockSize + 1;
		
		// invoke kernel
		VorticityKernel<<<blocks, VortKernelBlockSize>>> (nodes.data(), tripos.data(), strength.data(), reg2, cutoff2, mesh.numNodes(), trilen/VortKernelBlockSize, offset);
	}    
	cout << "Mesh buoyancy: " << (nodelen/diff)+1 << " calls" << endl;
	nodes.download();
	
	// readback new node positions
	for(size_t i=0; i<mesh.numNodes(); i++) {
		if (!mesh.isNodeFixed(i))
			mesh.nodes(i).pos = nodes[i];
	}
}


const int FvnBlockSize = 512;
const int FastEvalGridDivider = 2;
	
//! for each Node, integrate all triangles
//! use block cell structure for fast pruning
__global__ void FastVorticityKernelNoshare(CVec3Ptr nodes, const CVec3Ptr srcPos, const CVec3Ptr strength, const int* srcStart, const int* srcLen,
										   float reg2, float cutoff2, float safeRadius2, int intRadius, int3 gridRes, float maxres, const int len, const int divider, const int offset) 
{
	const int nodeIdx = threadIdx.x + FvnBlockSize * blockIdx.x + offset;
	if (nodeIdx >= len) return;
	
	// load current node position
	float3 pos = nodes.get(nodeIdx);
	int3 cell = make_int3((int)pos.x / divider, (int)pos.y / divider, (int)pos.z / divider);
	float3 u = pos;
	
	// query cells within block radius
	int3 minBlock = make_int3(max(cell.x - intRadius,0), max(cell.y - intRadius,0), max(cell.z - intRadius,0)); 
	int3 maxBlock = make_int3(min(cell.x + intRadius, gridRes.x - 1), min(cell.y + intRadius, gridRes.y - 1), min(cell.z + intRadius, gridRes.z - 1)); 
	for (int i=minBlock.x; i<=maxBlock.x; i++)
		for (int j=minBlock.y; j<=maxBlock.y; j++)
			for (int k=minBlock.z; k<=maxBlock.z; k++) {
				
				// test if block is within radius
				float3 d = make_float3(cell.x-i, cell.y-j, cell.z-k);
				if (normSqr(d) > safeRadius2) continue;
				
				// find source cell, and divide it into thread blocks
				int block = i + gridRes.x * (j + gridRes.y * k);
				int slen = srcLen[block]; 
				if (slen == 0) continue;
				int start = srcStart[block];
				
				// process sources
				for(int s=0; s<slen; s++) {                    
					// actual vorticity kernel
					float3 r = pos - srcPos.get(start+s);
					float r2 = normSqr(r);
					if (r2 < cutoff2) {
						float l = r2+reg2, l2=l*l;
						float div = rsqrtf(l2*l);
						u += cross (strength.get(start+s), r) * div;
					}                                
				}
			}
			
			
	// writeback
	nodes.set(nodeIdx, u);
}

inline int cIndex(const Vec3& pos, const Vec3i& s) {
	Vec3i p = toVec3i(pos) / FastEvalGridDivider;
	if (p.x < 0 || p.y < 0 || p.z < 0 || p.x >= s.x || p.y >= s.y || p.z >= s.z) return -1;
	return p.x + s.x * (p.y + s.y * p.z);
}

// TODO: don't reorder nodes -- performance ?
PYTHON void meshApplyBuoyancyLocalCuda(VortexSheetMesh& mesh, 
							 Real scale=1e-3, int cutoffCells=5, Real regularization=1)
{        
	Real dt = parent->getDt();
	
	// prepare center values and strength
	vector<Vec3> center(mesh.numTris());
	vector<Vec3> strength(mesh.numTris());
	for(size_t i=0; i<mesh.numTris(); i++) {
		Vec3 vort = mesh.sheet(i).vorticity - mesh.sheet(i).vorticitySmoothed;
		strength[i] = vort * (mesh.getFaceArea(i) * scale * dt);
		center[i] = mesh.getFaceCenter(i);        
	}
		
	// create triangles(sources) lookup grid    
	Vec3i gridRes = parent->getGridSize() / FastEvalGridDivider;
	const int numCells = gridRes.x * gridRes.y * gridRes.z;
	
	// 1. count sources per cell
	CArray<int> srcPerCell(numCells);
	Real maxres = gridRes.max();
	for (size_t i=0; i<center.size(); i++) {
		int cell = cIndex(center[i], gridRes);
		if (cell >= 0)
			srcPerCell[cell]++;
	}
	srcPerCell.upload();
	
	// 2. create start index lookup
	CArray<int> srcCellStart(numCells);
	int cnt=0;
	for (int i=0; i<numCells; i++) {
		srcCellStart[i] = cnt;
		cnt += srcPerCell[i];
	}
	srcCellStart.upload();
	
	// 3. reorder sources
	CVec3Array reorderStrength(center.size());
	CVec3Array reorderSourcePos(center.size());
	{
		vector<int> curSrcCell(numCells);
		for (int i=0; i<(int)center.size(); i++) {
			int cell = cIndex(center[i], gridRes);
			if (cell < 0) continue;
			int idx = srcCellStart[cell] + curSrcCell[cell];
			reorderStrength.set(idx, strength[i]);
			reorderSourcePos.set(idx, center[i]);
			curSrcCell[cell]++;
		}
	}
	reorderStrength.upload();
	reorderSourcePos.upload();
   
	// group nodes into blocks
	// 1. count nodes per cell
	vector<int> nodesPerCell(numCells);
	for (int i=0; i<mesh.numNodes(); i++) {
		int cell = cIndex(mesh.nodes(i).pos, gridRes);
		if (cell >= 0)
			nodesPerCell[cell]++;
	}    
	// 2. cluster blocks into execution plan
	vector<int> nodesCellStart(numCells);
	int offset = 0;
	for (int i=0; i<numCells; i++) {
		nodesCellStart[i] = offset;
		offset += nodesPerCell[i];
	}
	// 3. reorder nodes
	CVec3Array reorderNodes(mesh.numNodes());
	{
		vector<int> curNodeCell(numCells);
		for (int i=0; i<mesh.numNodes(); i++) {
			int cell = cIndex(mesh.nodes(i).pos, gridRes);
			if (cell < 0) continue;
			int idx = nodesCellStart[cell] + curNodeCell[cell];
			reorderNodes.set(idx, mesh.nodes(i).pos);
			curNodeCell[cell]++;
		}        
	}
	reorderNodes.upload();
	
	// construct parameters
	int cutoffInCells = cutoffCells / FastEvalGridDivider; // translate to lookup grid
	float safeRadius = (float)cutoffCells / (float)FastEvalGridDivider + sqrt(3.0);    
	float safeRadius2 = safeRadius*safeRadius;    
	float reg2 = regularization*regularization;
	float cutoff2 = cutoffCells*cutoffCells;
	cout << "cutoff int " <<cutoffInCells << " safe " << safeRadius << " cutoff " << cutoffCells << endl;
	
	// call in chunks for prevent timeout
	const int MAXBLOCK = 100;
	const int diff = MAXBLOCK*FvnBlockSize;
	int rcnt = 0;
	for (int offset=0; offset<(int)reorderNodes.size(); offset+= diff, rcnt++) {
		int blocks = (std::min(reorderNodes.size() - offset,diff)-1)/FvnBlockSize + 1;
	
		FastVorticityKernelNoshare<<<blocks,FvnBlockSize>>>
			(reorderNodes.data(), reorderSourcePos.data(), reorderStrength.data(), srcCellStart.data(), srcPerCell.data(),
			reg2, cutoff2, safeRadius2, cutoffInCells, make_int3(gridRes.x, gridRes.y, gridRes.z), maxres, reorderNodes.size(), FastEvalGridDivider, offset);
	}
	cout << "maxblocks " << MAXBLOCK << " total calls " << rcnt << endl;
	
	// download and reorder
	reorderNodes.download();    
	{
		vector<int> curNodeCell(numCells);        // redo ordering
		for (int i=0; i<mesh.numNodes(); i++) {
			int cell = cIndex(mesh.nodes(i).pos, gridRes);
			if (cell < 0) continue;
			int idx = nodesCellStart[cell] + curNodeCell[cell];
			if (!mesh.isNodeFixed(i))
				mesh.nodes(i).pos = reorderNodes[idx];
			curNodeCell[cell]++;
		}        
	}
}

#define MAXFRONT 500

__global__ void GaussKernel(CVec3Ptr distance, int* con, CVec3Ptr vort, CVec3Ptr vortSmooth, int blocksize, int stride, float cutoff, float mult, int offset) {
	const int Cidx = threadIdx.x + blocksize * blockIdx.x + offset;
	if (Cidx >= stride) return;
	
	int kernelIdx[MAXFRONT];
	float kernelDist[MAXFRONT];
	for(int i=0; i<MAXFRONT; i++) {
		kernelIdx[i] = -1;
		kernelDist[i] = 0;        
	}
	int kernelElements = 1, newKernelElements = 1, lastPos = 0;
	kernelIdx[0] = Cidx;
	
	//float a =0;
	int iter;
	for(iter=0; iter < 40; iter++) {
		for (int i=lastPos; i<kernelElements; i++) {
			
			int curidx = kernelIdx[i];
			float curdist = kernelDist[i];
			
			// add adjacent triangles to front list
			float nextdist[3];
			nextdist[0] = distance.x[curidx];
			nextdist[1] = distance.y[curidx];
			nextdist[2] = distance.z[curidx];
			for (int c=0; c<3; c++) {
				int nextidx = con[c+3*curidx];
				if (nextidx < 0) continue;
				float dist = nextdist[c] + curdist; 
				if (dist > cutoff) continue;
				
				// check if already in list
				bool found = false;
				for (int j=0; j<newKernelElements; j++) {
					if (kernelIdx[j] == nextidx) {
						found = true;
						if (kernelDist[j] > dist) kernelDist[j] = dist;
						break;                        
					}                    
				}
				if (!found) {
					if (newKernelElements >= MAXFRONT-1) goto finish;
					kernelIdx[newKernelElements] = nextidx;
					kernelDist[newKernelElements] = dist;
					newKernelElements++;
				}                
			}              
		}
		if (kernelElements == newKernelElements) break;
		lastPos = kernelElements;
		kernelElements = newKernelElements;
	}
	finish:
	
	// run gauss kernel over all nodes
	float sum = 0;
	float3 smooth = make_float3(0,0,0);
	for (int j=0; j<newKernelElements; j++) {
		const int idx = kernelIdx[j];
		float dist = kernelDist[j];
		float coeff = exp(dist*dist*mult);
		sum += coeff;
		smooth += make_float3(vort.x[idx] * coeff, vort.y[idx] * coeff, vort.z[idx] * coeff);        
	}        
	vortSmooth.set(Cidx, make_float3(smooth.x/sum, smooth.y/sum, smooth.z/sum));
}

PYTHON void filterVorticityCuda(VortexSheetMesh& mesh, Real sigma) {
	const int len = mesh.numTris();
	
	// upload mesh properties
	CArray<int> connectivity(len*3);
	CVec3Array faceCenter(len), distance(len), vorticity(len), vortSmooth(len);
	for (int i=0; i<len; i++) {
		faceCenter.set(i, mesh.getFaceCenter(i));
		vorticity.set(i, mesh.sheet(i).vorticity);
	}
	faceCenter.upload();
	vorticity.upload();
	vortSmooth.upload();
	for (int i=0; i<len; i++) {
		Vec3 dist;
		for (int c=0; c<3; c++) {
			int ot = mesh.corners(mesh.corners(i,c).opposite).tri;
			connectivity[i*3+c] = ot;
			dist[c] = norm(faceCenter[i] - faceCenter[ot]);
		}
		distance.set(i, dist);
	}
	distance.upload();
	connectivity.upload();
	
	const float cutoff = 2.0*sigma;
	const float mult = -0.5/sigma/sigma;
	
	// invoke kernel in blocks (to avoid GUI stall)
	const int blockSize = 128;
	const int numblocks = 100;
	const int maxProcess = blockSize * numblocks;
	int cnt=0;
	for (int offset=0; offset<len; offset+= maxProcess, cnt++) {
		int blocks = (min(len - offset, maxProcess)-1)/blockSize + 1;
		
		GaussKernel<<<blocks,blockSize>>>(distance.data(), connectivity.data(), vorticity.data(), vortSmooth.data(), blockSize, len, cutoff, mult, offset);
	}   
	cout << "Maxblocks: " << numblocks << ", Totalblocks: " << faceCenter.size()/blockSize << " -> " << cnt << " calls in total" << endl;
	
	// download and set
	vortSmooth.download();
	for (int i=0; i<len; i++) {
		mesh.sheet(i).vorticitySmoothed = vortSmooth[i];
	}    
}

} // namespace