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
#include "mesh.h"
#include "grid.h"
#include <stack>

using namespace std;
namespace Manta {
	
const int SDFBlockSize = 8;

__global__ void SDFKernel(const int* partStart, const int* partLen, CVec3Ptr pos, CVec3Ptr normal, float* sdf, int3 gridRes, int intRadius, float safeRadius2, float cutoff2, float isigma2)
{
	// cell index, center
	int3 cell = make_int3(threadIdx.x + blockDim.x*blockIdx.x, threadIdx.y + blockDim.y*blockIdx.y, threadIdx.z + blockDim.z*blockIdx.z);
	if (cell.x >= gridRes.x || cell.y >= gridRes.y || cell.z >= gridRes.z) return;    
	float3 cpos = make_float3(cell.x + 0.5f, cell.y + 0.5f, cell.z + 0.5f);
	float sum = 0.0f;
	float dist = 0.0f;
	
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
				int slen = partLen[block];
				if (slen == 0) continue;
				int start = partStart[block];
								
				// process sources
				for(int s=0; s<slen; s++) {                     
					
					// actual sdf kernel
					float3 r = cpos - pos.get(start+s);
					float r2 = normSqr(r);
					if (r2 < cutoff2) {
						float w = expf(-r2*isigma2);
						sum += w;
						dist += dot(normal.get(start+s), r) * w;                        
					}
				}
			}
			
	// writeback
	if (sum > 0.0f)
		sdf[cell.x + gridRes.x * (cell.y + gridRes.y * cell.z)] = dist / sum;    
}

inline int _cIndex(const Vec3& pos, const Vec3i& s) {
	Vec3i p = toVec3i(pos);
	if (p.x < 0 || p.y < 0 || p.z < 0 || p.x >= s.x || p.y >= s.y || p.z >= s.z) return -1;
	return p.x + s.x * (p.y + s.y * p.z);
}

//! Obtain levelset from mesh.
//! This only works for dense meshes -- simply uses normals and triangle centers, no triangle integration
PYTHON void meshSDFCuda(Mesh& mesh, Grid<Real>& levelset, Real sigma, Real cutoff=-1)
{        
	if (cutoff<0) cutoff = 2*sigma;
	
	Vec3i gridRes = levelset.getSize();
	Vec3 mult = toVec3(gridRes) / toVec3(mesh.getParent()->getGridSize());
	
	// prepare center values
	vector<Vec3> center(mesh.numTris());
	for(size_t i=0; i<mesh.numTris(); i++)
		center[i] = mesh.getFaceCenter(i) * mult;

	// prepare grid    
	const int numCells = gridRes.x * gridRes.y * gridRes.z;
	CArray<Real> gridDev(numCells);
	for (int i=0; i<numCells; i++) 
		gridDev[i] = -cutoff;
	gridDev.upload();
	
	// 1. count sources per cell
	CArray<int> srcPerCell(numCells);
	for (size_t i=0; i<center.size(); i++) {
		int cell = _cIndex(center[i], gridRes);
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
	
	// 3. reorder nodes
	CVec3Array reorderPos(center.size());
	CVec3Array reorderNormal(center.size());
	{
		vector<int> curSrcCell(numCells);
		for (int i=0; i<(int)center.size(); i++) {
			int cell = _cIndex(center[i], gridRes);
			if (cell < 0) continue;
			int idx = srcCellStart[cell] + curSrcCell[cell];
			reorderPos.set(idx, center[i]);
			reorderNormal.set(idx, mesh.getFaceNormal(i));
			curSrcCell[cell]++;
		}
	}
	reorderPos.upload();
	reorderNormal.upload();
	
	// construct parameters
	float safeRadius = cutoff + sqrt(3.0)*0.5;
	float safeRadius2 = safeRadius*safeRadius;
	float cutoff2 = cutoff*cutoff;
	float isigma2 = 1.0/(sigma*sigma);
	int intRadius = (int)(cutoff+0.5);
	
	dim3 blocksize(SDFBlockSize, SDFBlockSize, SDFBlockSize);
	dim3 blocks((gridRes.x-1)/SDFBlockSize+1, (gridRes.y-1)/SDFBlockSize+1, (gridRes.z-1)/SDFBlockSize+1);
	SDFKernel<<<blocks, blocksize>>>(srcCellStart.data(), srcPerCell.data(), reorderPos.data(), reorderNormal.data(), gridDev.data(), 
							  make_int3(gridRes.x, gridRes.y, gridRes.z), intRadius, safeRadius2, cutoff2, isigma2);

	gridDev.download();
	for (int i=0;i<numCells; i++)
		levelset[i] = gridDev[i];
	
	// floodfill outside
	stack<Vec3i> outside;
	FOR_IJK(levelset) {
		if (levelset(i,j,k) >= cutoff-1.0f) 
			outside.push(Vec3i(i,j,k));
	}
	while(!outside.empty()) {
		Vec3i c = outside.top();
		outside.pop();
		levelset(c) = cutoff;
		if (c.x > 0 && levelset(c.x-1, c.y, c.z) < 0) outside.push(Vec3i(c.x-1,c.y,c.z));
		if (c.y > 0 && levelset(c.x, c.y-1, c.z) < 0) outside.push(Vec3i(c.x,c.y-1,c.z));
		if (c.z > 0 && levelset(c.x, c.y, c.z-1) < 0) outside.push(Vec3i(c.x,c.y,c.z-1));
		if (c.x < levelset.getSizeX()-1 && levelset(c.x+1, c.y, c.z) < 0) outside.push(Vec3i(c.x+1,c.y,c.z));
		if (c.y < levelset.getSizeY()-1 && levelset(c.x, c.y+1, c.z) < 0) outside.push(Vec3i(c.x,c.y+1,c.z));
		if (c.z < levelset.getSizeZ()-1 && levelset(c.x, c.y, c.z+1) < 0) outside.push(Vec3i(c.x,c.y,c.z+1));
	};
}

} // namespace