/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Multigrid solver
 *
 ******************************************************************************/

// TODO
// - 2D specialization
// - active vertex lists
// - parallelization
// - finest level optimization
// - analyze performance
// - accuracy parameters configuration (coarsest CG)

#include "multigrid.h"

#define FOR_LVL(IDX,LVL) \
	for(int IDX=0, IDX##total=mSize[LVL].x*mSize[LVL].y*mSize[LVL].z; IDX<IDX##total; IDX++)

#define FOR_VEC_LVL(VEC,LVL) Vec3i VEC; \
	for(VEC.z=0; VEC.z<mSize[LVL].z; VEC.z++) \
	for(VEC.y=0; VEC.y<mSize[LVL].y; VEC.y++) \
	for(VEC.x=0; VEC.x<mSize[LVL].x; VEC.x++)

#define FOR_VEC_MINMAX(VEC,MIN,MAX) Vec3i VEC; \
	const Vec3i VEC##__min = (MIN), VEC##__max = (MAX); \
	for(VEC.z=VEC##__min.z; VEC.z<=VEC##__max.z; VEC.z++) \
	for(VEC.y=VEC##__min.y; VEC.y<=VEC##__max.y; VEC.y++) \
	for(VEC.x=VEC##__min.x; VEC.x<=VEC##__max.x; VEC.x++)

#define FOR_VECLIN_MINMAX(VEC,LIN,MIN,MAX) Vec3i VEC; int LIN = 0; \
	const Vec3i VEC##__min = (MIN), VEC##__max = (MAX); \
	for(VEC.z=VEC##__min.z; VEC.z<=VEC##__max.z; VEC.z++) \
	for(VEC.y=VEC##__min.y; VEC.y<=VEC##__max.y; VEC.y++) \
	for(VEC.x=VEC##__min.x; VEC.x<=VEC##__max.x; VEC.x++, LIN++)


using namespace std;
namespace Manta 
{

// ----------------------------------------------------------------------------
// Efficient min heap for <ID, key> pairs with 0<=ID<N and 0<=key<K
// (IDs are stored in K buckets, where each bucket is a list of IDs).
// - if K<<N, all ops are O(1) on avg (worst case O(K)).
// - memory usage O(K+N): (K+N) * 3 * sizeof(int).
class NKMinHeap
{
private:
	struct Entry {
		int key, prev, next;
		Entry() : key(-1), prev(-1), next(-1) {}
	};

	int mN, mK, mSize, mMinKey;

	// Double linked lists of IDs, one for each bucket/key.
	// The first K entries are the buckets' head pointers,
	// and the last N entries correspond to the IDs.
	std::vector<Entry> mEntries; 
	
public:
	NKMinHeap(int N, int K) : mN(N), mK(K), mSize(0), mMinKey(-1), mEntries(N+K) {}

	int size() { return mSize; }
	int getKey(int ID) { return mEntries[mK+ID].key; }
	
	// Insert, decrease or increase key (or delete by setting key to -1)
	void setKey(int ID, int key);

	// peek min key (returns ID/key pair)
	std::pair<int,int> peekMin();

	// pop min key (returns ID/key pair)
	std::pair<int,int> popMin();
	
	void print(); // for debugging
};

void NKMinHeap::setKey(int ID, int key)
{
	assertMsg(0 <=ID  && ID <mN, "NKMinHeap::setKey: ID out of range");
	assertMsg(-1<=key && key<mK, "NKMinHeap::setKey: key out of range");

	const int kid = mK + ID;

	if (mEntries[kid].key == key) return; // nothing changes

	// remove from old key-list if ID existed previously
	if (mEntries[kid].key != -1)
	{
		int pred = mEntries[kid].prev;
		int succ = mEntries[kid].next; // can be -1

		mEntries[pred].next = succ;
		if (succ != -1) mEntries[succ].prev = pred;

		// if removed key was minimum key, mMinKey may need to be updated
		int removedKey = mEntries[kid].key;
		if (removedKey==mMinKey)
		{
			if (mSize==1) { mMinKey = -1; }
			else {
				for (; mMinKey<mK; mMinKey++) {
					if (mEntries[mMinKey].next != -1) break;
				}
			}
		}

		mSize--;
	}

	// set new key of ID
	mEntries[kid].key = key;

	if (key==-1) {
		// finished if key was set to -1
		mEntries[kid].next = mEntries[kid].prev = -1;
		return;
	}

	// add key
	mSize++;
	if (mMinKey == -1) mMinKey = key;
	else mMinKey = std::min(mMinKey, key);

	// insert into new key-list (headed by mEntries[key])
	int tmp = mEntries[key].next;

	mEntries[key].next = kid;
	mEntries[kid].prev = key;

	mEntries[kid].next = tmp;
	if (tmp != -1) mEntries[tmp].prev = kid;
}

std::pair<int,int> NKMinHeap::peekMin()
{
	if (mSize==0) return std::pair<int,int>(-1,-1); // error
		
	const int ID = mEntries[mMinKey].next - mK;
	return std::pair<int,int>(ID, mMinKey);
}

std::pair<int,int> NKMinHeap::popMin()
{
	if (mSize==0) return std::pair<int,int>(-1,-1); // error

	const int kid = mEntries[mMinKey].next;
	const int ID = kid - mK;
	const int key = mMinKey;

	// remove from key-list
	int pred = mEntries[kid].prev;
	int succ = mEntries[kid].next; // can be -1

	mEntries[pred].next = succ;
	if (succ != -1) mEntries[succ].prev = pred;

	// remove entry
	mEntries[kid] = Entry();	
	mSize--;

	// update mMinKey
	if (mSize==0) { mMinKey = -1; }
	else {
		for (; mMinKey<mK; mMinKey++) {
			if (mEntries[mMinKey].next != -1) break;
		}
	}
		
	// return result
	return std::pair<int,int>(ID, key);
}

void NKMinHeap::print()
{
	std::cout << "Size: "<<mSize<<", MinKey: "<<mMinKey<< std::endl;
	for (int key=0; key<mK; key++) {
		if (mEntries[key].next != -1) {
			std::cout << "Key "<<key<<": ";
			int kid = mEntries[key].next;
			while (kid != -1) {
				std::cout << kid-mK<<" ";
				kid = mEntries[kid].next;
			}
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

// ----------------------------------------------------------------------------
// GridMg methods
//
// Illustration of 27-point stencil indices
// y     | z = -1    z = 0      z = 1
// ^     | 6  7  8,  15 16 17,  24 25 26
// |     | 3  4  5,  12 13 14,  21 22 23
// o-> x | 0  1  2,   9 10 11,  18 19 20
//
// Symmetric storage with only 14 entries per vertex
// y     | z = -1    z = 0      z = 1
// ^     | -  -  -,   2  3  4,  11 12 13
// |     | -  -  -,   -  0  1,   8  9 10
// o-> x | -  -  -,   -  -  -,   5  6  7

GridMg::GridMg(const Vec3i& gridSize)
  : mA(),
	mAccuracy(VECTOR_EPSILON),
	mNumPreSmooth(1),
	mNumPostSmooth(1)
{
	// Create level 0 (=original grid)
	mSize.push_back(gridSize);
	mPitch.push_back(Vec3i(1, mSize.back().x, mSize.back().x*mSize.back().y));
	int n = mSize.back().x * mSize.back().y * mSize.back().z;

	mA.push_back(std::vector<Real>(n * 14));
	mx.push_back(std::vector<Real>(n));
	mb.push_back(std::vector<Real>(n));
	mr.push_back(std::vector<Real>(n));
	mActive.push_back(std::vector<char>(n));

	debMsg("GridMg::GridMg level 0: "<<mSize[0].x<<" x " << mSize[0].y << " x " << mSize[0].z << " x ", 1);

	// Create coarse levels >0
	for (int l=1; l<=100; l++)
	{
		if (mSize[l-1].x <= 5 && mSize[l-1].y <= 5 && mSize[l-1].z <= 5)
			break;

		mSize.push_back((mSize[l-1] + 2) / 2);
		mPitch.push_back(Vec3i(1, mSize.back().x, mSize.back().x*mSize.back().y));
		int n = mSize.back().x * mSize.back().y * mSize.back().z;

		mA.push_back(std::vector<Real>(n * 14));
		mx.push_back(std::vector<Real>(n));
		mb.push_back(std::vector<Real>(n));
		mr.push_back(std::vector<Real>(n));
		mActive.push_back(std::vector<char>(n));
		
		debMsg("GridMg::GridMg level "<<l<<": " << mSize[l].x << " x " << mSize[l].y << " x " << mSize[l].z << " x ", 1);
	}
}

void GridMg::setA(Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk)
{
	// Copy level 0
	FOR_LVL(v,0) {
		for (int i=0; i<14; i++) { mA[0][v*14 + i] = Real(0); }

		mA[0][v*14 + 0] = (* A0)[v];
		mA[0][v*14 + 1] = (*pAi)[v];
		mA[0][v*14 + 3] = (*pAj)[v];
		mA[0][v*14 + 9] = (*pAk)[v];
			
		mActive[0][v] = char(mA[0][v*14 + 0] != Real(0));
	}

	// Create coarse grids and operators on levels >0
	for (int l=1; l<mA.size(); l++) {
		genCoarseGrid(l);	
		genCoraseGridOperator(l);	
	}
}

void GridMg::setRhs(Grid<Real>& rhs)
{
	FOR_IDX(rhs)
	{
		mb[0][idx] = rhs[idx];
	}
}

bool GridMg::doVCycle(Grid<Real>& dst)
{
	const int maxLevel = mA.size() - 1;

	for (int i=0; i<mx[0].size(); i++) {
		mx[0][i] = dst[i];
	}

	// Next two lines are debug code, remove later
	calcResidual(0);
	Real resOld = calcResidualNorm(0);

	for (int l=0; l<maxLevel; l++)
	{
		for (int i=0; i<mNumPreSmooth; i++) {
			smoothGS(l);
		}
		calcResidual(l);
		restrict(l+1, mr[l], mb[l+1]);

		for (int i=0; i<mx[l+1].size(); i++) {
			mx[l+1][i] = Real(0);
		}
	}

	solveCG(maxLevel);

	for (int l=maxLevel-1; l>=0; l--)
	{
		interpolate(l, mx[l+1], mr[l]);
		for (int i=0; i<mx[l].size(); i++) {
			mx[l][i] += mr[l][i];
		}

		for (int i=0; i<mNumPostSmooth; i++) {
			smoothGS(l);
		}
	}

	calcResidual(0);
	Real res = calcResidualNorm(0);

	for (int i=0; i<mx[0].size(); i++) {
		dst[i] = mx[0][i];
	}

	debMsg("VCycle Residual: "<<resOld<<" -> "<<res, 1);

	return res >= mAccuracy;
}

// Determine active cells on coarse level l from active cells on fine level l-1
// This is an implementation of the algorithm described as part of Section 3.3 in
//     Solving the Fluid Pressure Poisson Equation Using Multigrid-Evaluation
//     and Improvements, C. Dick, M. Rogowsky, R. Westermann, IEEE TVCG 2015
// to ensure a full-rank interpolation operator.
void GridMg::genCoarseGrid(int l)
{
	//    AF_Free: unused/untouched vertices
	//    AF_Zero: vertices selected for coarser level
	// AF_Removed: vertices removed from coarser level
	enum activeFlags : char {AF_Removed = 0, AF_Zero = 1, AF_Free = 2, };

	// initialize all coarse vertices with 'free'
	FOR_LVL(v,l) { mActive[l][v] = AF_Free; }

	// initialize min heap of (ID: fine grid vertex, key: #free interpolation vertices) pairs
	NKMinHeap heap(mb[l-1].size(), 9); // max 8 free interpolation vertices
		
	FOR_LVL(v,l-1) {
		if (mActive[l-1][v]) {
			Vec3i V = vecIdx(v,l-1);
			int fiv = 1 << ((V.x % 2) + (V.y % 2) + (V.z % 2));
			heap.setKey(v, fiv);
		}
	}

	// process fine vertices in heap consecutively, always choosing the vertex with 
	// the currently smallest number of free interpolation vertices
	while (heap.size() > 0)
	{
		int v = heap.popMin().first;
		Vec3i V = vecIdx(v,l-1);


		// loop over associated interpolation vertices of V on coarse level l:
		// the first encountered 'free' vertex is set to 'zero',
		// all remaining 'free' vertices are set to 'removed'.
		bool vdone = false;

		FOR_VEC_MINMAX(I, V/2, (V+1)/2) {
			int i = linIdx(I,l);

			if (mActive[l][i] == AF_Free) {
				if (vdone) {
					mActive[l][i] = AF_Removed; 
				} else {
					mActive[l][i] = AF_Zero; 
					vdone = true;
				}

				// update #free interpolation vertices in heap:
				// loop over all associated restriction vertices of I on fine level l-1
				FOR_VEC_MINMAX(R, vmax(0, I*2-1), vmin(mSize[l-1]-1, I*2+1)) {
					int r = linIdx(R,l-1);
					int key = heap.getKey(r); 

					if      (key > 1) { heap.setKey(r, key-1); } // decrease key of r
					else if (key >-1) { heap.setKey(r, -1); } // removes r from heap
				}
			}
		}
	}

	// set all remaining 'free' vertices to 'removed'
	FOR_LVL(v,l) { if (mActive[l][v] == AF_Free) mActive[l][v] = AF_Removed; }
}

// Calculate A_l on coarse level l from A_{l-1} on fine level l-1 using 
// Galerkin-based coarsening, i.e., compute A_l = R * A_{l-1} * I.
void GridMg::genCoraseGridOperator(int l)
{
	// loop over coarse grid vertices V
	FOR_LVL(v,l) {
		if (!mActive[l][v]) continue;

		for (int i=0; i<14; i++) { mA[l][v*14+i] = Real(0); } // clear stencil

		Vec3i V = vecIdx(v,l);

		// Calculate the stencil of A_l at V by considering all vertex paths of the form:
		// (V) <--restriction-- (U) <--A_{l-1}-- (W) <--interpolation-- (N)
		// V and N are vertices on the coarse grid level l, 
		// U and W are vertices on the fine grid level l-1.

		// loop over restriction vertices U on level l-1 associated with V
		FOR_VEC_MINMAX(U, vmax(0, V*2-1), vmin(mSize[l-1]-1, V*2+1)) {
			int u = linIdx(U,l-1);
			if (!mActive[l-1][u]) continue;

			// restriction weight			
			Real rw = Real(1) / Real(1 << ((U.x % 2) + (U.y % 2) + (U.z % 2))); 

			// loop over all stencil neighbors N of V level l that can be reached via restriction to U
			FOR_VEC_MINMAX(N, (U-1)/2, vmin(mSize[l]-1, (U+2)/2)) {
				int n = linIdx(N,l);
				if (!mActive[l][n]) continue;
							
				// stencil entry at V associated to N (coarse grid level l)
				Vec3i SC = N-V+1;
				int sc = SC.x + 3*SC.y + 9*SC.z;
				if (sc < 13) continue;

				// loop over all vertices W which are in the stencil of A_{l-1} at U 
				// and which interpolate from N
				FOR_VEC_MINMAX(W, vmax(           0, vmax(U-1,N*2-1)),
				                  vmin(mSize[l-1]-1, vmin(U+1,N*2+1))) {
					int w = linIdx(W,l-1);
					if (!mActive[l-1][w]) continue;

					// stencil entry at U associated to W (fine grid level l-1)
					Vec3i SF = W-U+1;
					int sf = SF.x + 3*SF.y + 9*SF.z;

					Real iw = Real(1) / Real(1 << ((W.x % 2) + (W.y % 2) + (W.z % 2))); // interpolation weight

					if (sf < 13) {
						mA[l][v*14 + sc-13] += rw * mA[l-1][w*14 + 13-sf] *iw;
					} else {
						mA[l][v*14 + sc-13] += rw * mA[l-1][u*14 + sf-13] *iw;
					}
				}
			}
		}
	}		
}

void GridMg::smoothGS(int l)
{
	FOR_LVL(v,l) {
		if (!mActive[l][v]) continue;

		Vec3i V = vecIdx(v,l);

		Real sum = mb[l][v];

		FOR_VECLIN_MINMAX(S, s, -1, 1) {
			if (s==13) continue;

			Vec3i N = V + S;
			int n = linIdx(N,l);

			if (inGrid(N,l) && mActive[l][n]) {
				if (s < 13) {
					sum -= mA[l][n*14 + 13-s] * mx[l][n];
				} else {
					sum -= mA[l][v*14 + s-13] * mx[l][n];
				}
			}
		}

		mx[l][v] = sum / mA[l][v*14 + 0];
	}
}

void GridMg::calcResidual(int l)
{
	FOR_LVL(v,l) {
		if (!mActive[l][v]) continue;
		
		Vec3i V = vecIdx(v,l);

		Real sum = mb[l][v];

		FOR_VECLIN_MINMAX(S, s, -1, 1) {
			Vec3i N = V + S;
			int n = linIdx(N,l);

			if (inGrid(N,l) && mActive[l][n]) {
				if (s < 13) {
					sum -= mA[l][n*14 + 13-s] * mx[l][n];
				} else {
					sum -= mA[l][v*14 + s-13] * mx[l][n];
				}
			}
		}

		mr[l][v] = sum;
	}
}

Real GridMg::calcResidualNorm(int l)
{
	Real res = Real(0);

	FOR_LVL(v,l) {
		if (!mActive[l][v]) continue;

		res += mr[l][v] * mr[l][v];
	}

	return std::sqrt(res);
}

// Standard conjugate gradients with Jacobi preconditioner
void GridMg::solveCG(int l)
{
	// TODO: preallocate
	std::vector<Real> z(mb[l].size());
	std::vector<Real> p(mb[l].size());

	std::vector<Real>& x = mx[l];
	std::vector<Real>& r = mr[l];

	// Initialization:
	Real alphaTop = Real(0);

	FOR_LVL(v,l) {
		if (!mActive[l][v]) continue;
		
		Vec3i V = vecIdx(v,l);

		Real sum = mb[l][v];

		FOR_VECLIN_MINMAX(S, s, -1, 1) {
			Vec3i N = V + S;
			int n = linIdx(N,l);

			if (inGrid(N,l) && mActive[l][n]) {
				if (s < 13) {
					sum -= mA[l][n*14 + 13-s] * x[n];
				} else {
					sum -= mA[l][v*14 + s-13] * x[n];
				}
			}
		}

		r[v] = sum;
		z[v] = r[v] / mA[l][v*14 + 0];
		p[v] = z[v];
		alphaTop += r[v] * z[v];
	}

	int iter = 0;
	const int maxIter = 10000;
	Real residual = Real(-1);

	// CG iterations
	for (; iter<maxIter; iter++)
	{
		Real alphaBot = Real(0);

		FOR_LVL(v,l) {
			if (!mActive[l][v]) continue;
		
			Vec3i V = vecIdx(v,l);

			z[v] = Real(0);

			FOR_VECLIN_MINMAX(S, s, -1, 1) {
				Vec3i N = V + S;
				int n = linIdx(N,l);

				if (inGrid(N,l) && mActive[l][n]) {
					if (s < 13) {
						z[v] += mA[l][n*14 + 13-s] * p[n];
					} else {
						z[v] += mA[l][v*14 + s-13] * p[n];
					}
				}
			}

			alphaBot += p[v] * z[v];
		}

		Real alpha = alphaTop / alphaBot;
		
		Real alphaTopNew = Real(0);
		residual = Real(0);

		FOR_LVL(v,l) {
			if (!mActive[l][v]) continue;
		
			x[v] += alpha * p[v];
			r[v] -= alpha * z[v];
			residual += r[v] * r[v];
			z[v] = r[v] / mA[l][v*14 + 0];
			alphaTopNew += r[v] * z[v];
		}

		residual = std::sqrt(residual);

		if (residual < 1E-8) break;

		Real beta = alphaTopNew / alphaTop;
		alphaTop = alphaTopNew;

		FOR_LVL(v,l) {
			p[v] = z[v] + beta * p[v];
		}
	}

	if (iter == maxIter) { debMsg("GridMg::solveCG Warning: Reached maximum number of CG iterations", 1); }
	else { debMsg("GridMg::solveCG Info: Reached residual "<<residual<<" in "<<iter<<" iterations", 1); }
}

void GridMg::restrict(int l_dst, std::vector<Real>& src, std::vector<Real>& dst)
{
	const int l_src = l_dst - 1;
	
	FOR_LVL(v,l_dst) {
		if (!mActive[l_dst][v]) continue;

		// Coarse grid vertex
		Vec3i V = vecIdx(v,l_dst);
		
		Real sum = Real(0);

		FOR_VEC_MINMAX(R, vmax(0, V*2-1), vmin(mSize[l_src]-1, V*2+1)) {
			int r = linIdx(R,l_src);
			if (!mActive[l_src][r]) continue;

			// restriction weight			
			Real rw = Real(1) / Real(1 << ((R.x % 2) + (R.y % 2) + (R.z % 2))); 
			
			sum += rw * src[r]; 
		}

		dst[v] = sum;
	}
}

void GridMg::interpolate(int l_dst, std::vector<Real>& src, std::vector<Real>& dst)
{
	const int l_src = l_dst + 1;

	FOR_LVL(v,l_dst) {
		if (!mActive[l_dst][v]) continue;
		
		Vec3i V = vecIdx(v,l_dst);

		Real sum = Real(0);

		FOR_VEC_MINMAX(I, V/2, (V+1)/2) {
			int i = linIdx(I,l_src);
			if (mActive[l_src][i]) sum += src[i]; 
		}

		// restriction weight			
		Real iw = Real(1) / Real(1 << ((V.x % 2) + (V.y % 2) + (V.z % 2)));

		dst[v] = iw * sum;
	}
}



}; // DDF
