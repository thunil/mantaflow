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

#ifndef _MULTIGRID_H
#define _MULTIGRID_H

#include "vectorbase.h"
#include "grid.h"
// #include "kernel.h"


namespace Manta { 

//! Multigrid solver
class GridMg {
	public:
		//! constructor
		GridMg(const Vec3i& gridSize);
		~GridMg() {};

		// solving functions
		void setA(Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk);
		void setRhs(Grid<Real>& rhs);
		
		// returns false if finished
		Real doVCycle(Grid<Real>& dst, Grid<Real>* src = nullptr); // if src is null, then a zero vector is used instead
		
		// access
		//Real getIterations() const = 0;
		//Real getResNorm() const = 0;
		void setCoarsestLevelAccuracy(Real accuracy) { mCoarsestLevelAccuracy = accuracy; }
		Real getCoarsestLevelAccuracy() const { return mCoarsestLevelAccuracy; }
		void setSmoothing(int numPreSmooth, int numPostSmooth) { mNumPreSmooth = numPreSmooth; mNumPostSmooth = numPostSmooth; }
		int getNumPreSmooth() { return mNumPreSmooth; }
		int getNumPostSmooth() { return mNumPostSmooth; }

	private:
		Vec3i vecIdx(int   v, int l) { return Vec3i(v%mSize[l].x, (v%(mSize[l].x*mSize[l].y))/mSize[l].x, v/(mSize[l].x*mSize[l].y)); }
		int   linIdx(Vec3i V, int l) { return V.x + V.y*mPitch[l].y + V.z*mPitch[l].z; }
		bool  inGrid(Vec3i V, int l) { return V.x>=0 && V.y>=0 && V.z>=0 && V.x<mSize[l].x && V.y<mSize[l].y && V.z<mSize[l].z; }

		void genCoarseGrid(int l);
		void genCoraseGridOperator(int l);

		void smoothGS(int l);
		void calcResidual(int l);
		Real calcResidualNorm(int l);
		void solveCG(int l);

		void restrict(int l_dst, std::vector<Real>& src, std::vector<Real>& dst);
		void interpolate(int l_dst, std::vector<Real>& src, std::vector<Real>& dst);

	private:
		struct CoarseningPath {
			Vec3i U, W, N;
			int sc, sf;
			Real rw, iw;
			bool inUStencil;
		};

		//! accuracy of solver (max. residuum)
		int mNumPreSmooth;
		int mNumPostSmooth;
		Real mCoarsestLevelAccuracy;

		// A has a 7-point stencil on level 0, and a full 27-point stencil on levels >0
		std::vector<std::vector<Real>> mA; // A[level][vertex*14 + stencilentry]
		std::vector<std::vector<Real>> mx; // x[level][vertex]
		std::vector<std::vector<Real>> mb; // b[level][vertex]
		std::vector<std::vector<Real>> mr; // residual[level][vertex]
		std::vector<std::vector<char>> mActive; // active[level][vertex]
		std::vector<Vec3i> mSize, mPitch;
		std::vector<CoarseningPath> mCoarseningPaths0;

		bool mIs3D;
		int mDim;
		int mStencilSize;
		int mStencilSize0;
		Vec3i mStencilMin;
		Vec3i mStencilMax;
}; // GridCg

} // namespace

#endif 
