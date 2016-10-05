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
 * Author: Florian Ferstl (florian.ferstl.ff@gmail.com)
 *
 * This is an implementation of the solver developed by Dick et al. [1]
 * without topology awareness (= vertex duplication on coarser levels). This 
 * simplification allows us to use regular grids for all levels of the multigrid
 * hierarchy and works well for moderately complex domains.
 *
 * [1] Solving the Fluid Pressure Poisson Equation Using Multigrid-Evaluation
 *     and Improvements, C. Dick, M. Rogowsky, R. Westermann, IEEE TVCG 2015
 *
 ******************************************************************************/

#ifndef _MULTIGRID_H
#define _MULTIGRID_H

#include "vectorbase.h"
#include "grid.h"

namespace Manta { 

//! Multigrid solver
class GridMg {
	public:
		//! constructor: preallocates most of required memory for multigrid hierarchy
		GridMg(const Vec3i& gridSize);
		~GridMg() {};

		//! update system matrix A from symmetric 7-point stencil
		void setA(Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk);
		
		//! set right-hand side
		void setRhs(Grid<Real>& rhs);
		
		//! perform VCycle iteration
		// - if src is null, then a zero vector is used instead
		// - returns norm of residual after VCylcle
		Real doVCycle(Grid<Real>& dst, Grid<Real>* src = nullptr); 
		
		// access
		void setCoarsestLevelAccuracy(Real accuracy) { mCoarsestLevelAccuracy = accuracy; }
		Real getCoarsestLevelAccuracy() const { return mCoarsestLevelAccuracy; }
		void setSmoothing(int numPreSmooth, int numPostSmooth) { mNumPreSmooth = numPreSmooth; mNumPostSmooth = numPostSmooth; }
		int getNumPreSmooth() { return mNumPreSmooth; }
		int getNumPostSmooth() { return mNumPostSmooth; }

	private:
		Vec3i vecIdx(int   v, int l) { return Vec3i(v%mSize[l].x, (v%(mSize[l].x*mSize[l].y))/mSize[l].x, v/(mSize[l].x*mSize[l].y)); }
		int   linIdx(Vec3i V, int l) { return V.x + V.y*mPitch[l].y + V.z*mPitch[l].z; }
		bool  inGrid(Vec3i V, int l) { return V.x>=0 && V.y>=0 && V.z>=0 && V.x<mSize[l].x && V.y<mSize[l].y && V.z<mSize[l].z; }

		bool isEquationTrivial(int v);

		void genCoarseGrid(int l);
		void genCoraseGridOperator(int l);

		void smoothGS(int l, bool reversedOrder);
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

		int mNumPreSmooth;
		int mNumPostSmooth;
		Real mCoarsestLevelAccuracy;

		std::vector<std::vector<Real>> mA;
		std::vector<std::vector<Real>> mx;
		std::vector<std::vector<Real>> mb;
		std::vector<std::vector<Real>> mr;
		std::vector<std::vector<char>> mActive;
		std::vector<std::vector<Real>> mCGtmp1, mCGtmp2;
		std::vector<Vec3i> mSize, mPitch;
		std::vector<CoarseningPath> mCoarseningPaths0;

		bool mIs3D;
		int mDim;
		int mStencilSize;
		int mStencilSize0;
		Vec3i mStencilMin;
		Vec3i mStencilMax;
}; // GridMg

} // namespace

#endif 
