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

#include "multigrid.h"

using namespace std;
namespace Manta {

GridMg::GridMg(const Grid<Real>& sizeRef)
  : mA(),
	mAccuracy(VECTOR_EPSILON)
{
	// Create level 0 (=original grid)
	mSizeX.push_back(sizeRef.getSizeX());
	mSizeY.push_back(sizeRef.getSizeY());
	mSizeZ.push_back(sizeRef.getSizeZ());
	int n = mSizeX[0] * mSizeY[0] * mSizeZ[0];

	mA.push_back(std::vector<Real>(n * 14));
	mx.push_back(std::vector<Real>());
	mb.push_back(std::vector<Real>(n));
	mr.push_back(std::vector<Real>(n));

	debMsg("GridMg::GridMg level 0: "<<mSizeX[0]<<" x " << mSizeY[0] << " x " << mSizeZ[0] << " x ", 1);

	// Create coarse levels >0
	for (int l=1; ;l++)
	{
		if (mSizeX[l-1] <= 5 && mSizeY[l-1] <= 5 && mSizeZ[l-1] <= 5)
			break;

		mSizeX.push_back((mSizeX[l-1] + 2) / 2);
		mSizeY.push_back((mSizeY[l-1] + 2) / 2);
		mSizeZ.push_back((mSizeZ[l-1] + 2) / 2);
		int n = mSizeX[l] * mSizeY[l] * mSizeZ[l];

		mA.push_back(std::vector<Real>(n * 14));
		mx.push_back(std::vector<Real>(n));
		mb.push_back(std::vector<Real>(n));
		mr.push_back(std::vector<Real>(n));
		
		debMsg("GridMg::GridMg level "<<l<<": " << mSizeX[l] << " x " << mSizeY[l] << " x " << mSizeZ[l] << " x ", 1);
	}
}


// 27-Point stencil indices
// y     | z = -1    z = 0      z = 1
// ^     | x  x  x,   2  3  4,  11 12 13
// |     | x  x  x,   x  0  1,   8  9 10
// o-> x | x  x  x,   x  x  x,   5  6  7


void GridMg::setA(FlagGrid& flags, Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk)
{
	// Copy level 0
	FOR_IDX(rhs)
	{
		for (int i=0; i<14; i++) mA[0][idx * 14 + i] = Real(0);

		mA[0][idx * 14 +  0] = (* A0)[idx];
		mA[0][idx * 14 +  1] = (*pAi)[idx];
		mA[0][idx * 14 +  3] = (*pAj)[idx];
		mA[0][idx * 14 +  9] = (*pAk)[idx];
	}

	// Create coarse operators on levels >0
	for (int l=1; l<mA.size(); l++)
	{
		
	}
}

void GridMg::setRhs(Grid<Real>& rhs)
{
	FOR_IDX(rhs)
	{
		mb[0][idx] = rhs[idx];
	}
}

// returns false if finished
bool GridMg::doVCycle(Grid<Real>& dst)
{
	return true; // not finished
}


}; // DDF
