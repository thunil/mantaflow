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
	mAccuracy(VECTOR_EPSILON),
	mNumPreSmooth(1),
	mNumPostSmooth(1)
{
	// Create level 0 (=original grid)
	mSizeX.push_back(sizeRef.getSizeX());
	mSizeY.push_back(sizeRef.getSizeY());
	mSizeZ.push_back(sizeRef.getSizeZ());
	int n = mSizeX[0] * mSizeY[0] * mSizeZ[0];

	mA.push_back(std::vector<Real>(n * 14));
	mx.push_back(std::vector<Real>(n));
	mb.push_back(std::vector<Real>(n));
	mr.push_back(std::vector<Real>(n));

	debMsg("GridMg::GridMg level 0: "<<mSizeX[0]<<" x " << mSizeY[0] << " x " << mSizeZ[0] << " x ", 1);

	// Create coarse levels >0
	for (int l=1; ;l++)
	{
		break;
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
//
// y     | z = -1    z = 0      z = 1
// ^     | 6  7  8,  15 16 17,  24 25 26
// |     | 3  4  5,  12 13 14,  21 22 23
// o-> x | 0  1  2,   9 10 11,  18 19 20


void GridMg::setA(FlagGrid& flags, Grid<Real>* A0, Grid<Real>* pAi, Grid<Real>* pAj, Grid<Real>* pAk)
{
	// Copy level 0
	FOR_IDX(*A0)
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
		// loop over coarse grid vertices
		for (int v=0; v<mb[l].size(); v++)
		{
			Vec3i V(v % mSizeX[l], (v % (mSizeX[l]*mSizeY[l])) / mSizeX[l], v / (mSizeX[l]*mSizeY[l]));

			// loop over stencil entries
			for (int s = 0; s<14; s++)
			{
				Vec3i S((s+13)%3-1, ((s+13)%9)/3-1, (s+13)/9-1);
				
				Real sum = Real(0);
				
				for (int r=0; r<27; r++)
				{
					Vec3i R(r%3-1, (r%9)/3-1, r/9-1);
					Real rw = Real(1) / Real(1 << (std::abs(R.x)+std::abs(R.y)+std::abs(R.z)));

					for (int i=0; i<27; i++)
					{
						Vec3i I(i%3-1, (i%9)/3-1, i/9-1);
						Real iw = Real(1) / Real(1 << (std::abs(I.x)+std::abs(I.y)+std::abs(I.z)));

						Vec3i A1 = 2*V+R;
						Vec3i A2 = 2*(V+S)+I;
						int a1 = dot(A1, Vec3i(1,mSizeX[l-1],mSizeX[l-1]*mSizeY[l-1]));
						int a2 = dot(A2, Vec3i(1,mSizeX[l-1],mSizeX[l-1]*mSizeY[l-1]));

						if (A1.x>=0 && A1.x<mSizeX[l-1] && A2.x>=0 && A2.x<mSizeX[l-1] &&
							A1.y>=0 && A1.y<mSizeY[l-1] && A2.y>=0 && A2.y<mSizeY[l-1] &&
							A1.z>=0 && A1.z<mSizeZ[l-1] && A2.z>=0 && A2.z<mSizeZ[l-1])
						{
							Vec3i d = A2-A1;
							if (d.x>=-1 && d.x<=1 && d.y>=-1 && d.y<=1 && d.z>=-1 && d.z<=1)
							{
								int stencil = dot(d+1, Vec3i(1,3,9));

								if (stencil >= 13)
								{
									sum += rw * mA[l-1][a1*14 + stencil - 13] * iw;
								}
								else
								{
									sum += rw * mA[l-1][a2*14 + 13 - stencil] * iw;
								}
							}
						}
					}
				}

				mA[l][v*14 + s] = sum;
			}
		}		
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
		mx[0][i] = Real(0);
	}

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

	return res >= mAccuracy;
}


void GridMg::smoothGS(int l)
{
	for (int v=0; v<mb[l].size(); v++) {		
		if (mA[l][v*14 + 0] == Real(0)) continue;

		Vec3i V(v % mSizeX[l], (v % (mSizeX[l]*mSizeY[l])) / mSizeX[l], v / (mSizeX[l]*mSizeY[l])); 	

		Real sum = mb[l][v];

		for (int s=0; s<27; s++) {
			if (s==13) continue;

			Vec3i S(s%3-1, (s%9)/3-1, s/9-1);
			Vec3i N = V + S;
			int n = dot(N, Vec3i(1,mSizeX[l],mSizeX[l]*mSizeY[l]));

			if (N.x>=0 && N.x<mSizeX[l] && N.y>=0 && N.y<mSizeY[l] && N.z>=0 && N.z<mSizeZ[l]) {
				if (s < 13) {
					sum -= mA[l][n*14 + 13 - s] * mx[l][n];
				} else {
					sum -= mA[l][v*14 + s - 13] * mx[l][n];
				}
			}
		}

		mx[l][v] = sum / mA[l][v*14 + 0];
	}
}

void GridMg::calcResidual(int l)
{
	for (int v=0; v<mb[l].size(); v++) {
		if (mA[l][v*14 + 0] == Real(0)) continue;
		
		Vec3i V(v % mSizeX[l], (v % (mSizeX[l]*mSizeY[l])) / mSizeX[l], v / (mSizeX[l]*mSizeY[l])); 	

		Real sum = mb[l][v];

		for (int s=0; s<27; s++) {
			Vec3i S(s%3-1, (s%9)/3-1, s/9-1);
			Vec3i N = V + S;
			int n = dot(N, Vec3i(1,mSizeX[l],mSizeX[l]*mSizeY[l]));

			if (N.x>=0 && N.x<mSizeX[l] && N.y>=0 && N.y<mSizeY[l] && N.z>=0 && N.z<mSizeZ[l]) {
				if (s < 13) {
					sum -= mA[l][n*14 + 13 - s] * mx[l][n];
				} else {
					sum -= mA[l][v*14 + s - 13] * mx[l][n];
				}
			}
		}

		mr[l][v] = sum;
	}
}

Real GridMg::calcResidualNorm(int l)
{
	Real res = Real(0);

	for (int v=0; v<mb[l].size(); v++) {
		res += mr[l][v] * mr[l][v];
	}

	return std::sqrt(res);
}

void GridMg::solveCG(int l)
{
	// temporarily replace with full GS solve
	while (true)
	{
		calcResidual(l);
		Real res = calcResidualNorm(l);
		if (res < Real(1E-6)) break;

		for (int i=0; i<5; i++)
		{
			smoothGS(l);
		}
	}
}

void GridMg::restrict(int l_dst, std::vector<Real>& src, std::vector<Real>& dst)
{

}

void GridMg::interpolate(int l_dst, std::vector<Real>& src, std::vector<Real>& dst)
{

}



}; // DDF
