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
	mSize.push_back(sizeRef.getSize());
	mPitch.push_back(Vec3i(1, mSize.back().x, mSize.back().x*mSize.back().y));
	int n = mSize.back().x * mSize.back().y * mSize.back().z;

	mA.push_back(std::vector<Real>(n * 14));
	mx.push_back(std::vector<Real>(n));
	mb.push_back(std::vector<Real>(n));
	mr.push_back(std::vector<Real>(n));
	mActive.push_back(std::vector<bool>(n));

	debMsg("GridMg::GridMg level 0: "<<mSize[0].x<<" x " << mSize[0].y << " x " << mSize[0].z << " x ", 1);

	// Create coarse levels >0
	for (int l=1; l<=1; l++)
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
		mActive.push_back(std::vector<bool>(n));
		
		debMsg("GridMg::GridMg level "<<l<<": " << mSize[l].x << " x " << mSize[l].y << " x " << mSize[l].z << " x ", 1);
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
			
		mActive[0][idx] = (mA[0][idx*14 + 0] != Real(0));
	}

	// Create coarse operators on levels >0
	for (int l=1; l<mA.size(); l++)
	{
		// loop over coarse grid vertices
		for (int v=0; v<mb[l].size(); v++)
		{
			Vec3i V(v % mSize[l].x, (v % (mSize[l].x*mSize[l].y)) / mSize[l].x, v / (mSize[l].x*mSize[l].y));

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
						int a1 = dot(A1, mPitch[l-1]);
						int a2 = dot(A2, mPitch[l-1]);

						if (A1.x>=0 && A1.x<mSize[l-1].x && A2.x>=0 && A2.x<mSize[l-1].x &&
							A1.y>=0 && A1.y<mSize[l-1].y && A2.y>=0 && A2.y<mSize[l-1].y &&
							A1.z>=0 && A1.z<mSize[l-1].z && A2.z>=0 && A2.z<mSize[l-1].z)
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

			mActive[l][v] = (mA[l][v*14 + 0] != Real(0));
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
		mx[0][i] = dst[i];
	}

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

	return res >= mAccuracy;
}


void GridMg::smoothGS(int l)
{
	for (int v=0; v<mb[l].size(); v++) {		
		if (!mActive[l][v]) continue;

		Vec3i V(v % mSize[l].x, (v % (mSize[l].x*mSize[l].y)) / mSize[l].x, v / (mSize[l].x*mSize[l].y)); 	

		Real sum = mb[l][v];

		for (int s=0; s<27; s++) {
			if (s==13) continue;

			Vec3i S(s%3-1, (s%9)/3-1, s/9-1);
			Vec3i N = V + S;
			int n = dot(N, mPitch[l]);

			if (N.x>=0 && N.x<mSize[l].x && N.y>=0 && N.y<mSize[l].y && N.z>=0 && N.z<mSize[l].z) {
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
		if (!mActive[l][v]) continue;
		
		Vec3i V(v % mSize[l].x, (v % (mSize[l].x*mSize[l].y)) / mSize[l].x, v / (mSize[l].x*mSize[l].y)); 	

		Real sum = mb[l][v];

		for (int s=0; s<27; s++) {
			Vec3i S(s%3-1, (s%9)/3-1, s/9-1);
			Vec3i N = V + S;
			int n = dot(N, mPitch[l]);

			if (N.x>=0 && N.x<mSize[l].x && N.y>=0 && N.y<mSize[l].y && N.z>=0 && N.z<mSize[l].z) {
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
		if (!mActive[l][v]) continue;

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
	const int l_src = l_dst - 1;

	for (int v=0; v<mb[l_dst].size(); v++) {
		if (!mActive[l_dst][v]) continue;

		// Coarse grid vertex
		Vec3i V(v % mSize[l_dst].x, (v % (mSize[l_dst].x*mSize[l_dst].y)) / mSize[l_dst].x, v / (mSize[l_dst].x*mSize[l_dst].y)); 	
		Vec3i Vfine = V*2; // coordinates on fine grid


		// Box of fine grid vertices to restrict from
		Vec3i RMin = Vfine - 1;
		Vec3i RMax = Vfine + 1;
		for (int d=0; d<3; d++) {
			if (RMin[d] < 0) { RMin[d] = 0; } 
			if (RMax[d] >= mSize[l_src][d]) { RMax[d] = mSize[l_src][d] - 1; }
		}

		Real sum = Real(0);

		Vec3i R;
		for (R.z=RMin.z; R.z<=RMax.z; R.z++)
		for (R.y=RMin.y; R.y<=RMax.y; R.y++)
		for (R.x=RMin.x; R.x<=RMax.x; R.x++)
		{
			int r = dot(R,mPitch[l_src]);
			if (!mActive[l_src][r]) continue;
			Vec3i D = (R - Vfine);
			Real w = Real(1) / Real(1 << (std::abs(D.x)+std::abs(D.y)+std::abs(D.z)));
			sum += w * src[r]; 
		}

		dst[v] = sum;
	}
}

void GridMg::interpolate(int l_dst, std::vector<Real>& src, std::vector<Real>& dst)
{
	const int l_src = l_dst + 1;

	for (int v=0; v<mb[l_dst].size(); v++) {
		if (!mActive[l_dst][v]) continue;
		
		Vec3i V(v % mSize[l_dst].x, (v % (mSize[l_dst].x*mSize[l_dst].y)) / mSize[l_dst].x, v / (mSize[l_dst].x*mSize[l_dst].y)); 	

		Vec3i IMin = V / 2;
		Vec3i IMax = (V+1) / 2;

		Real sum = Real(0);

		Vec3i I;
		for (I.z=IMin.z; I.z<=IMax.z; I.z++)
		for (I.y=IMin.y; I.y<=IMax.y; I.y++)
		for (I.x=IMin.x; I.x<=IMax.x; I.x++)
		{
			int i = dot(I,mPitch[l_src]);
			if (mActive[l_src][i]) sum += src[i]; 
		}

		Real w = Real(1) / Real(1 << dot(IMax-IMin,Vec3i(1)));
		dst[v] = w * sum;
	}
}



}; // DDF
