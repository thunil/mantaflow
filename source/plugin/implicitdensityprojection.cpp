/******************************************************************************
 *
 * MantaFlow fluid solver framework 
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * IDP-FLIP (Implicit-Density-Projection FLIP / APIC)
 * for use with particle data fields
 * Author: Tassilo Kugelstadt
 *
 ******************************************************************************/

#include "particle.h"
#include "grid.h"
#include "commonkernels.h"
#include "randomstream.h"
#include "shapes.h"

using namespace std;
namespace Manta {

//******************************************************************************
// marking of fluid and boundary cells that contain a particle
KERNEL() void knClearFluidFlags(FlagGrid& flags, int dummy = 0) {
	if (flags.isFluid(i, j, k)) {
		flags(i, j, k) = (flags(i, j, k) | FlagGrid::TypeEmpty) & ~FlagGrid::TypeFluid;
	}
}

PYTHON() void markFluidAndBoundaryCells(const BasicParticleSystem& particles, FlagGrid& flags, MACGrid& deltaX, const Grid<Real>& phiObs, const ParticleDataImpl<int>* ptype = NULL, const int exclude = 0) {
	// remove all fluid cells
	knClearFluidFlags(flags, 0);
	deltaX.clear();

	// mark all particles in flaggrid as fluid
	for (IndexInt idx = 0; idx< particles.size(); idx++) {
		if (!particles.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) continue;
		Vec3 pos = particles.getPos(idx);
		Vec3i p = toVec3i(pos);
		if (flags.isInBounds(p) && flags.isEmpty(p))
			flags(p) = (flags(p) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
		else if (flags.isInBounds(p) && flags.isObstacle(p))	//particle is in a boundary cell. Set displacements on the boundary cell so that the particle gets pushed out.
		{
			Real dist = phiObs.getInterpolated(pos);	//distance to the closest point on the surface
			if(dist > 0) continue;
			//compute gradient direction by using central FD
			Vec3 dir;
			Real eps = 1.0e-3;
			dir.x = (phiObs.getInterpolated(pos + Vec3(eps, 0, 0)) - phiObs.getInterpolated(pos + Vec3(-eps, 0, 0))) / ((Real)2.0 * eps);
			dir.y = (phiObs.getInterpolated(pos + Vec3(0, eps, 0)) - phiObs.getInterpolated(pos + Vec3(0, -eps, 0))) / ((Real)2.0 * eps);
			if (flags.is3D())
				dir.z = (phiObs.getInterpolated(pos + Vec3(0, 0, eps)) - phiObs.getInterpolated(pos + Vec3(0, 0, -eps))) / ((Real)2.0 * eps);
				
			if (dist < -1.0) dist = -1.0;	//clamp displacement
			dir = -(dist + 1.0e-2) * dir;

			//set displacement on all adjacent faces of the MAC grid. Check if the displacement is already set by another particle. 
			//If it is already larger than the displacement of the current particle then we don't set it.
			int i = p.x; int j = p.y; int k = p.z;
			if (fabs(dir.x) > fabs(deltaX(i, j, k).x)) deltaX(i, j, k).x = dir.x;
			if (fabs(dir.x) > fabs(deltaX(i + 1, j, k).x)) deltaX(i + 1, j, k).x = dir.x;
			if (fabs(dir.y) > fabs(deltaX(i, j, k).y)) deltaX(i, j, k).y = dir.y;
			if (fabs(dir.y) > fabs(deltaX(i, j + 1, k).y)) deltaX(i, j + 1, k).y = dir.y;
			if (flags.is3D()) {
				if (fabs(dir.z) > fabs(deltaX(i, j, k).z)) deltaX(i, j, k).z = dir.z;
				if (fabs(dir.z) > fabs(deltaX(i, j, k + 1).z)) deltaX(i, j, k + 1).z = dir.z;
			}
		}
	}
}

//******************************************************************************
// mapping of particle masses to the grid
KERNEL(pts, single) template<class T>
void knMapLinear(const BasicParticleSystem& p, const FlagGrid& flags, const Grid<T>& target, Grid<Real>& gtmp,
        const ParticleDataImpl<T>& psource )
{
	unusedParameter(flags);
	if (!p.isActive(idx)) return;
	target.setInterpolated( p[idx].pos, psource[idx], gtmp );
} 

KERNEL(pts) template<class T>
void knMapFromGrid(const BasicParticleSystem& p, const Grid<T>& gsrc, ParticleDataImpl<T>& target)
{
	if (!p.isActive(idx)) return;
	target[idx] = gsrc.getInterpolated(p[idx].pos);
}

// compute density and set boundary conditions for obstacles
KERNEL(bnd=0) template<class T>
void knComputeDensity(Grid<T>& density, const FlagGrid& flagsTmp, FlagGrid& flags, const MACGrid& deltaX, Real dt, Real mass, bool noDensityClamping = false)
{
	IndexInt idx = flags.index(i, j, k);
	if (flags.isFluid(idx))
	{
		density[idx] = (static_cast<Real>(1.0) - density[idx] * mass);	//grid cell size = 1 in manta, rho_0 = 1.0
						
		//add deltaX from the solid cells that contain particles like in [Bridson - Fluid Simulation for Computer Graphics - 2015 - Fig. 5.4]. (for all non-boundary cells deltaX is 0, so no special cases needed)
		density[idx] -= deltaX[idx].x - deltaX(i + 1, j, k).x + deltaX[idx].y - deltaX(i, j + 1, k).y;
		if (flags.is3D()) density[idx] -= deltaX[idx].z - deltaX(i, j, k + 1).z;

		//clamp only if at least one neighbor cell is empty
		bool isSurface = flagsTmp.isEmpty(i - 1, j, k) || flagsTmp.isEmpty(i + 1, j, k) || flagsTmp.isEmpty(i, j - 1, k) || flagsTmp.isEmpty(i, j + 1, k);
		if (flags.is3D()) isSurface = isSurface || (flagsTmp.isEmpty(i, j, k - 1) || flagsTmp.isEmpty(i, j, k + 1));

		//account for the missing density due to particle deficiency in the boundary -> pretend that obstacle cells are uniformly sampled with particles 
		//todo: better solution would be to compute a density grid for the static boundary. This can be done by simply samping the boundary level-set with particles and map them to the density_0 grid.
		//todo: then in the solver the density should not get cleared but reset to density_0. This works for arbitrary samplings and interpolation kernels.
		if(flags.is3D())	//todo: add for 2d, and add arbitrary particle per cell samplings. This only works for 8 particles per cell.
		{					
			Real N[3] = { 0.25, 0.75, 0.25 };

			for (int l = -1; l < 2; l++)
				for (int m = -1; m < 2; m++)
					for (int n = -1; n < 2; n++)
						if(flags.isObstacle(i + l, j + m, k + n) || flags.isEmpty(i + l, j + m, k + n))
							if(l==0 && m==0 || l==0 && k==0 || m==0 && k==0)
								density[idx] -= N[l+1] * N[m + 1] * N[n + 1] * mass * 4.0;		//face neighbors
							else if(l != 0 && m != 0 || l != 0 && k != 0 || m != 0 && k != 0)
								density[idx] -= N[l + 1] * N[m + 1] * N[n + 1] * mass * 2.0;	//edge neighbors
							else
								density[idx] -= N[l + 1] * N[m + 1] * N[n + 1] * mass;			//vertex neighbors
		}
		
		if(isSurface && density[idx] > static_cast<Real>(0.0))
		{
			flags[idx] = FlagGrid::TypeEmpty;
			density[idx] = 0;
		}

		if(!noDensityClamping)
		{
			if (density[idx] < -0.5)
				density[idx] = -0.5;
			
			if (density[idx] > 0.5)
				density[idx] = 0.5;

			density[idx] = density[idx] / dt;
		}
	}
	else  //setting the density to 0 and excluding the cell from the pressure solve is the same boundary condition as in the divegence free or normal pressure solve.
		density[idx] = 0;
}

template<class T>
void mapMassRealHelper(FlagGrid& flags, Grid<T>& density, const BasicParticleSystem& parts, ParticleDataImpl<T>& source, 
	MACGrid& deltaX, const Grid<Real>& phiObs, Real dt, Real particleMass, bool noDensityClamping = false)
{
	markFluidAndBoundaryCells(parts, flags, deltaX, phiObs);
	Grid<Real> tmp(flags.getParent());
	FlagGrid flagsTmp(flags);
	density.clear();

	knMapLinear<T>( parts, flags, tmp, density, source ); 	//map weights to target grid, source is not needed because all particles have the same mass.
	knComputeDensity<T>(density, flagsTmp, flags, deltaX, dt, particleMass, noDensityClamping);
}

PYTHON() void mapMassToGrid(FlagGrid& flags, Grid<Real>& density, const BasicParticleSystem& parts , ParticleDataImpl<Real>& source, MACGrid& deltaX, 
	                        const Grid<Real>& phiObs, Real dt, Real particleMass, bool noDensityClamping = false) {
	mapMassRealHelper<Real>(flags, density,parts,source, deltaX, phiObs, dt, particleMass, noDensityClamping);
}

//******************************************************************************
// computation of displacements
KERNEL(bnd = 1)
void knRemoveEmptyLambdas(const FlagGrid& flags, Grid<Real>& Lambda)
{
	if (flags.isEmpty(i,j,k)) Lambda(i,j,k) = 0;
}

KERNEL()
void knComputeDeltaX(const FlagGrid& flags, MACGrid& deltaX, const Grid<Real>& Lambda)
{
	if (!flags.isObstacle(i, j, k))		//no displacement into obstacles. The displacement for obstacles with particle inside is already set.
	{
		if (!flags.isObstacle(i - 1, j, k)) deltaX(i, j, k).x = Lambda(i, j, k) - Lambda(i - 1, j, k);
		if (!flags.isObstacle(i, j - 1, k)) deltaX(i, j, k).y = Lambda(i, j, k) - Lambda(i, j - 1, k);
		if (flags.is3D() && !flags.isObstacle(i, j, k - 1)) deltaX(i, j, k).z = Lambda(i, j, k) - Lambda(i, j, k - 1);
	}
}

PYTHON() void computeDeltaX(MACGrid& deltaX, Grid<Real>& Lambda, const FlagGrid& flags)
{
	knRemoveEmptyLambdas(flags, Lambda);
	knComputeDeltaX(flags, deltaX, Lambda);
}

inline void clamp(Vec3& v, const Vec3& min, const Vec3& max)
{
	if (v[0] > max[0]) v[0] = max[0];
	if (v[0] < min[0]) v[0] = min[0];
	if (v[1] > max[1]) v[1] = max[1];
	if (v[1] < min[1]) v[1] = min[1];
	if (v[2] > max[2]) v[2] = max[2];
	if (v[2] < min[2]) v[2] = min[2];
}

//******************************************************************************
// mapping of grid displacements back to the particles
KERNEL(pts)
void knMapLinearMACGridToVec3_Position(BasicParticleSystem& p, const FlagGrid& flags, const MACGrid& deltaX, 
		const ParticleDataImpl<int>* ptype, const int exclude, const Vec3& min, const Vec3& max, Real dt)
{
	if (!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;
	Vec3 dx = deltaX.getInterpolated(p[idx].pos);
	
	p[idx].pos += dx * dt;
	clamp(p[idx].pos, min, max);
}

PYTHON() void mapMACToPartPositions(const FlagGrid& flags, const MACGrid& deltaX,
	BasicParticleSystem& parts, Real dt, const ParticleDataImpl<int>* ptype = NULL, const int exclude = 0, bool mapQuadratic = false) 
{
	Vec3 min, max;
	if (flags.is3D())
	{
		min = Vec3(1.001f, 1.001f, 1.001f);
		max = Vec3(flags.getSizeX() - 1.001f, flags.getSizeY() - 1.001f, flags.getSizeZ() - 1.001f);
	}
	else  // don't clamp z component in 2d 
	{
		min = Vec3(1.001f, 1.001f, -10.001f);
		max = Vec3(flags.getSizeX() - 1.001f, flags.getSizeY() - 1.001f, 10.001f);
	}
	knMapLinearMACGridToVec3_Position(parts, flags, deltaX, ptype, exclude, min, max, dt);
}

//******************************************************************************
// handling of degenerate configurations (too many particles in one cell so that it is likely that particles have the same position).
// We resample the particles by dividing the cell in at least as many subcells as we have particles in the cell so that each particle can be placed in its own subcell. 
// todo: implement this more efficiently. First detect cells where resampling is needed. 
//       Then build gridParticleIndex and map velocities to grid only where needed. Currently this is done for all cells.
PYTHON() void resampeOverfullCells(const MACGrid& vel, Grid<Real>& density, const Grid<int> &index, 
	                                const ParticleIndexSystem& indexSys, BasicParticleSystem& part, ParticleDataImpl<Vec3> &pVel, Real dt)
{
	FOR_IJK(density)
	{
		IndexInt idx = density.index(i, j, k);

		if(density(i,j,k) < -1.0 )	//too many particles in cell -> needs resampling
		{
			// loop for particles in cell
			const int isysIdxS = index.index(i, j, k);
			const int pStart = index(isysIdxS);
			const int pEnd = index.isInBounds(isysIdxS + 1) ? index(isysIdxS + 1) : indexSys.size();
			const int nParts = pEnd - pStart;

			int nCellsPerDir = ceil(cbrt(nParts));
			if (!density.is3D()) nCellsPerDir = ceil(sqrt(nParts));
			int nCells = nCellsPerDir * nCellsPerDir * nCellsPerDir;
			if (!density.is3D()) nCells = nCellsPerDir * nCellsPerDir;

			vector<int> cellIndicesI(nCells);
			vector<int> cellIndicesJ(nCells);
			vector<int> cellIndicesK(nCells);
			for (int l = 0; l < nCellsPerDir; l++)
			{
				if(density.is3D())
					for (int m = 0; m < nCellsPerDir * nCellsPerDir; m++)
					{
						cellIndicesI[l * nCellsPerDir * nCellsPerDir + m] = l;
						cellIndicesJ[l * nCellsPerDir * nCellsPerDir + m] = l;
						cellIndicesK[l * nCellsPerDir * nCellsPerDir + m] = l;
					}
				else
					for (int m = 0; m < nCellsPerDir; m++)
					{
						cellIndicesI[l * nCellsPerDir + m] = l;
						cellIndicesJ[l * nCellsPerDir + m] = l;
					}
			}
				
			std::random_shuffle(cellIndicesI.begin(), cellIndicesI.end());
			std::random_shuffle(cellIndicesJ.begin(), cellIndicesJ.end());
			if(density.is3D()) std::random_shuffle(cellIndicesK.begin(), cellIndicesK.end());
			int l = 0;	//position in the cellIndices vector
			Vec3 cellPos(i, j, k);
				
			for (int p = pStart; p < pEnd; ++p)
			{
				const int psrc = indexSys[p].sourceIndex;
				Vec3 &xj = part[psrc].pos;
					
				if(density.is3D())
					xj = cellPos + Vec3( (cellIndicesI[l] + 0.5f) / static_cast<Real>(nCellsPerDir),
										 (cellIndicesJ[l] + 0.5f) / static_cast<Real>(nCellsPerDir), 
										 (cellIndicesK[l] + 0.5f) / static_cast<Real>(nCellsPerDir));
				else
					xj = cellPos + Vec3((cellIndicesI[l] + 0.5f) / static_cast<Real>(nCellsPerDir),
										(cellIndicesJ[l] + 0.5f) / static_cast<Real>(nCellsPerDir),
											0);
				pVel[psrc] = vel.getInterpolated(xj);

				l++;
			}
			density[idx] = -1.0;

		}
		else if (density(i, j, k) < -0.5)
			density[idx] = -0.5;

		if (density(i, j, k) > 0.5)	//too few particles in cell -> no resampling needed, just clamp
		{
			if (density[idx] > 0.5)
				density[idx] = 0.5;			
		}
		density[idx] = density[idx] / dt;
	}
}

//******************************************************************************
// helper for initialization
PYTHON() void copyFlagsToFlags(FlagGrid &source, FlagGrid &target)
{
	FOR_IJK(target) {
		target(i, j, k) = source(i, j, k);
	}
}

} // namespace

