/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Common grid kernels
 *
 ******************************************************************************/

#ifndef _COMMONKERNELS_H
#define _COMMONKERNELS_H

#include "general.h"
#include "kernel.h"
#include "grid.h"

namespace Manta {
   
//! Kernel: Invert real values, if positive and fluid
KERNEL(idx) 
void InvertCheckFluid (FlagGrid& flags, Grid<Real>& grid) 
{
    if (flags.isFluid(idx) && grid[idx] > 0)
        grid[idx] = 1.0 / grid[idx];
}

//! Kernel: Squared sum over grid
KERNEL(idx, reduce) 
struct GridSumSqr (Grid<Real>& grid)
{
    double sum;
    
    void operator()(int idx) {
        sum += square((double)grid[idx]);
    }
    void setup() { sum = 0; }
    void join(const GridSumSqr& a) { sum += a.sum; }
};

//! Kernel: rotation operator \nabla x v for centered vector fields
KERNEL(bnd=1) 
void CurlOp (const Grid<Vec3>& grid, Grid<Vec3>& dst) {
    dst(i,j,k) = Vec3(0.5*((grid(i,j+1,k).z - grid(i,j-1,k).z) - (grid(i,j,k+1).y - grid(i,j,k-1).y)),
                      0.5*((grid(i,j,k+1).x - grid(i,j,k-1).x) - (grid(i+1,j,k).z - grid(i-1,j,k).z)),
                      0.5*((grid(i+1,j,k).y - grid(i-1,j,k).y) - (grid(i,j+1,k).x - grid(i,j-1,k).x)));
};

//! Kernel: divergence operator (from MAC grid)
KERNEL(bnd=1) 
void DivergenceOpMAC(Grid<Real>& div, const MACGrid& grid) {
    Vec3 del = Vec3(grid(i+1,j,k).x, grid(i,j+1,k).y, grid(i,j,k+1).z) - grid(i,j,k);
    div(i,j,k) = del.x + del.y + del.z;
}

//! Kernel: gradient operator (create MAC grid)
KERNEL(bnd=1)void GradientOpMAC(MACGrid& gradient, const Grid<Real>& grid) {
    gradient(i,j,k) = (Vec3(grid(i,j,k)) - Vec3(grid(i-1,j,k), grid(i,j-1,k), grid(i,j,k-1)));
}

//! Kernel: gradient operator 
KERNEL(bnd=1) void GradientOp(Grid<Vec3>& gradient, const Grid<Real>& grid) {
    gradient(i,j,k) = 0.5 * Vec3( grid(i+1,j,k)-grid(i-1,j,k), grid(i,j+1,k)-grid(i,j-1,k), grid(i,j,k+1)-grid(i,j,k-1));
}

//! Kernel: get component at MAC positions
KERNEL(bnd=1) void GetShiftedComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) {
    Vec3i ishift(i,j,k);
    ishift[dim]--;
    comp(i,j,k) = 0.5*(grid(i,j,k)[dim] + grid(ishift)[dim]);
};

//! Kernel: get component (not shifted)
KERNEL(idx) void GetComponent(const Grid<Vec3>& grid, Grid<Real>& comp, int dim) {
    comp[idx] = grid[idx][dim];
};

//! Kernel: get norm of centered grid
KERNEL(idx) void GridNorm(Grid<Real>& n, const Grid<Vec3>& grid) {
    n[idx] = norm(grid[idx]);
};

//! Kernel: set component (not shifted)
KERNEL(idx) void SetComponent(Grid<Vec3>& grid, const Grid<Real>& comp, int dim) {
    grid[idx][dim] = comp[idx];
};

//! Kernel: compute centered velocity field from MAC
KERNEL(bnd=1) void GetCentered(Grid<Vec3>& center, const MACGrid& vel) {
    center(i,j,k) = 0.5*(vel(i,j,k)+Vec3(vel(i+1,j,k).x, vel(i,j+1,k).y, vel(i,j,k+1).z));
};

//! Kernel: compute MAC from centered velocity field
KERNEL(bnd=1) void GetMAC(MACGrid& vel, const Grid<Vec3>& center) {
    vel(i,j,k) = 0.5*(center(i,j,k)+Vec3(center(i-1,j,k).x, center(i,j-1,k).y, center(i,j,k-1).z));
};

} // namespace
#endif