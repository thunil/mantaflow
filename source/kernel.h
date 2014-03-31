/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Function and macros for defining compution kernels over grids
 *
 ******************************************************************************/

#ifndef _KERNEL_H
#define _KERNEL_H

#ifdef TBB
#   include <tbb/blocked_range3d.h>
#   include <tbb/blocked_range.h>
#   include <tbb/parallel_for.h>
#   include <tbb/parallel_reduce.h>
#endif

#ifdef OPENMP
#	include <omp.h>
#endif

namespace Manta {
// fwd decl
class GridBase;
class ParticleBase;
    
    
// simple iteration
#define FOR_IJK_BND(grid, bnd) \
    for(int k=((grid).is3D() ? bnd : 0),__kmax=((grid).is3D() ? ((grid).getSizeZ()-bnd) : 1); k<__kmax; k++) \
        for(int j=bnd; j<(grid).getSizeY()-bnd; j++) \
            for(int i=bnd; i<(grid).getSizeX()-bnd; i++)

#define FOR_IJK_REVERSE(grid) \
    for(int k=(grid).getSizeZ()-1; k>=0; k--) \
        for(int j=(grid).getSizeY()-1; j>=0; j--) \
            for(int i=(grid).getSizeX()-1; i>=0; i--)

#define FOR_IDX(grid) \
    for(int idx=0, total=(grid).getSizeX()*(grid).getSizeY()*(grid).getSizeZ(); idx<total; idx++)

#define FOR_IJK(grid) FOR_IJK_BND(grid, 0)
               
    
struct KernelBase {
    int maxX, maxY, maxZ, maxCells, minZ;
    int X, Y, Z;
	//! store thread info for this kernel 
	int threadId, threadNum;
    
    KernelBase(const GridBase* base, int bnd);
    KernelBase(int _maxX, int _maxY, int _maxZ, int _maxC, int _minZ, int _X, int _Y, int _Z);
    
    // specify in your derived classes:
    
    // kernel operators    
    // ijk mode: void operator() (size_t idx)
    // idx mode: void operator() (size_t i, size_t j, size_t k)
    
    // reduce mode: 
    // void join(classname& other)
    // void setup()    
};

struct ParticleKernelBase {
    int size;
	//! store thread info for this kernel 
	int threadId, threadNum;
    
    ParticleKernelBase(int sz);
    
    // specify in your derived classes:
    
    // kernel operators    
    // void operator() (size_t idx)
    
    // reduce mode: 
    // void join(classname& other)
    // void setup()    
};

} // namespace

// Define plugin documentation group
// all kernels will automatically be added to this group
/*! @defgroup Kernels Computation Kernels
 */

#endif
