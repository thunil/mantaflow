/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Levelset
 *
 ******************************************************************************/

#ifndef _LEVELSET_H_
#define _LEVELSET_H_

#include "grid.h"

namespace Manta {
class Mesh;

//! Special function for levelsets
PYTHON template<int DIM> 
class LevelsetGrid : public Grid<DIM, Real> {
public:
    PYTHON LevelsetGrid(FluidSolver* parent, bool show = true);
    
    //! reconstruct the levelset using fast marching
    PYTHON void reinitMarching(FlagGrid<DIM>& flags, Real maxTime=4.0, MACGrid<DIM>* velTransport=NULL, bool ignoreWalls=false, bool correctOuterLayer=true);
    //! create a triangle mesh from the levelset isosurface
    PYTHON void createMesh(Mesh& mesh);
    
    static Real invalidTimeValue();
    using Grid<DIM, Real>::get;
};

// named instances
PYTHON alias LevelsetGrid<2> LevelsetGrid2;
PYTHON alias LevelsetGrid<3> LevelsetGrid3;

} //namespace
#endif