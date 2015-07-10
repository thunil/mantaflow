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
PYTHON() class LevelsetGrid : public Grid<Real> {
public:
	PYTHON() LevelsetGrid(FluidSolver* parent, bool show = true);
	
	//! reconstruct the levelset using fast marching
	PYTHON() void reinitMarching(FlagGrid& flags, Real maxTime=4.0, 
			MACGrid* velTransport=NULL, bool ignoreWalls=false, bool correctOuterLayer=true, 
			int obstacleType = FlagGrid::TypeObstacle, Grid<Real>* scalarTransport = NULL );

    PYTHON() void reinitExact(FlagGrid& flags);

	//! create a triangle mesh from the levelset isosurface
	PYTHON() void createMesh(Mesh& mesh);
	
	//! union with another levelset
	PYTHON() void join(const LevelsetGrid& o);

    //! intersection with another levelset
	PYTHON() void intersect(const LevelsetGrid& o);

    //! difference with another levelset
	PYTHON() void difference(const LevelsetGrid& o);
	
	//! initialize levelset from flags (+/- 0.5 heaviside)
	PYTHON() void initFromFlags(FlagGrid& flags, bool ignoreWalls=false);

    //! set grid to  value
    PYTHON() void applyToGrid(GridBase* grid, FlagGrid* respectFlags=0, Real isoval=0.);	

	static Real invalidTimeValue();
};

} //namespace
#endif
