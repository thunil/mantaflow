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

//! Special function for levelsets
PYTHON class LevelsetGrid : public Grid<Real> {
public:
    PYTHON LevelsetGrid(FluidSolver* parent, bool show = true);
    
    PYTHON void reinitMarching(FlagGrid& flags, Real maxTime=4.0, MACGrid* velTransport=NULL, bool ignoreWalls=false, bool correctOuterLayer=true);
    
    static Real invalidTimeValue();
};

} //namespace
#endif