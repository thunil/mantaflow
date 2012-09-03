/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Tools to setup fields and inflows
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "shapes.h"
#include "commonkernels.h"
#include "noisefield.h"

using namespace std;

namespace Manta {
    
//! Apply noise to grid
KERNEL 
void KnApplyNoise(FlagGrid& flags, Grid<Real>& dens, WaveletNoiseField& noise, Grid<Real>& sdf, Real scale, Real sigma) 
{
    if (!flags.isFluid(i,j,k) || sdf(i,j,k) > sigma) return;
    Real factor = clamp(1.0f-0.5f/sigma * (sdf(i,j,k)+sigma), 0.0f, 1.0f);
    
    Real target = noise.evaluate(Vec3(i,j,k)) * scale * factor;
    if (dens(i,j,k) < target)
        dens(i,j,k) = target;
}

//! Init noise-moduled density inside shape
PYTHON void densityInflow(FlagGrid& flags, Grid<Real>& density, WaveletNoiseField& noise, Shape* shape, Real scale=1.0, Real sigma=0)
{
    Grid<Real> sdf = shape->computeLevelset();    
    KnApplyNoise(flags, density, noise, sdf, scale, sigma);
}

    
} // namespace