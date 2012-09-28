/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Set boundary conditions, gravity
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "commonkernels.h"

using namespace std;

namespace Manta { 

//! add Forces between fl/fl and fl/em cells
KERNEL(bnd=1) void KnAddForceField(FlagGrid& flags, MACGrid& vel, Grid<Vec3>& force) {
    bool curFluid = flags.isFluid(i,j,k);
    bool curEmpty = flags.isEmpty(i,j,k);
    if (!curFluid && !curEmpty) return;
    
    if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
        vel(i,j,k).x += 0.5*(force(i-1,j,k).x + force(i,j,k).x);
    if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
        vel(i,j,k).y += 0.5*(force(i,j-1,k).y + force(i,j,k).y);
    if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
        vel(i,j,k).z += 0.5*(force(i,j,k-1).z + force(i,j,k).z);
}

//! add Forces between fl/fl and fl/em cells
KERNEL(bnd=1) void KnAddForce(FlagGrid& flags, MACGrid& vel, Vec3 force) {
    bool curFluid = flags.isFluid(i,j,k);
    bool curEmpty = flags.isEmpty(i,j,k);
    if (!curFluid && !curEmpty) return;
    
    if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
        vel(i,j,k).x += force.x;
    if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
        vel(i,j,k).y += force.y;
    if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
        vel(i,j,k).z += force.z;
}

//! add gravity forces to all fluid cells
PYTHON void addGravity(FlagGrid& flags, MACGrid& vel, Vec3 gravity) {    
    Vec3 f = gravity * parent->getDt() / flags.getDx();
    KnAddForce(flags, vel, f);
}

//! add Buoyancy force based on smoke density
KERNEL(bnd=1) void KnAddBuoyancy(FlagGrid& flags, Grid<Real>& density, MACGrid& vel, Vec3 strength) {    
    if (!flags.isFluid(i,j,k)) return;
    if (flags.isFluid(i-1,j,k))
        vel(i,j,k).x += (0.5 * strength.x) * (density(i,j,k)+density(i-1,j,k));
    if (flags.isFluid(i,j-1,k))
        vel(i,j,k).y += (0.5 * strength.y) * (density(i,j,k)+density(i,j-1,k));
    if (vel.is3D() && flags.isFluid(i,j,k-1))
        vel(i,j,k).z += (0.5 * strength.z) * (density(i,j,k)+density(i,j,k-1));    
}

//! add Buoyancy force based on smoke density
PYTHON void addBuoyancy(FlagGrid& flags, Grid<Real>& density, MACGrid& vel, Vec3 gravity) {
    Vec3 f = - gravity * parent->getDt() / parent->getDx();
    KnAddBuoyancy(flags,density, vel, f);
}

//! DDF style no-slip (what does it actually do?)
KERNEL(bnd=1) void KnSetNoSlipBcs(FlagGrid& flags, MACGrid& vel) {
    if (flags.isObstacle(i,j,k) && !flags.isInflow(i,j,k))  { 
        vel(i,j,k) = 0; 
        return;        
    }
    if (flags.isEmpty(i,j,k)) 
        return;
    if (flags.isFluid(i,j,k)) {
        if (flags.isObstacle(i-1,j,k) && !flags.isInflow(i,j,k)) vel(i,j,k).x = 0;
        if (flags.isObstacle(i,j-1,k) && !flags.isInflow(i,j,k)) vel(i,j,k).y = 0;
        if (flags.is3D() && flags.isObstacle(i,j,k-1) && !flags.isInflow(i,j,k)) vel(i,j,k).z = 0;
    }
}
        
// Nils: unless I'm misunderstanding something
// noslip in DDF is actually no-stick, DDF freeselip does something weird ?!

//! set no-stick wall boundary condition between ob/fl and ob/ob cells
KERNEL void KnSetWallBcs(FlagGrid& flags, MACGrid& vel) {
    bool curFluid = flags.isFluid(i,j,k);
    bool curObstacle = flags.isObstacle(i,j,k);
    if (!curFluid && !curObstacle) return;
    
    // we use i>0 instead of bnd=1 to check outer wall
    if (i>0 && (flags.isObstacle(i-1,j,k) || (curObstacle && flags.isFluid(i-1,j,k))))
        vel(i,j,k).x = 0;
    if (j>0 && (flags.isObstacle(i,j-1,k) || (curObstacle && flags.isFluid(i,j-1,k))))
        vel(i,j,k).y = 0;
    if (vel.is2D() || (k>0 && (flags.isObstacle(i,j,k-1) || (curObstacle && flags.isFluid(i,j,k-1)))))
        vel(i,j,k).z = 0;
		
	if (curFluid) {
		if ((i>0 && flags.isStick(i-1,j,k)) || (i<flags.getSizeX()-1 && flags.isStick(i+1,j,k)))
			vel(i,j,k).y = vel(i,j,k).z = 0;
		if ((j>0 && flags.isStick(i,j-1,k)) || (j<flags.getSizeY()-1 && flags.isStick(i,j+1,k)))
			vel(i,j,k).x = vel(i,j,k).z = 0;
		if (vel.is3D() && ((k>0 && flags.isStick(i,j,k-1)) || (k<flags.getSizeZ()-1 && flags.isStick(i,j,k+1))))
			vel(i,j,k).x = vel(i,j,k).y = 0;
	}
}

//! set no-stick boundary condition on walls
PYTHON void setWallBcs(FlagGrid& flags, MACGrid& vel) {
    //KnSetWallBcs(flags, vel);
    KnSetNoSlipBcs(flags, vel);
} 

//! set boundary conditions at empty cells
KERNEL(bnd=1) void KnSetLiquidBcs(FlagGrid& flags, MACGrid& vel) {
    if (!flags.isFluid(i,j,k)) return;
    
    // init empty cells from fluid
    if (flags.isEmpty(i+1,j,k)) 
        vel(i+1,j,k).x = vel(i,j,k).x;
    if (flags.isEmpty(i,j+1,k)) 
        vel(i,j+1,k).y = vel(i,j,k).y;
    if (flags.isEmpty(i,j,k+1)) 
        vel(i,j,k+1).z = vel(i,j,k).z;
    
    // "left" sides of fluid - fluid cells neighboring
    // empty cells along neg. dir, get velocities from fluid 
    if (flags.isEmpty(i-1,j,k) && flags.isFluid(i+1,j,k)) 
        vel(i,j,k).x = vel(i+1,j,k).x;
    if (flags.isEmpty(i,j-1,k) && flags.isFluid(i,j+1,k)) 
        vel(i,j,k).y = vel(i,j+1,k).y;
    if (flags.isEmpty(i,j,k-1) && flags.isFluid(i,j,k+1)) 
        vel(i,j,k).z = vel(i,j,k+1).z;
}

//! set boundary conditions at empty cells
PYTHON void setLiquidBcs(FlagGrid& flags, MACGrid& vel) {
    FOR_IDX(flags) {
        if (flags.isEmpty(idx)) vel[idx]=_0;
    }
    KnSetLiquidBcs(flags, vel);
} 

//! Kernel: gradient norm operator
KERNEL(bnd=1) void KnConfForce(Grid<Vec3>& force, const Grid<Real>& grid, const Grid<Vec3>& curl, Real str) {
    Vec3 grad = 0.5 * Vec3( grid(i+1,j,k)-grid(i-1,j,k), grid(i,j+1,k)-grid(i,j-1,k), grid(i,j,k+1)-grid(i,j,k-1));
    normalize(grad);
    force(i,j,k) = str*cross(grad, curl(i,j,k));
}

PYTHON void vorticityConfinement(MACGrid& vel, FlagGrid& flags, Real strength) {
    assertMsg(vel.is3D(), "Only 3D grids supported so far");
    Grid<Vec3> velCenter(parent), curl(parent), force(parent);
    Grid<Real> norm(parent);
    
    GetCentered(velCenter, vel);
    CurlOp(velCenter, curl);
    GridNorm(norm, curl);
    KnConfForce(force, norm, curl, strength);
    KnAddForceField(flags, vel, force);
}


} // namespace