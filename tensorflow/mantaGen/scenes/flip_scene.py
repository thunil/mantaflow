#******************************************************************************
#
# MantaGen
# Copyright 2018 Steffen Wiewel, Moritz Becher, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
#******************************************************************************

from manta import *

from random import randint
from scenes.scene import Scene
from scenes.volumes import *
from scenes.functions import *
from util.logger import *


def instantiate_scene(**kwargs): # instantiate independent of name , TODO replace?
    info(kwargs)
    return FLIPScene(**kwargs) 


class FLIPScene(Scene):
    #----------------------------------------------------------------------------------
    def __init__(self, **kwargs):
        super(FLIPScene,self).__init__(**kwargs)

        # const
        self.max_iter_fac       = 2
        self.accuracy           = 5e-4
        self.narrow_band        = 3 

        # arguments
        self.max_drop_count   = int(kwargs.get("max_drop_count", 3))
        self.max_obstacle_count   = int(kwargs.get("max_obstacle_count", 0))

        # grids and pp
        self.vel_org = self.solver.create(MACGrid, name="VelOrg")
        self.vel_parts = self.solver.create(MACGrid, name="VelParts")
        self.tmp_vec3 = self.solver.create(MACGrid, name="temp_vec3")
        self.phi_parts = self.solver.create(LevelsetGrid, name="PhiParticles")

        self.pp     = self.solver.create(BasicParticleSystem) 
        self.pVel   = self.pp.create(PdataVec3) 
        self.pindex = self.solver.create(ParticleIndexSystem) 
        self.gpi    = self.solver.create(IntGrid)

        info("FLIPScene initialized")

    #----------------------------------------------------------------------------------
    def set_velocity(self, volume, velocity):
        if self.dimension == 2:
            velocity.z = 0.0
        volume.applyToGrid(solver=self.solver, grid=self.vel, value=velocity)
        mapGridToPartsVec3(source=self.vel, parts=self.pp, target=self.pVel)

    #----------------------------------------------------------------------------------
    def add_fluid(self, volume):
        vol_ls = volume.computeLevelset(self.solver)
        self.phi_fluid.join(vol_ls)

    #----------------------------------------------------------------------------------
    def add_obstacle(self, volume):
        vol_ls = volume.computeLevelset(self.solver)
        self.phi_obs.join(vol_ls)

    #----------------------------------------------------------------------------------
    # def add_source(self, volume):
    #     print("WARNING - sources not yet supported for FLIP scene")
    #     vol_ls = volume.computeLevelset(self.solver)
    #     self.phi_source.join(vol_ls)

    #----------------------------------------------------------------------------------
    def _create_scene(self):
        super(FLIPScene, self)._create_scene()

        # reset fields
        self.vel_org.setConst(vec3(0,0,0))
        self.vel_parts.setConst(vec3(0,0,0))
        self.tmp_vec3.setConst(vec3(0,0,0))
        self.phi_parts.setConst(0) 

        self.pp.clear()
        self.pVel.setConst(vec3(0,0,0))
        self.pindex.clear()
        self.gpi.setConst(0)
        is3d = (self.dimension > 2)

        ##### spawn scene
        # low basin is the large volume of liquid at the bottom of the scene. in general adds most of the liquid in the scene
        low_basin = random_box(center_min=[0.5, 0.0, 0.5], center_max=[0.5, 0.0, 0.5], size_min=[1, 0.3, 1], size_max=[1, 0.6, 1], is3d=is3d)
        self.add_fluid(low_basin)
        low_basin_velo = random_velocity(vel_min=[-1.5, 0, -1.5], vel_max=[1.5, 0, 1.5])
        self.set_velocity(low_basin, low_basin_velo)

        # the high basin is a tall but narrow block of liquid, adding motion to the scene
        high_basin = random_box(center_min=[0.0, 0.3, 0.0], center_max=[1.0, .3, 1.0], size_min=[0.1, 0.2, 0.1], size_max=[0.4, 0.5, 0.4], is3d=is3d)
        self.add_fluid(high_basin)
        high_basin_velo =  random_velocity(vel_min=[-2.5, -2.5, -2.5], vel_max=[2.5, 0, 2.5])
        self.set_velocity(high_basin, high_basin_velo)

        # drops: spawn at most 3 drops or meshes
        drop_count = randint(0, self.max_drop_count)
        for i in range(drop_count):
            box = random_box(center_min=[0, 0.7, 0], center_max=[1, 0.9, 1], size_min=[0.005, 0.005, 0.005], size_max=[0.25, 0.25, 0.25], is3d=is3d)
            drop_velo = random_velocity(vel_min=[-2.5, -2.5, -2.5], vel_max=[2.5, 0, 2.5])
            self.add_fluid(box)
            self.set_velocity(box, drop_velo)

        # optional: add solid box to scene
        if self.max_obstacle_count > 0:
            obstacle_count = randint(0, self.max_obstacle_count)
            for i in range(obstacle_count):
                obstacle = random_box(center_min=[0.1,0.1,0.1], center_max=[1.0,0.5,1.0], size_min=[0.2,0.2,0.2], size_max=[0.3,1.0,0.3], is3d=is3d)
                self.add_obstacle(obstacle)

        # extrapolate velocities from 1-cell inside towards empty region for particles
        self.phi_fluid.addConst( 1.)
        self.flags.updateFromLevelset(self.phi_fluid)
        self.phi_fluid.addConst(-1.)
        extrapolateMACSimple(flags=self.flags, vel=self.vel , distance=3, intoObs=True)
        sampleLevelsetWithParticles(phi=self.phi_fluid, flags=self.flags, parts=self.pp, discretization=2, randomness=0.05)
        mapGridToPartsVec3(source=self.vel, parts=self.pp, target=self.pVel)

    #==================================================================================
    # SIMULATION
    #----------------------------------------------------------------------------------
    def _compute_simulation_step(self):
        self._advect(ls_order=1)
        self._enforce_boundaries(distance=2)
        self._solve_pressure(max_iter_fac=self.max_iter_fac, accuracy=self.accuracy)

    #----------------------------------------------------------------------------------
    def _advect(self, extrapol_dist=3, ls_order=2):
        # extrapolate the grids into empty cells
        #self.phi_fluid.reinitMarching(flags=self.flags, velTransport=self.vel, maxTime=4.0)

        self.pp.advectInGrid(flags=self.flags, vel=self.vel, integrationMode=2, deleteInObstacle=False, stopInObstacle=False )
        pushOutofObs( parts=self.pp, flags=self.flags, phiObs=self.phi_obs )

        # make sure nothings sticks to the top... (helper in test.cpp)
        #deleteTopParts( parts=self.pp, phi=self.phi_fluid, maxHeight=self.resolution-1-(self.boundary+2) ) # needs 2 more to make sure it's out of the setBoundNeumann range 

        # advect the levelset
        advectSemiLagrange(flags=self.flags, vel=self.vel, grid=self.phi_fluid, order=ls_order) 

        # velocity self-advection
        advectSemiLagrange(flags=self.flags, vel=self.vel, grid=self.vel, order=2)

        # particle SDF
        gridParticleIndex( parts=self.pp , flags=self.flags, indexSys=self.pindex, index=self.gpi )
        unionParticleLevelset( self.pp, self.pindex, self.flags, self.gpi, self.phi_parts )

        # TODO: source & sink , not yet working...
        #self.phi_fluid.subtract(self.phi_sink) 
        #self.phi_fluid.join(self.phi_source)

        # combine level set of particles with grid level set
        self.phi_fluid.addConst(1.) # shrink slightly
        self.phi_fluid.join( self.phi_parts )
        extrapolateLsSimple(phi=self.phi_fluid, distance=self.narrow_band+2, inside=True ) 
        extrapolateLsSimple(phi=self.phi_fluid, distance=3 )

        # enforce boundaries
        self.phi_fluid.setBoundNeumann(self.boundary-1 if self.boundary>0 else 0) # 1 cell less...
        self.flags.updateFromLevelset(self.phi_fluid)

        # combine particles velocities with advected grid velocities
        mapPartsToMAC(vel=self.vel_parts, flags=self.flags, velOld=self.vel_org, parts=self.pp, partVel=self.pVel, weight=self.tmp_vec3)
        extrapolateMACFromWeight( vel=self.vel_parts , distance=2, weight=self.tmp_vec3 )
        combineGridVel(vel=self.vel_parts, weight=self.tmp_vec3 , combineVel=self.vel, phi=self.phi_fluid, narrowBand=(self.narrow_band-1), thresh=0)
        self.vel_org.copyFrom(self.vel)

        addGravity(flags=self.flags, vel=self.vel, gravity=self.gravity)

    #----------------------------------------------------------------------------------
    def _enforce_boundaries(self, distance):
        extrapolateMACSimple( flags=self.flags, vel=self.vel , distance=distance, intoObs=True )
        setWallBcs(flags=self.flags, vel=self.vel, fractions=self.fractions, phiObs=self.phi_obs)

    #----------------------------------------------------------------------------------
    def _solve_pressure(self, max_iter_fac=2, accuracy=5e-4):
        solvePressure(vel=self.vel, pressure=self.pressure, flags=self.flags, fractions=self.fractions, cgAccuracy=accuracy, cgMaxIterFac=max_iter_fac, phi=self.phi_fluid)

        # remove pressure discontinuity at boundary - note: only outer boundary, does not influence the required pressure gradients in any way
        if self.boundary>0:
            self.pressure.setBoundNeumann(self.boundary-1)

        self._enforce_boundaries(4)

        minParticles  = pow(2,self.dimension)
        self.pVel.setSource( self.vel, isMAC=True )
        adjustNumber( parts=self.pp, vel=self.vel, flags=self.flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=self.phi_fluid, exclude=self.phi_obs, narrowBand=self.narrow_band ) 
        flipVelocityUpdate(vel=self.vel, velOld=self.vel_org, flags=self.flags, parts=self.pp, partVel=self.pVel, flipRatio=0.97 )

