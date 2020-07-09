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

import os
import numpy
import argparse
import datetime
import time
from random import randint, seed
from scenes.scene import Scene
import scenes.volumes as volumes

def random_vec3(vmin=[-1,-1,-1], vmax=[1,1,1]):
    return vec3(
        numpy.random.uniform(low=vmin[0], high=vmax[0]),
        numpy.random.uniform(low=vmin[1], high=vmax[1]),
        numpy.random.uniform(low=vmin[2], high=vmax[2])
    )

def random_vec3s(vmin=-1, vmax=1): # scalar params
    return vec3(
        numpy.random.uniform(low=vmin, high=vmax),
        numpy.random.uniform(low=vmin, high=vmax),
        numpy.random.uniform(low=vmin, high=vmax)
    )


class SimpleScene(Scene):
    #file_num = 0
    open_bound = True
    sources = []
    source_strengths = []

#     #----------------------------------------------------------------------------------
#     def set_velocity(self, volume, velocity):
#         super(SmokeScene,self).set_velocity(volume, velocity)
# 
#     #----------------------------------------------------------------------------------
#     def add_sink(self, volume):
#         print("WARNING - sinks not yet supported for smoke scene")
# 
#     #----------------------------------------------------------------------------------
#     # sources used as smoke inflow in the following
#     def add_source(self, volume):
#         shape = volume.shape(self.solver)
#         self.sources.append(shape)
# 
#     #----------------------------------------------------------------------------------
#     # this is kind of a hack, since for smoke sources are much more desirable than just adding a fluid once
#     def add_fluid(self, volume):
#         self.add_source(volume)
# 
    def _init(self): 
        # solver
        self.max_iter_fac = 2
        self.accuracy =5e-4

        # grids
        self.density = self.solver.create(RealGrid, name="Density")
        self.pressure_static = self.solver.create(RealGrid, name="StatPressure_dummy_") # not used for smoke

        noise = self.solver.create(NoiseField, loadFromFile=True)
        noise.posScale = vec3(40) * numpy.random.uniform(low=0.25, high=1.)
        noise.posOffset = random_vec3s( vmin=0. ) * 100.
        noise.clamp = True
        noise.clampNeg = 0
        noise.clampPos = 1.
        noise.valOffset = 0.15
        noise.timeAnim  = 0.4 * numpy.random.uniform(low=0.2, high=1.)
        self.noise = noise

        self.source_strengths = []
        for i in range(100): # some more for safety
            self.source_strengths.append( numpy.random.uniform(low=0.5, high=1.) )

    #----------------------------------------------------------------------------------
    def _create_scene(self): 
        # from reset
        self.sources = [] 
        self.density.setConst(0) 
        self.vel.setConst(vec3(0)) 
        is3d = (self.dimension > 2)

        # TODO , needed?
        # dont reset! multiple sims written into single file with increasing index...
        #self.file_num = 0 

        source_count = randint(1, 2)  #  from scene_setup.scene_selection
        for i in range(source_count):
            print("_create_scene rand box ")
            box = volumes.random_box(center_min=[0.2, 0.1, 0.2], center_max=[0.8, 0.6, 0.8], size_min=[0.005, 0.005, 0.005], size_max=[0.2, 0.2, 0.2], is3d=is3d)

        super(SimpleScene, self)._create_scene()
        self.flags.initDomain(boundaryWidth=self.boundary) 
        self.flags.fillGrid()
        if self.open_bound:
            setOpenBound(self.flags, self.boundary, 'yY', 16 | 4) # FlagOutflow|FlagEmpty) 

        print("_create_scene done")

    #==================================================================================
    # SIMULATION
    #----------------------------------------------------------------------------------

    def _compute_simulation_step(self):
        # Add source
        # randomize noise offset , note - sources are turned off earlier, the more there are
        for i in range(len(self.sources)):
            if self.solver.frame<i*(100./len(self.sources)):
                src,sstr = self.sources[i], self.source_strengths[i]
                densityInflow(flags=self.flags, density=self.density, noise=self.noise, shape=src, scale=2.0*sstr, sigma=0.5)

        print("Simulation step simple")
        advectSemiLagrange(flags=self.flags, vel=self.vel, grid=self.density, order=2, clampMode=2)
        advectSemiLagrange(flags=self.flags, vel=self.vel, grid=self.vel    , order=2, clampMode=2)

        vorticityConfinement( vel=self.vel, flags=self.flags, strength=0.1 )
        addBuoyancy(density=self.density, vel=self.vel, gravity=0.2*self.gravity, flags=self.flags)

        setWallBcs(flags=self.flags, vel=self.vel)

        solvePressure(flags=self.flags, vel=self.vel, pressure=self.pressure, cgMaxIterFac=self.max_iter_fac, cgAccuracy=self.accuracy)
        if self.boundary>0:
            self.pressure.setBoundNeumann(self.boundary-1)

        # done in solveP, correctVelocities(vel=self.vel, pressure=self.pressure, flags=self.flags)
        self.vel.setBoundNeumann(self.boundary)

