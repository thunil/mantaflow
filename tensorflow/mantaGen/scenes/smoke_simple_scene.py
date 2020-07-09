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
    return SmokeSimpleScene(**kwargs) 


class SmokeSimpleScene(Scene):
    #----------------------------------------------------------------------------------
    def __init__(self, **kwargs):
        super(SmokeSimpleScene,self).__init__(**kwargs)
        # optionally, init more grids etc.

        self.max_iter_fac       = 2
        self.accuracy           = 5e-4
        self.max_source_count   = int(kwargs.get("max_source_count", 2))
        self.velocity_scale     = float(kwargs.get("velocity_scale", self.resolution.y * 0.05))
        self.use_inflow_sources = kwargs.get("use_inflow_sources", "True") == "True"
        self.open_bound         = kwargs.get("use_open_bound", "True") == "True"
        self.sources            = []

        info("SmokeSimpleScene initialized")

    #----------------------------------------------------------------------------------
    def set_velocity(self, volume, velocity):
        if self.dimension == 2:
            velocity.z = 0.0
        volume.applyToGrid(solver=self.solver, grid=self.vel, value=velocity)

    #----------------------------------------------------------------------------------
    def _create_scene(self):
        super(SmokeSimpleScene, self)._create_scene()

        self.sources = []
        self.vel.setConst(vec3(0))

        self.flags.initDomain(boundaryWidth=self.boundary) 
        self.flags.fillGrid()

        if self.open_bound:
            setOpenBound(self.flags, self.boundary, 'yY', CellType_TypeOutflow|CellType_TypeEmpty)

        is3d = (self.dimension > 2)
        source_count = randint(1, self.max_source_count)
        for i in range(source_count):
            volume = random_box(center_min=[0.2, 0.1, 0.2], center_max=[0.8, 0.6, 0.8], size_min=[0.005, 0.005, 0.005], size_max=[0.2, 0.2, 0.2], is3d=is3d)
            velo = random_vec3( vmin=[-self.velocity_scale, -self.velocity_scale, -self.velocity_scale],
                                vmax=[self.velocity_scale, self.velocity_scale, self.velocity_scale])
            self.sources.append( (volume, velo) )
            # set velocity on startup if needed
            if not self.use_inflow_sources:
                self.set_velocity( volume, velo )

        if self.show_gui:
            # central view is more interesting for smoke
            self._gui.setPlane( self.resolution.z // 2 )

        info("SmokeSimpleScene created with {} sources".format(len(self.sources)))

    #==================================================================================
    # SIMULATION
    #----------------------------------------------------------------------------------
    def _compute_simulation_step(self):
        advectSemiLagrange(flags=self.flags, vel=self.vel, grid=self.vel, order=2, clampMode=2)
        # apply velocity source
        if self.use_inflow_sources:
            for vol, vel in self.sources:
                self.set_velocity( vol, vel )

        vorticityConfinement( vel=self.vel, flags=self.flags, strength=0.1 )

        setWallBcs(flags=self.flags, vel=self.vel)
        solvePressure(flags=self.flags, vel=self.vel, pressure=self.pressure, cgMaxIterFac=self.max_iter_fac, cgAccuracy=self.accuracy)