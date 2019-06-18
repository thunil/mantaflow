#******************************************************************************
#
# MantaGen
# Copyright 2018 Steffen Wiewel, Moritz Becher, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# Basic simulation scene class
# Contains only simulation code, I/O specific code is external, typically create_dataset.py
# See "on_simulation_step" call below
#
#******************************************************************************

from __future__ import print_function
import time
from manta import *
#import scenes.volumes

class Scene(object):
    """ Superclass for all simulation scenes """
    #----------------------------------------------------------------------------------
    def __init__(self, resolution=vec3(64,64,64), dimension=2, timestep=1.0, boundary=1, gravity=vec3(0, -0.01, 0), name="Scene", show_gui=True, pause_on_start=False, **kwargs):
        print("INIT base")
        self.resolution = resolution
        self.dimension = dimension
        self.timestep = timestep
        self.boundary = boundary
        self.gravity = gravity
        self.name = name
        self.show_gui = show_gui
        self.pause_on_start = pause_on_start

        self._parse_arguments(**kwargs)

        # previously _init_solver(self):
        self._grid_size = self.resolution
        if (self.dimension == 2):
            self._grid_size.z = 1

        self.solver = FluidSolver(name=self.name, gridSize=self._grid_size, dim=self.dimension)
        self.solver.timestep = self.timestep
        self.max_iter_fac = 10
        self.accuracy     = 5e-5

        # adaptive timestepping settings, not used at the moment
        self.solver.frameLength = self.timestep
        self.solver.timestepMin = self.timestep * 0.25 
        self.solver.timestepMax = self.timestep * 4.0  
        self.solver.cfl         = 2.0

        setDebugLevel(0)

        # previously _init_grids(self):
        self.flags = self.solver.create(FlagGrid, name="Flags")
        self.vel = self.solver.create(MACGrid, name="Velocity")
        self.pressure = self.solver.create(RealGrid, name="Pressure")
        self.phi_fluid = self.solver.create(LevelsetGrid, name="Fluid")
        self.phi_obs = self.solver.create(LevelsetGrid, name="Obstacles")
        #self.phi_sink = self.solver.create(LevelsetGrid, name="Sinks")
        #self.phi_source = self.solver.create(LevelsetGrid, name="Sources")
        self.fractions = self.solver.create(MACGrid, name="Fractions")
        
        if self.show_gui and self.dimension > 2:
            self.debugmesh = self.solver.create(Mesh)
        self.flags.initDomain(boundaryWidth=self.boundary, phiWalls=self.phi_obs)

    #----------------------------------------------------------------------------------
    # Interface
    #----------------------------------------------------------------------------------
    def _parse_arguments(self, **kwargs):
        """ Parse additional arguments and store them in internal variables """
        pass

    #----------------------------------------------------------------------------------
    def _create_scene(self):
        """ Reset and initialize randomized initial conditions for simulation, override in derived classes """

        # TODO, should be done in scene init
        self.flags.setConst(0)
        self.vel.setConst(vec3(0,0,0))
        self.pressure.setConst(0)
        self.phi_fluid.setConst(0)
        self.phi_obs.setConst(0)
        #self.phi_sink.setConst(0)
        #self.phi_source.setConst(0)
        self.fractions.setConst(vec3(0,0,0))
        self.flags.initDomain(boundaryWidth=self.boundary, phiWalls=self.phi_obs)
        
        # TODO potentially dangerous - make optional?
        self.phi_fluid.subtract(self.phi_obs)
        updateFractions(flags=self.flags, phiObs=self.phi_obs, fractions=self.fractions, boundaryWidth=self.boundary)
        setObstacleFlags(flags=self.flags, phiObs=self.phi_obs, fractions=self.fractions)
        
        extrapolateLsSimple(phi=self.phi_fluid, distance=5)
        self.flags.updateFromLevelset(self.phi_fluid)

        if self.show_gui and "_gui" not in self.__dict__:
            self._gui = Gui()
            self._gui.show(True)
            self._gui.setCamRot(40.0,0,0)
            self._gui.setCamPos(0,0,-1.5)
            self._gui.setPlane(2) 
            if self.pause_on_start:
                self._gui.pause()

    #----------------------------------------------------------------------------------
    def _compute_simulation_step(self):
        """ Compute simulation step for simulations, override in derived classes """
        # ... implement simulation step ...
        pass

    #----------------------------------------------------------------------------------
    # Functions
    #----------------------------------------------------------------------------------
#    def get_settings(self):
#        settings = {
#            "resolution": self.resolution,
#            "dimension": self.dimension,
#            "timestep": self.timestep,
#            "boundary": self.boundary,
#            "name": self.name
#        }
#        return settings
#
#     #----------------------------------------------------------------------------------
#     def add_obstacle(self, volume):
#         vol_ls = volume.computeLevelset(self.solver)
#         self.phi_obs.join(vol_ls)
# 
#     #----------------------------------------------------------------------------------
#     def add_fluid(self, volume):
#         vol_ls = volume.computeLevelset(self.solver)
#         self.phi_fluid.join(vol_ls)
# 
#     #----------------------------------------------------------------------------------
#     def set_velocity(self, volume, velocity):
#         if self.dimension == 2:
#             velocity.z = 0.0
#         volume.applyToGrid(solver=self.solver, grid=self.vel, value=velocity)
# 
#     #----------------------------------------------------------------------------------
#     def add_sink(self, volume):
#         vol_ls = volume.computeLevelset(self.solver)
#         self.phi_sink.join(vol_ls)
# 
#     #----------------------------------------------------------------------------------
#     def add_source(self, volume):
#         vol_ls = volume.computeLevelset(self.solver)
#         self.phi_source.join(vol_ls)

    #----------------------------------------------------------------------------------
    def simulate(self, num_steps=100, on_simulation_step=None):
        self._create_scene()

        self.solver.frame = 0 
        last_frame = -1
        while self.solver.frame < num_steps:
            # TODO, make optional via derived class?
            #maxVel = self.vel.getMax()
            #self.solver.adaptTimestep( maxVel )
            print("\r{} Step {:3d}, Time {:3.3f}, dt {:0.3f}".format(self.name, self.solver.frame + 1, self.solver.timeTotal, self.solver.timestep), end='\r')

            # Execute simulation step
            self._compute_simulation_step()

            # Update GUI
            if self.show_gui and self.dimension > 2:
                #self.phi_fluid.setBound(0.5, 0) # optionally, close sides for display
                self.phi_fluid.createMesh(self.debugmesh)

            # Advance solver and call callback function
            self.solver.step()
            if callable(on_simulation_step) and (last_frame != self.solver.frame):
                assert on_simulation_step.__code__.co_argcount == 2, "on_simulation_step must be a function with 2 arguments (scene and timestep)!"
                on_simulation_step(self, self.solver.frame-1) # solver already progressed one frame
            last_frame = self.solver.frame




#=========================================================================================================================
# TODO, NaiveScene NT_DEBUG update?
class NaiveScene_dontUse_(Scene):
    #----------------------------------------------------------------------------------
    def _init_solver(self):
        super(NaiveScene,self)._init_solver()

    #----------------------------------------------------------------------------------
    def _parse_arguments(self, **kwargs):
        self.merge_ghost_fluid = kwargs.get("merge_ghost_fluid", False) 

    #----------------------------------------------------------------------------------
    def _compute_simulation_step(self):
        self._advect(ls_order=1)
        self._enforce_boundaries(distance=2)
        self._solve_pressure(max_iter_fac=self.max_iter_fac, accuracy=self.accuracy)

    #==================================================================================
    # SIMULATION
    #----------------------------------------------------------------------------------
    def _advect(self, extrapol_dist=3, ls_order=2):
        # extrapolate the grids into empty cells
        self.phi_fluid.reinitMarching(flags=self.flags, velTransport=self.vel, maxTime=32.0)

        # advect the levelset
        advectSemiLagrange(flags=self.flags, vel=self.vel, grid=self.phi_fluid, order=ls_order)

        # source & sink
        self.phi_fluid.subtract(self.phi_sink)
        #self.phi_fluid.join(self.phi_source)

        # enforce boundaries
        self.phi_fluid.setBoundNeumann(self.boundary)
        self.flags.updateFromLevelset(self.phi_fluid)

        # velocity self-advection
        advectSemiLagrange(flags=self.flags, vel=self.vel, grid=self.vel, order=2)
        addGravity(flags=self.flags, vel=self.vel, gravity=self.gravity)

    #----------------------------------------------------------------------------------
    def _enforce_boundaries(self, distance):
        # enforce boundaries
        setWallBcs(flags=self.flags, vel=self.vel, fractions=self.fractions, phiObs=self.phi_obs)

    #----------------------------------------------------------------------------------
    def _solve_pressure(self, max_iter_fac, accuracy):
        solvePressure(vel=self.vel, pressure=self.pressure, flags=self.flags, fractions=self.fractions, cgAccuracy=accuracy, cgMaxIterFac=max_iter_fac, phi=self.phi_fluid)

        if self.boundary>0:
            self.pressure.setBoundNeumann(self.boundary-1)

        self._enforce_boundaries(distance=4)
