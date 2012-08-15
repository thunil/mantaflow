#
# Simple example scene (hello world)
# Simulation of a buoyant smoke density plume

from manta import *

# solver params
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 1.0

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
phi = s.create(LevelsetGrid)
pressure = s.create(RealGrid)

# scene setup
flags.initDomain()
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.3,0.6,1))
fluidbox.computeLevelset(phi)
flags.updateFromLevelset(phi)

if (GUI):
    gui = Gui()
    gui.show()
    #bgr = s.create(Mesh)
    #bgr.fromShape(fluidbox)
    
#main loop
for t in range(200):
    
    phi.reinitMarching(flags)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
    flags.updateFromLevelset(phi)
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))
    
    setWallBcs(flags=flags, vel=vel)    
    setLiquidBcs(flags=flags, vel=vel)    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setLiquidBcs(flags=flags, vel=vel)    
    setWallBcs(flags=flags, vel=vel)
    
    s.step()
