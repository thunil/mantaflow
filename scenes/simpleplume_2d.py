#
# Simple example scene for a 2D simulation
# Simulation of a buoyant smoke density plume

from manta import *

# solver params
res = 64
gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=2)
s.timestep = 1.0

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field
#noise = s.create(NoiseField)
#noise.posScale = vec3(45)
#noise.clamp = True
#noise.clampNeg = 0
#noise.clampPos = 1
#noise.valScale = 1
#noise.valOffset = 0.75
#noise.timeAnim = 0.2

flags.initDomain()
flags.fillGrid()

if (GUI):
    gui = Gui()
    gui.show()

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0), radius=res*0.14, z=gs*vec3(0, 0.02, 0))
    
#main loop
for t in range(2000):
    if t<100:
        #densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
        source.applyToGrid(grid=density, value=1)
        
    #source.applyToGrid(grid=vel, value=velInflow)
    advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    
    setWallBcs(flags=flags, vel=vel)    
    addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-4,0), flags=flags)
    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)
    
    s.step()
