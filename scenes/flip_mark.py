#
# Simple example for free-surface simulation
# with FLIP advection and simple cell marking

from manta import *

# solver params
meshing = True
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
flip = s.create(FlipSystem)
mesh = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=1)
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1))
fluidbox.applyToGrid(grid=flags, value=FlagFluid, respectFlags=flags)
flip.adjustNumber(vel=vel, flags=flags, minParticles=8, maxParticles=30)
    
if (GUI):
    gui = Gui()
    gui.show()
    gui.pause()
    
#main loop
for t in range(200):
    
    # FLIP advect and writeback
    flip.advectInGrid(flaggrid=flags, vel=vel, integrationMode=IntRK4)
    flip.velocitiesToGrid(vel=vel, flags=flags)
    flip.markFluidCells(flags=flags)
    #advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    
    addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    setLiquidBcs(flags=flags, vel=vel)
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setLiquidBcs(flags=flags, vel=vel)    
    setWallBcs(flags=flags, vel=vel)
    
    # FLIP load
    flip.velocitiesFromGrid(vel=vel, flags=flags, flipRatio=0.96)
        
    s.step()
