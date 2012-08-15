#
# Simple example for free-surface simulation
# with FLIP advection

from manta import *

# solver params
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.3

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
phi = s.create(LevelsetGrid)
pressure = s.create(RealGrid)
flip = s.create(FlipSystem)

# scene setup
flags.initDomain()
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.3,0.6,1))
fluidbox.computeLevelset(phi)
flags.updateFromLevelset(phi)

if (GUI):
    gui = Gui()
    gui.show()
    
#main loop
for t in range(200):
    
    # update and advect levelset
    phi.reinitMarching(flags)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
    flags.updateFromLevelset(phi)
    
    # FLIP advect and writeback
    flip.advectInGrid(flaggrid=flags, vel=vel, integrationMode=IntRK4)
    flip.velocitiesToGrid(vel=vel, flags=flags)
    flip.adjustNumber(vel=vel, flags=flags, minParticles=6, maxParticles=12)
    #advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    
    addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    setLiquidBcs(flags=flags, vel=vel)
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setLiquidBcs(flags=flags, vel=vel)    
    setWallBcs(flags=flags, vel=vel)
    
    # FLIP load
    flip.velocitiesFromGrid(vel=vel, flags=flags, flipRatio=0)
        
    s.step()
