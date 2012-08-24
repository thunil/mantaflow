#
# Simple example for free-surface simulation
# with FLIP advection

from manta import *

# solver params
meshing = False
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5

# prepare grids and particles
flags = s.create(FlagGrid3)
vel = s.create(MACGrid3)
phi = s.create(LevelsetGrid3)
pressure = s.create(RealGrid3)
flip = s.create(FlipSystem)
mesh = s.create(Mesh)

# scene setup
flags.initDomain()
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.8,1))
fluidbox.computeLevelset(phi)
flags.updateFromLevelset(phi)
flip.adjustNumber(vel=vel, flags=flags, minParticles=8, maxParticles=30)
    
if (GUI):
    gui = Gui()
    gui.show()
    
#main loop
for t in range(200):
    
    # FLIP advect and writeback
    flip.advectInGrid(flaggrid=flags, vel=vel, integrationMode=IntRK4)
    flip.velocitiesToGrid(vel=vel, flags=flags)
    flip.adjustNumber(vel=vel, flags=flags, minParticles=4, maxParticles=30)
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
        
    # update and advect levelset
    phi.reinitMarching(flags=flags, ignoreWalls=False)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
    flags.updateFromLevelset(phi)
    
    # note: these meshes are created by fast marching only, should smooth
    #       geometry and normals before rendering
    if (meshing):
        phi.createMesh(mesh)
        mesh.save('phi%04d.bobj.gz' % t)
    
    s.step()
