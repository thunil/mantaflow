#
# Simple example for free-surface simulation
# with MacCormack advection

from manta import *

# solver params
gs = vec3(96,40,40)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)
phi = s.create(LevelsetGrid)
vorty = s.create(RealGrid)

# scene setup
flags.initDomain(boundaryWidth=0)
step = s.create(Box, p0=gs*vec3(0.3,0,0), p1=gs*vec3(0.4,0.2,1))
step.applyToGrid(grid=flags, value=FlagObstacle)

h = 0.6
riverInit = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,h,1))
#riverInit.applyToGrid(grid=flags, value=FlagFluid, respectFlags=flags)
#riverInit.applyToGrid(grid=vel, value=(0.2,0,0), respectFlags=flags)
phi.initFromFlags(flags=flags, ignoreWalls = True)

if (GUI):
    gui = Gui()
    gui.show()
    gui.pause()
    
#main loop
for t in range(2000):
    
    # update and advect levelset
    phi.reinitMarching(flags=flags, velTransport=vel, ignoreWalls = True) #, ignoreWalls=False)
    getCurl(vel=vel, vort=vorty, comp=2)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
    flags.updateFromLevelset(phi)
    
    # velocity self-advection
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.005,0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    setLiquidBcs(flags=flags, vel=vel)
    setinflow(flags=flags, vel=vel,phi=phi,h=h)
    solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=1.5, cgAccuracy=0.001, useResNorm=False, phi=phi, ghostAccuracy=0.0) 
    setLiquidBcs(flags=flags, vel=vel)    
    setWallBcs(flags=flags, vel=vel)
    setinflow(flags=flags, phi=phi, vel=vel,h=h)
        
    # note: these meshes are created by fast marching only, should smooth
    #       geometry and normals before rendering
    phi.createMesh(mesh)
    
    s.step()
    #gui.pause()
