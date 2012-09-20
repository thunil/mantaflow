#
# Filament test
# don't use yet --- work in progress

from manta import *

# solver params
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.25

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=1)
drop = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.15)
basin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,0.2,1))
phi = basin.computeLevelset()
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)

if (GUI):
    gui = Gui()
    gui.show()
    
#main loop
for t in range(200):
    
    # update and advect levelset
    phi.reinitMarching(flags=flags, velTransport=vel) #, ignoreWalls=False)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
    flags.updateFromLevelset(phi)
    
    # velocity self-advection
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.025,0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    setLiquidBcs(flags=flags, vel=vel)
    solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=0.005, useResNorm=False)
    setLiquidBcs(flags=flags, vel=vel)    
    setWallBcs(flags=flags, vel=vel)
    
    # note: these meshes are created by fast marching only, should smooth
    #       geometry and normals before rendering
    phi.createMesh(mesh)
    #mesh.save('phi%04d.bobj.gz' % t)
    
    s.printTimings()
    s.step()
    #gui.pause()
