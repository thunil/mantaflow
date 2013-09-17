#
# Simple example for free-surface simulation, falling drop
# that shouldnt hit the floor and give straight down velocities

from manta import *
from helperInclude import *

# solver params
dim = 3
res = 45
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.6
accuracy = 5e-5

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=0)
liqDrop = s.create(Box, p0=gs*vec3(0.4,0.75,0.4), p1=gs*vec3(0.6,0.95,0.6))
phi = liqDrop.computeLevelset()
flags.updateFromLevelset(phi)
        
if 0 and (GUI):
    gui = Gui()
    gui.show()
    #gui.pause()
    
#main loop
for t in range(18):
    
    # update and advect levelset
    phi.reinitMarching(flags=flags, velTransport=vel) #, ignoreWalls=False)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
    flags.updateFromLevelset(phi)
    
    # velocity self-advection
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.0125,0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)
    solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, useResNorm=False, phi=phi)
    setWallBcs(flags=flags, vel=vel)
    
    if (dim==3):
        phi.createMesh(mesh)
    
    s.step()
	#gui.screenshot( 'screenOn_%04d.png' % t );


# check final state
doTestGrid( __file__,"phi"  , s, phi    )
doTestGrid( __file__,"vel"  , s, vel    )

