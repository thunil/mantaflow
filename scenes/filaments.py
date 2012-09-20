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
filaments = s.create(VortexFilamentSystem)

# scene setup
flags.initDomain(boundaryWidth=1)
flags.fillGrid()

filaments.addRing(position=gs*vec3(0.5,0.1,0.5), circulation=1, radius=0.2*res, normal=(0,1,0), number=100)
filaments.addRing(position=gs*vec3(0.5,0.14,0.5), circulation=1, radius=0.2*res, normal=(0,1,0), number=100)

if (GUI):
    gui = Gui()
    gui.show()
    
#main loop
for t in range(2000):
    filaments.advectSelf(scale=1e-3, regularization=1, integrationMode=IntRK4)
    
    # velocity self-advection
    #advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    
    # pressure solve
    #setWallBcs(flags=flags, vel=vel)    
    #solvePressure(flags=flags, vel=vel, pressure=pressure)
    #setWallBcs(flags=flags, vel=vel)
    
    s.step()
    
