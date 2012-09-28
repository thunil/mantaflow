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

filaments.addRing(position=gs*vec3(0.5,0.1,0.5), circulation=50, radius=15, normal=(0,1,0), number=128)
#filaments.addRing(position=gs*vec3(0.5,0.14,0.5), circulation=10, radius=0.2*res, normal=(0,1,0), number=120)

if (GUI):
    gui = Gui()
    gui.show()
#    gui.pause()
    
#main loop
for t in range(10):
    #filaments.doublyDiscreteUpdate(regularization=0.1)
    filaments.advectSelf(scale=1, regularization=0.1, integrationMode=IntRK4)
    filaments.remesh(maxLen=1)
    
    s.step()
    
