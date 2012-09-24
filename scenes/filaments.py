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

for i in range(20):
    filaments.addRing(position=gs*vec3(0.5,0.1+i*0.03,0.5), circulation=4, radius=0.2*res, normal=(0,1,0), number=30)

if (GUI):
    gui = Gui()
    gui.show()
    
#main loop
for t in range(200):
    filaments.advectSelf(scale=1, regularization=1, integrationMode=IntRK4)
    
    s.step()
    
#main loop
for t in range(200):
    filaments.advectSelf(scale=-1, regularization=1, integrationMode=IntRK4)
    
    s.step()
    
