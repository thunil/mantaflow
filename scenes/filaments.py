#
# Filament test
# don't use yet --- work in progress

from manta import *
from math import *
import random

# solver params
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.1

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)
filaments = s.create(VortexFilamentSystem)
tracer = s.create(TracerParticleSystem)

# scene setup
flags.initDomain(boundaryWidth=1)
flags.fillGrid()

#filaments.addRing(position=gs*vec3(0.5,0.1,0.5), circulation=10, radius=0.2*res, normal=(0,1,0), number=40)
#filaments.addRing(position=gs*vec3(0.5,0.14,0.5), circulation=10, radius=0.2*res, normal=(0,1,0), number=1000)
#filaments.remesh(maxLen=3, minLen=1)
    
#filaments.setPos(10,filaments.getPos(10)+vec3(0,1.5,0))
source = s.create(Cylinder, center=gs*vec3(0.5,0.09,0.5), radius=res*0.17, z=gs*vec3(0, 0.03, 0))
fixedRegion = s.create(Box, center=gs*vec3(0.5,0.07,0.5), size=gs*vec3(0.4,0.03,0.4))
mesh.fromShape(source)

if (GUI):
    gui = Gui()
    gui.show()
    gui.pause()

r = 0.2
scale = 1
    
#main loop
for t in range(1000):
    markAsFixed(mesh=mesh, shape=fixedRegion)

    # seed rings
    if random.random()<0.05:
        theta = random.gauss(0,0.01*pi)
        phi = random.gauss(0,0.01*pi)
        n = vec3(sin(theta)*cos(phi), cos(theta), sin(theta)*sin(phi))
        filaments.addRing(position=gs*vec3(0.5,0.14,0.5), circulation=random.uniform(10,20), radius=res*random.gauss(r,0.01), normal=n, number=20)
    
    # seed tracer particles
    #for i in range(10):
    #    x=100
    #    z=100
    #    while x*x+z*z > r*r :
    #        x=random.uniform(-r,r)
    #        y=random.uniform(0,0.1)
    #        z=random.uniform(-r,r)
    #    tracer.addParticle(vec3(x+0.5,y+0.1,z+0.5)*gs)
    
    # mesh cosmetics
    smoothMesh(mesh=mesh, strength=1e-4, steps=1)
    subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=2*scale, cutTubes=True)
    killSmallComponents(mesh=mesh, elements=20)
        
    filaments.remesh(maxLen=3, minLen=1)
    #gui.pause()
    #filaments.ddTest(d=0.25,phi=0.1)
    #filaments.doublyDiscreteUpdate(regularization=2)
    filaments.advectSelf(scale=1, regularization=1, integrationMode=IntRK4)
    #filaments.advectParticles(sys=tracer, scale=1, regularization=1, integrationMode=IntRK2)
    filaments.advectMesh(mesh=mesh, scale=1, regularization=6, integrationMode=IntRK2)
    
    s.printTimings()
    s.step()
    
