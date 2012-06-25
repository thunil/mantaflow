from manta import *
import os, shutil, math

scene='plume'
res = 64
scale = 0.2

# solver params
vc = vec3(1, 1, 1) * 1.5 * res
gs = vec3(1,1.5,1) * res
dx = 1/(1.5*res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5

gravity = vec3(0,-0.05,0)*3
velInflow = vec3(0,0.26,0)*2

# prepare grids
flags = s.create(FlagGrid)
pressure = s.create(RealGrid, show=False)
vel = s.create(MACGrid)
density = s.create(RealGrid)
densityInflow = s.create(RealGrid, show=False)

# noise field
noise = s.create(NoiseField)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

mesh = s.create(VortexSheetMesh)
backgr = s.create(Mesh)

vp = s.create(VortexParticleSystem)

flags.initDomain()
flags.fillGrid()

if (GUI):
    gui = Gui()
    gui.show()
    gui.setBackgroundMesh(backgr)

source = s.create(Cylinder, center=gs*vec3(0.5,0.13,0.5), radius=res*0.1*1.4, z=gs*vec3(0, 0.03, 0))
mesh.fromShape(source)
subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=2*scale)

fixedRegion = s.create(Box, center=gs*vec3(0.5,0.09,0.5), size=gs*vec3(0.4,0.03,0.4))
    
#main loop
for t in range(250):
    densityNoiseInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
    
    source.applyToGrid(grid=vel, value=velInflow)
    advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)
    
    markAsFixed(mesh=mesh, shape=fixedRegion)

    mesh.calcCirculation()        
    mesh.advectInGrid(vel=vel, flaggrid=flags, integrationMode=IntRK4)  
    mesh.calcVorticity()
    
    #smoothVorticity(mesh=mesh, iter=50, alpha=1.0, sigma=0.6)        
    filterVorticityCuda(mesh=mesh, sigma=1.6)
    meshApplyBuoyancyLocalCuda(mesh=mesh, scale=0.6*1e-2, cutoffCells=5, regularization=0.5)
    mesh.calcVorticity()
    
    smoothMesh(mesh=mesh, strength=1e-4, steps=1)
    subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=2*scale, cutTubes=True)
    killSmallComponents(mesh=mesh, elements=20)
    
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    setWallBcs(flags=flags, vel=vel)
    
    addBuoyancy(density=density, vel=vel, gravity=gravity*dx, flags=flags)
    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)
    
    vorticitySource(mesh=mesh, gravity=gravity, scale=0.1, maxAmount=200)
    
    mesh.save('/home/tpfaff/tmp/plume/e%04d.bobj.gz' % t)
    s.printTimings()
    s.step()
    