#
# Example scene for mesh-based smoke simulations using Vortex Sheets
# ( http://graphics.ethz.ch/publications/papers/paperPfa12.php )
#
# This example requires CUDA to run

from manta import *

scale = 0.2

# solver params
res = 64
gs = vec3(res,res/2,res/2)
s = Solver(name='main', gridSize = gs)
s.timestep = 1

velInflow = vec3(0.52,0,0)

# prepare grids
flags = s.create(FlagGrid)
pressure = s.create(RealGrid, show=False)
vel = s.create(MACGrid)
density = s.create(RealGrid, show=False)

k = s.create(RealGrid)
eps = s.create(RealGrid)
prod = s.create(RealGrid)
nuT= s.create(RealGrid)
strain= s.create(RealGrid)
vc=s.create(MACGrid)

# noise field
noise = s.create(NoiseField)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

flags.initDomain()
flags.fillGrid()

for i in range(4):
    for j in range(4):
        obs = s.create(Sphere, center=gs*vec3(0.2,(i+1)/5.0,(j+1)/5.0), radius=res*0.025)
        obs.applyToGrid(grid=flags,value=FlagObstacle)

if (GUI):
    gui = Gui()
    gui.show()

KEpsilonInit(flags=flags,k=k,eps=eps,intensity=0.1,nu=0.1,fillArea=True)

#main loop
for t in range(10000):
    # seed smoke within the source region and apply inflow velocity condition
    #densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
    #source.applyToGrid(grid=vel, value=velInflow)
    #advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)
    
    KEpsilonInit(flags=flags,k=k,eps=eps,intensity=0.1,nu=0.1,fillArea=False)
    advectSemiLagrange(flags=flags, vel=vel, grid=k, order=1)
    advectSemiLagrange(flags=flags, vel=vel, grid=eps, order=1)
    KEpsilonInit(flags=flags,k=k,eps=eps,intensity=0.1,nu=0.1,fillArea=False)
    KEpsilonComputeProduction(vel=vel, k=k, eps=eps, prod=prod, nuT=nuT, strain=strain, pscale=1.0) 
    KEpsilonSources(k=k, eps=eps, prod=prod)
    
    KEpsilonGradientDiffusion(k=k, eps=eps, vel=vel, nuT=nuT);

    # base solver
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    setWallBcs(flags=flags, vel=vel)
    setInflowBcs(vel=vel,dir='xXyYzZ',value=velInflow)
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)
    setInflowBcs(vel=vel,dir='xXyYzZ',value=velInflow)
    
    s.step()
    