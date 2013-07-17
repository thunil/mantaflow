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

k = s.create(RealGrid)
eps = s.create(RealGrid)
prod = s.create(RealGrid)
nuT= s.create(RealGrid)
strain= s.create(RealGrid)
vc=s.create(MACGrid)
temp=s.create(RealGrid)
vel = s.create(MACGrid)

# noise field
noise = s.create(NoiseField)
noise.timeAnim = 0
#noise.posScale = vec3(1)

# turbulence particles
#turb = s.create(TurbulenceParticleSystem, noise=noise)

flags.initDomain()
flags.fillGrid()

# obstacle grid
for i in range(4):
    for j in range(4):
        obs = s.create(Sphere, center=gs*vec3(0.2,(i+1)/5.0,(j+1)/5.0), radius=res*0.025)
        obs.applyToGrid(grid=flags,value=FlagObstacle)

# particle inflow
cyl = s.create(Cylinder, center = gs*vec3(0.07,0.5,0.5), radius=res*0.1, z=vec3(res*0.03,0,0))

if (GUI):
    gui = Gui()
    gui.show()

KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=0.1,nu=0.1,fillArea=True)

#box=s.create(Box, p0=vec3(0,0,0), p1=gs*vec3(1,1,1))
#turb.seed(box,10000)

#main loop
for t in range(10000):
    #turb.seed(cyl,20)
    #turb.advectInGrid(flaggrid=flags, vel=vel, integrationMode=IntRK4)
    #turb.synthesize(flags=flags, octaves=1, k=k, switchLength=10, L0=0.001, scale=1e-1)
    applyK41(grid=vc,noise=noise, L0=0.1, scale=1,octaves=1)
    
    KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=0.1,nu=0.1,fillArea=False)
    advectSemiLagrange(flags=flags, vel=vel, grid=k, order=2)
    advectSemiLagrange(flags=flags, vel=vel, grid=eps, order=2)
    KEpsilonBcs(flags=flags,k=k,eps=eps,intensity=0.1,nu=0.1,fillArea=False)
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
    