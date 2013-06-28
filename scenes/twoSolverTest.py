#
# Simple example scene (hello world)
# Simulation of a buoyant smoke density plume

from manta import *
import os, shutil, math, sys

# dimension two/three d
dim = 2
# how much to upres the XL sim?
upres = 3

# solver params
res = 60
gs = vec3(res,int(1.5*res),res)
# gs = vec3(res,res,res)
if (dim==2): gs.z = 1  # 2D

sm = Solver(name='main', gridSize = gs, dim=dim)


# sm.timestep = 1.0
sm.timestep = 2.0

velInflow = vec3(1, 0, 0)
#velInflow = vec3(1, 1, 1)

# prepare grids
flags    = sm.create(FlagGrid)
vel      = sm.create(MACGrid)
density  = sm.create(RealGrid)
pressure = sm.create(RealGrid)
energy   = sm.create(RealGrid)

# inflow noise field
noise = sm.create(NoiseField, fixedSeed=765)
noise.posScale = vec3(20)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 2
noise.valScale = 1
noise.valOffset = 0.075
noise.timeAnim = 0.2

flags.initDomain()
flags.fillGrid()

source = sm.create(Cylinder, center=gs*vec3(0.2,0.2,0.5), radius=res*0.1, z=gs*vec3(0.1, 0, 0))


# larger solver, recompute sizes...

xl_gs = vec3(upres*gs.x,upres*gs.y,upres*gs.z)
if (dim==2): xl_gs.z = 1  # 2D
xl = Solver(name='larger', gridSize = xl_gs, dim=dim)
xl.timestep = upres*sm.timestep
xl.timestep = sm.timestep

xl_flags   = xl.create(FlagGrid)
xl_vel     = xl.create(MACGrid)
xl_density = xl.create(RealGrid)
xl_weight  = xl.create(RealGrid)

xl_flags.initDomain()
xl_flags.fillGrid()

xl_source = xl.create(Cylinder, center=xl_gs*vec3(0.2,0.2,0.5), radius=xl_gs.x*0.1, z=xl_gs*vec3(0.1, 0, 0))

xl_noise = xl.create(NoiseField, fixedSeed=765)
xl_noise.posScale = vec3(20)
xl_noise.clamp = noise.clamp
xl_noise.clampNeg = noise.clampNeg
xl_noise.clampPos = noise.clampPos
xl_noise.valScale = noise.valScale
xl_noise.valOffset = noise.valOffset
xl_noise.timeAnim  = noise.timeAnim

# wavelet turbulence noise
wltnoise = sm.create(NoiseField)
wltnoise.posScale = vec3(gs.x) # scale according to lowres sim
wltnoise.timeAnim = 0.1



if (GUI):
    gui = Gui()
    gui.show()

#main loop
for t in range(200):
    
    curt = t * sm.timestep
    #sys.stdout.write( "Curr t " + str(curt) +" \n" )
        
    #source.applyToGrid(grid=vel, value=velInflow)
    advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    
    if (curt>=0 and curt<30) or (curt>60 and curt<90):
        densityInflow( flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5 )
        source.applyToGrid( grid=vel , value=velInflow )
    
    setWallBcs(flags=flags, vel=vel)    
    addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-3,0), flags=flags)
    
    #applyNoiseVec3( flags=flags, target=vel, noise=noise, scale=1 ) # just to test, add everywhere...
    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)
    
    computeEnergy(flags=flags, vel=vel, energy=energy)
    computeWaveletCoeffs(energy)
    
#    density.save('densitySm_%04d.vol' % t)
    
    sm.step()
    
    # xl ...
    # same inflow
    
    interpolateGrid( target=xl_weight, source=energy )
    interpolateMACGrid( source=vel, target=xl_vel )
    
    applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=wltnoise, scale=2.5, weight=xl_weight)
    
    for substep in range(upres):
        advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=2)    
    # advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_vel, order=2)
    
    if (curt>=0 and curt<30) or (curt>60 and curt<90):
         densityInflow( flags=xl_flags, density=xl_density, noise=xl_noise, shape=xl_source, scale=1, sigma=0.5 )
         # source.applyToGrid( grid=xl_vel , value=velInflow )
    
#    xl_density.save('densityXl_%04d.vol' % t)
    
    xl.step()

