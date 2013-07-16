#
# Simple example scene (hello world)
# Simulation of a buoyant smoke density plume

from manta import *
import os, shutil, math, sys

# dimension two/three d
dim = 2
# how much to upres the XL sim?
upres = 2

# solver params
res = 120
gs = vec3(res,int(1.5*res),res)
# gs = vec3(res,res,res)
if (dim==2): gs.z = 1  # 2D

sm = Solver(name='main', gridSize = gs, dim=dim)


# sm.timestep = 1.0
sm.timestep = 1.5

velInflow = vec3(2, 0, 0)
#velInflow = vec3(1, 1, 1)

# prepare grids
flags    = sm.create(FlagGrid)
vel      = sm.create(MACGrid)
density  = sm.create(RealGrid)
pressure = sm.create(RealGrid)
energy   = sm.create(RealGrid)

# inflow noise field
noise = sm.create(NoiseField, fixedSeed=265, loadFromFile=True)
noise.posScale = vec3(20)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 2
noise.valScale = 1
noise.valOffset = 0.075
noise.timeAnim = 0.6

flags.initDomain()
flags.fillGrid()

source = sm.create(Cylinder, center=gs*vec3(0.3,0.2,0.5), radius=res*0.081, z=gs*vec3(0.081, 0, 0))
sourceVel = sm.create(Cylinder, center=gs*vec3(0.3,0.2,0.5), radius=res*0.15, z=gs*vec3(0.15, 0, 0))


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

xl_source = xl.create(Cylinder, center=xl_gs*vec3(0.3,0.2,0.5), radius=xl_gs.x*0.081, z=xl_gs*vec3(0.081, 0, 0))

xl_noise = xl.create(NoiseField, fixedSeed=265, loadFromFile=True)
xl_noise.posScale = vec3(20)
xl_noise.clamp = noise.clamp
xl_noise.clampNeg = noise.clampNeg
xl_noise.clampPos = noise.clampPos
xl_noise.valScale = noise.valScale
xl_noise.valOffset = noise.valOffset
xl_noise.timeAnim  = noise.timeAnim

# wavelet turbulence noise
wltnoise = sm.create(NoiseField, loadFromFile=True)
# scale according to lowres sim , smaller numbers mean larger vortices
wltnoise.posScale = vec3( int(0.5*gs.x) )
wltnoise.timeAnim = 0.1

wltnoise2 = sm.create(NoiseField, loadFromFile=True)
wltnoise2.posScale = vec3( int(1*gs.x) )
wltnoise2.timeAnim = 0.1



if (GUI):
    gui = Gui()
    gui.show()

#main loop
for t in range(150):
    
    curt = t * sm.timestep
    #sys.stdout.write( "Curr t " + str(curt) +" \n" )
        
    #source.applyToGrid(grid=vel, value=velInflow)
    advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    
    applyInflow=False
    if (curt>=0 and curt<75):
        densityInflow( flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5 )
        sourceVel.applyToGrid( grid=vel , value=velInflow )
        applyInflow=True
    
    setWallBcs(flags=flags, vel=vel)    
    addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-3,0), flags=flags)
    
    #applyNoiseVec3( flags=flags, target=vel, noise=noise, scale=1 ) # just to test, add everywhere...
    
    solvePressure(flags=flags, vel=vel, pressure=pressure , openBound='Y', \
        cgMaxIterFac=0.5, cgAccuracy=0.01 )
    setWallBcs(flags=flags, vel=vel)
    
    computeEnergy(flags=flags, vel=vel, energy=energy)
    computeWaveletCoeffs(energy)
    #computeVorticity( vel=vel, vorticity=vort, norm=energy);
    #computeStrainRateMag( vel=vel, vorticity=vort, mag=energy);
    
    #density.save('densitySm_%04d.vol' % t)
    
    sm.step()
    
    # xl ...
    # same inflow
    
    interpolateGrid( target=xl_weight, source=energy )
    interpolateMACGrid( source=vel, target=xl_vel )
    
    applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=wltnoise, scale=1. , weight=xl_weight)
    # manually apply second octave for now
    applyNoiseVec3( flags=xl_flags, target=xl_vel, noise=wltnoise2, scale=0.6 , weight=xl_weight)
    
    for substep in range(upres):
        advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_density, order=2)    
    # advectSemiLagrange(flags=xl_flags, vel=xl_vel, grid=xl_vel, order=2)
    
    if (applyInflow):
         densityInflow( flags=xl_flags, density=xl_density, noise=xl_noise, shape=xl_source, scale=1, sigma=0.5 )
         # source.applyToGrid( grid=xl_vel , value=velInflow )
    
    #xl_density.save('densityXl08_%04d.vol' % t)
    
    xl.step()    
    gui.screenshot( 'densXl_wlt_%04d.png' % t );

