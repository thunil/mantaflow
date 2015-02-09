#
# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

# solver params
dim = 2
particleNumber = 2
res = 64
gs = vec3(res,res,res)
if (dim==2):
    gs.z=1
    particleNumber = 3 # use more particles in 2d
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5

usePhi = True

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 
phi     = s.create(LevelsetGrid) 
if usePhi:
    pindex = s.create(ParticleIndexSystem) 
    gpi    = s.create(IntGrid) 

# scene setup
flags.initDomain(boundaryWidth=2) # MLE: should not use default 0, leads to asymmetry and other inaccuracies
# enable one of the following
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
#fluidbox = s.create(Box, p0=gs*vec3(0.4,0.72,0.4), p1=gs*vec3(0.6,0.92,0.6)) # centered falling block
phi = fluidbox.computeLevelset()
flags.updateFromLevelset(phi)

setOpenBound(flags,'y')

sampleLevelsetWithParticles(phi=phi,flags=flags,parts=pp,discretization=particleNumber,randomness=0.05) if usePhi else sampleFlagsWithParticles( flags=flags, parts=pp, discretization=particleNumber, randomness=0.2 )
    
if (GUI):
    gui = Gui()
    gui.show()
    #gui.pause()
    
#main loop
for t in range(2500):
    
    # FLIP 
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
    resetOpenBoundFLIP(flags,pp) 
    mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
    extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
    markFluidCells( parts=pp, flags=flags )

    if usePhi:
        gridParticleIndex(parts=pp,flags=flags,indexSys=pindex,index=gpi)
        unionParticleLevelset(pp,pindex,flags,gpi,phi,1.0) 
        # set source grids for resampling, used in adjustNumber!
        pVel.setSource(vel,isMAC=True)
        adjustNumber(parts=pp,vel=vel,flags=flags,minParticles=1*pow(2,dim),maxParticles=2*pow(2,dim),phi=phi,radiusFactor=1.0) 
        
    addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))

    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi) if usePhi else solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)

    # we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
    extrapolateMACSimple( flags=flags, vel=vel )
    
    # FLIP velocity update
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )
    
    #gui.screenshot( 'flipt_%04d.png' % t );
    s.step()

