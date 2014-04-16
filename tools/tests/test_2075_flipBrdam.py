#
# More complex flip setup, breaking dam with resampling and
# additional particle data fields
# 
import sys
from manta import *
from helperInclude import *

# solver params
dim = 3
res = 52
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.75
minParticles = pow(2,dim)

# use slightly larger radius to make the sim a bit harder
radiusFactor = 1.5

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
tstGrid  = s.create(RealGrid)
dens     = s.create(RealGrid)
dens2    = s.create(RealGrid)

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
# test real value, not necessary for simulation
pInt     = pp.create(PdataInt) # todo, no effect so far
pDens    = pp.create(PdataReal) 
pDens2   = pp.create(PdataReal) # read density back in from grid
#mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup
flags.initDomain(boundaryWidth=0)
# enable one of the following
fluidbox = s.create(Box, p0=gs*vec3(0.6,0.2,0.1), p1=gs*vec3(0.95,0.7,0.8)) # breaking dam, asymm.
phi = fluidbox.computeLevelset()
flags.updateFromLevelset(phi)

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.2 )

# testing the real channel while resampling - original particles
# will have a value of 0.1, new particle will get a value from the tstGrid
testInitGridWithPos(tstGrid)
setConstPdata( pDens , 0.1 )
setConstPdata( pDens2, 0.8 ) # should be overwritten
    
if 0 and (GUI):
    gui = Gui()
    gui.show()
    #gui.pause()
    
#main loop
for t in range(25):
    
    # FLIP  ,  as a test, delete particles in obstacles here
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=True )

    # make sure we have velocities throught liquid region
    mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
    extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
    markFluidCells( parts=pp, flags=flags )

	# create approximate surface level set, resample particles
    gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
    unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) 
    phi.reinitMarching(flags=flags, maxTime=int(2*radiusFactor) )
    pVel.setSource( vel, isMAC=True )
    pDens.setSource( tstGrid );
    adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

    mapPartsToGrid(target=dens,  flags=flags, parts=pp, source=pDens )

    mapGridToParts(source=dens,               parts=pp, target=pDens2 )
    # to check the result, we have to write to a grid again...
    mapPartsToGrid(target=dens2, flags=flags, parts=pp, source=pDens2 )

	# forces & pressure solve
    addGravity(flags=flags, vel=vel, gravity=(0,-0.01,0))
    setWallBcs(flags=flags, vel=vel)    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)

    # make sure we have proper velocities
    extrapolateMACSimple( flags=flags, vel=vel )
    
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

    #if (dim==3):
        #phi.createMesh(mesh)
    
    s.step()

doTestGrid( sys.argv[0],"dens" , s, dens   )
doTestGrid( sys.argv[0],"dens2", s, dens2  )
doTestGrid( sys.argv[0],"vel"  , s, vel    )


