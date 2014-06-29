#
# Flip scene with adaptive time stepping (otherwise similar to flip02)
# 
from manta import *

# how many frames to calculate 
frames    = 200
# length of one frame (in "world time")
frameTime = 0.60

# maximal velocity per cell, adaptDt variables
cflFac = 1.0
dtMax  = 2.00
dtMin  = 0.10
lockDt = False
maxvel = 0
#cflFac = 999.0 # use this to turn adaptive time steps off, use dtMax

# solver params
dim = 3
res = 50
#res = 100
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
minParticles = pow(2,dim)
timings = Timings()

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
tstGrid  = s.create(RealGrid)
phiObs   = s.create(LevelsetGrid)

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
# test real value, not necessary for simulation
pTest    = pp.create(PdataReal) 
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup, 0=breaking dam, 1=drop into pool
# geometry in world units (to be converted to grid space upon init)
setup = 0
flags.initDomain(boundaryWidth=0)
fluidVel = 0
fluidSetVel = 0

if setup==0:
	# breaking dam
	fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
	phi = fluidbox.computeLevelset()
	#fluidbox2 = s.create(Box, p0=gs*vec3(0.2,0.7,0.3), p1=gs*vec3(0.5,0.8,0.6)) # breaking dam
	#phi.join( fluidbox2.computeLevelset() )
elif setup==1:
	# falling drop
	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.2,1.0)) # basin
	dropCenter = vec3(0.5,0.5,0.5)
	dropRadius = 0.15
	fluidSetVel= vec3(0,-0.03,0)
	fluidDrop  = s.create(Sphere, center=gs*dropCenter, radius=res*dropRadius)
	fluidVel   = s.create(Sphere, center=gs*dropCenter, radius=res*(dropRadius+0.05) )
	phi = fluidBasin.computeLevelset()
	phi.join( fluidDrop.computeLevelset() )

flags.updateFromLevelset(phi)

# obstacle init needs to go after updateFromLs
obsBox = s.create(Box, p0=gs*vec3(0.7,0.0,0.5), p1=gs*vec3(0.8,1.0,0.8)) 
obsBox.applyToGrid(grid=flags, value=(FlagObstacle) )
#obsBox.applyToGrid(grid=flags, value=(FlagObstacle|FlagStick) )

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05 )
mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

if fluidVel!=0:
	# set initial velocity
	fluidVel.applyToGrid( grid=vel , value=gs*fluidSetVel )
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

# testing the real channel while resampling - original particles
# will have a value of 0.1, new particle will get a value from the tstGrid
testInitGridWithPos(tstGrid)
setConstPdata( pTest , 0.1 )

if 1 and (GUI):
	gui = Gui()
	gui.show( dim==2 )
	gui.pause()
	  
	# show all particles shaded by velocity
	gui.nextPdata()
	gui.nextPartDisplay()
	gui.nextPartDisplay()


# recompute time step
def adaptDt( solver ):
	global lockDt, dtMin, dtMax, maxvel

	if not lockDt:
		maxvel = (vel.getMaxValue()+1e-05)
		# calculate current timestep from maxvel, clamp range
		s.timestep = max( min( cflFac/maxvel, dtMax) , dtMin ) 
		if (t+s.timestep) > frameTime:
			# add epsilon to prevent roundoff errors...
			s.timestep = ( frameTime - t ) + 1e-04
		elif (t+s.timestep + dtMin) > frameTime or (t+(s.timestep*1.25)) > frameTime:
			# avoid tiny timesteps and strongly varying ones, do 2 medium size ones if necessary...
			s.timestep = (frameTime-t)*0.5
			lockDt = True
	
	# progress info
	print "Frame %f current max velocity: %f , dt: %f, %f/%f lock:%s" % (frame, maxvel, s.timestep, t, frameTime, lockDt)
	
	# NT_DEBUG , sanity check
	if (s.timestep < (dtMin/2) ):
		print "Error - Invalid dt encountered! Shouldnt happen..."; exit(1);


#main loop
for frame in range(frames):
	
	t = 0.
	lockDt = False
	while t<frameTime: 
		adaptDt(s)
		t += s.timestep
		
		# FLIP 
		pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )

		# make sure we have velocities throught liquid region
		mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
		extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )  # note, tmpVec3 could be free'd now...
		markFluidCells( parts=pp, flags=flags )

		# create approximate surface level set, resample particles
		gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
		unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) 
		phi.reinitMarching(flags=flags, velTransport=vel, correctOuterLayer=False ) # optionally, beautify levelset, needed for adjustNumber

		# set source grids for resampling, used in adjustNumber!
		pVel.setSource( vel, isMAC=True )
		pTest.setSource( tstGrid );
		adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

		# forces & pressure solve
		addGravity(flags=flags, vel=vel, gravity=(0,-0.003,0))
		setWallBcs(flags=flags, vel=vel)	
		solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
		setWallBcs(flags=flags, vel=vel)

		# make sure we have proper velocities
		extrapolateMACSimple( flags=flags, vel=vel , distance=(int(maxvel)+2+10), phiObs=phiObs ) # with knUnprojectNormalComp 
		
		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

		if 0 and (dim==3):
			phi.createMesh(mesh)
		
		#s.printMemInfo()
		#timings.display()
		s.step()

	# save particle data
	#pp.save( 'flipParts_%04d.uni' % frame );

	if 0 and (GUI):
		gui.screenshot( 'addt04_%04d.png' % frame );


