#
# Simple flip with level set and basic resampling
# 
from manta import *

# solver params
new_BC = False
dim    = 2
res    = 64
gs     = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.8
minParticles = pow(2,dim)

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags     = s.create(FlagGrid)
phi       = s.create(LevelsetGrid)

vel       = s.create(MACGrid)
velOld    = s.create(MACGrid)
pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)
tmpVec3   = s.create(VecGrid)
tstGrid   = s.create(RealGrid)
phiWalls  = s.create(LevelsetGrid)

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
# test real value, not necessary for simulation
pTest    = pp.create(PdataReal) 
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup, 0=breaking dam, 1=drop into pool
setup = 0
bWidth=1
flags.initDomain(boundaryWidth=bWidth, phiWalls=phiWalls )
fluidVel = 0
fluidSetVel = 0

if setup==0:
	# breaking dam
	fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
	#fluidbox = s.create(Box, p0=gs*vec3(0.4,0.72,0.4), p1=gs*vec3(0.6,0.92,0.6)) # centered falling block
	phi = fluidbox.computeLevelset()
	sphere = s.create(Sphere, center=gs*vec3(0.66,0.3,0.5), radius=res*0.2)
	phiObs = sphere.computeLevelset()
	phiObs.join(phiWalls)
elif setup==1:
	# falling drop
	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.1,1.0)) # basin
	dropCenter = vec3(0.5,0.3,0.5)
	dropRadius = 0.1
	fluidDrop  = s.create(Sphere, center=gs*dropCenter, radius=res*dropRadius)
	fluidVel   = s.create(Sphere, center=gs*dropCenter, radius=res*(dropRadius+0.05) )
	fluidSetVel= vec3(0,-1,0)
	phi = fluidBasin.computeLevelset()
	phi.join( fluidDrop.computeLevelset() )

flags.updateFromLevelset(phi)
#setOpenBound(flags,bWidth,'xX',FlagOutflow|FlagEmpty)
sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05 )

if fluidVel!=0:
	# set initial velocity
	fluidVel.applyToGrid( grid=vel , value=fluidSetVel )
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

# testing the real channel while resampling - original particles
# will have a value of 0.1, new particle will get a value from the tstGrid
testInitGridWithPos(tstGrid)
pTest.setConst( 0.1 )

updateFractions(flags=flags, phiObs=phiObs, fractions=fractions)

if 1 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
   

#main loop
for t in range(2500):
	
	# FLIP 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=True)

	# make sure we have velocities throught liquid region
	mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )  # note, tmpVec3 could be free'd now...
	markFluidCells( parts=pp, flags=flags )

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) 
	#resetOutflow(flags=flags,parts=pp,index=gpi,indexSys=pindex) 
	# extend levelset somewhat, needed by particle resampling in adjustNumber
	extrapolateLsSimple(phi=phi, distance=4, inside=True); 

	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=(0,-0.001,0))
	if(new_BC):
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)	
		solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi, fractions=fractions)
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
	else:
		setWallBcs(flags=flags, vel=vel)	
		solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
		setWallBcs(flags=flags, vel=vel)
	# set source grids for resampling, used in adjustNumber!
	pVel.setSource( vel, isMAC=True )
	pTest.setSource( tstGrid );
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

	# make sure we have proper velocities
	extrapolateMACSimple( flags=flags, vel=vel )
	
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

	if (dim==3):
		phi.createMesh(mesh)
	
	#s.printMemInfo()
	s.step()

	# generate data for flip03_gen.py surface generation scene
	#pp.save( 'flipParts_%04d.uni' % t );

	if 0 and (GUI):
		gui.screenshot( 'flip02_%04d.png' % t );


