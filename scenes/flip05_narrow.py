#
# Flip scene with particle seeding only in a narrow band around the surface
# 
from manta import *

# how many frames to calculate 
frames    = 200

# no. of cells around surface
narrowBand = 5;

# solver params
dim = 3
res = 64
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# adaptive time stepping
s.frameLength = 0.6   # length of one frame (in "world time")
s.timestepMin = 0.2   # time step range
s.timestepMax = 2.0
s.cfl         = 2.5   # maximal velocity per cell
s.timestep    = (s.timestepMax+s.timestepMin)*0.5

minParticles = pow(2,dim)
timings = Timings()

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)
phiParts = s.create(LevelsetGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
velParts = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)

pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup, 0=breaking dam, 1=drop into pool
# geometry in world units (to be converted to grid space upon init)
setup = 1
flags.initDomain(boundaryWidth=0)
fluidVel = 0
fluidSetVel = 0

if setup==0:
	# breaking dam
	fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
	phi = fluidbox.computeLevelset()
	#fluidbox2 = s.create(Box, p0=gs*vec3(0.2,0.7,0.3), p1=gs*vec3(0.5,0.8,0.6)) # breaking dam
	#phi.join( fluidbox2.computeLevelset() )

	# obstacle init needs to go after updateFromLs
	obsBox = s.create(Box, p0=gs*vec3(0.7,0.0,0.5), p1=gs*vec3(0.8,1.0,0.8)) 
	obsBox.applyToGrid(grid=flags, value=(FlagObstacle) )
	#obsBox.applyToGrid(grid=flags, value=(FlagObstacle|FlagStick) )
	flags.updateFromLevelset(phi)

elif setup==1:
	# falling drop
	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.2,1.0)) # basin
	dropCenter = vec3(0.5,0.5,0.5)
	dropRadius = 0.15
	dropRadius = 0.075 #sm
	fluidSetVel= vec3(0,-0.03,0)
	fluidDrop  = s.create(Sphere, center=gs*dropCenter, radius=res*dropRadius)
	fluidVel   = s.create(Sphere, center=gs*dropCenter, radius=res*(dropRadius+0.05) )
	phi = fluidBasin.computeLevelset()
	phi.join( fluidDrop.computeLevelset() ) 
	flags.updateFromLevelset(phi)

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05 )
mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

if fluidVel!=0:
	# set initial velocity
	fluidVel.applyToGrid( grid=vel , value=gs*fluidSetVel )
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

if 1 and (GUI):
	gui = Gui()
	gui.show( dim==2 )
	gui.pause()
	  
	# show all particles shaded by velocity
	if(dim==2):
		gui.nextPdata()
		gui.nextPartDisplay()
		gui.nextPartDisplay()



#main loop
while s.frame < frames:
	
	maxVel = vel.getMaxValue()
	s.adaptTimestep( maxVel )
	
	# velocities are extrapolated at the end of each step
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

	# make sure we have velocities throught liquid region
	if narrowBand>0:
		mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 );
		combineGridVel(vel=velParts, weight=tmpVec3 , combineVel=vel );
	else:
		mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 );
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )  # note, tmpVec3 could be free'd/reused now...
	#markFluidCells( parts=pp, flags=flags )

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts , radiusFactor ) 

	if narrowBand>0:
		phi.addConst(1.); # shrink slightly
		phi.join( phiParts );
		extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True)
	else:
		phi.copyFrom( phiParts );
		extrapolateLsSimple(phi=phi, distance=4, inside=True)

	extrapolateLsSimple(phi=phi, distance=3) 
	#phi.reinitMarching(flags=flags, velTransport=vel, correctOuterLayer=False ) # optionally, beautify levelset, needed for adjustNumber
	flags.updateFromLevelset(phi)
	extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)) ) 

	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=(0,-0.003,0))
	setWallBcs(flags=flags, vel=vel)	
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
	setWallBcs(flags=flags, vel=vel)

	# make sure we have proper velocities
	extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)) )
	
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )

	# set source grids for resampling, used in adjustNumber!
	pVel.setSource( vel, isMAC=True )
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor , narrowBand=narrowBand ) 

	if 1 and (dim==3):
		phi.createMesh(mesh)
	
	#timings.display()
	#s.printMemInfo()
	s.step()

	# optionally particle data , or screenshot 
	if 0 and (GUI):
		nbid = 0;
		if narrowBand>0:
			nbid = int(narrowBand);
		gui.screenshot( 'flip05nb%02d_%04d.png' % (nbid,s.frame) );
	#pp.save( 'flipParts_%04d.uni' % s.frame );


