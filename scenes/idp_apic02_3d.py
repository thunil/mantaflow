#
# Very simple 3d dam break scene with apic
#
from manta import *

# solver params
dim = 3
particleNumber = 2
res = 64
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
	particleNumber = 3      # use more particles in 2d

s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 1.0

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
velOld = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3 = s.create(VecGrid)
pp = s.create(BasicParticleSystem)
# add velocity data to particles
pVel = pp.create(PdataVec3)
phiObs = s.create(LevelsetGrid, name='phiObs')
# apic part
apic_mass = s.create(MACGrid)
apic_pCx = pp.create(PdataVec3)
apic_pCy = pp.create(PdataVec3)
apic_pCz = pp.create(PdataVec3)

#position solver stuff
usePositionSolver = True
density = s.create(RealGrid)
Lambda = s.create(RealGrid)		#pressure of the position solver (p_2 in the paper)
deltaX = s.create(MACGrid)
flagsPos = s.create(FlagGrid)
pMass = pp.create(PdataReal)
mass = 1.0 / (particleNumber * particleNumber * particleNumber) 
if (dim==2):
	mass = 1.0 / (particleNumber * particleNumber) 

resampleParticles = False

if (resampleParticles):
	pindex = s.create(ParticleIndexSystem) 
	gpi = s.create(IntGrid)
	gCnt = s.create(IntGrid)

s.timestep = 1
s.frameLength = 10000000.0   # length of one frame (in "world time")
s.timestepMin = 0.01   # time step range
s.timestepMax = 1.0
s.cfl = 5.0   # maximal velocity per cell

# scene setup
flags.initDomain(boundaryWidth=1)
fluidbox = Box(parent=s, p0=gs*vec3(0,0,0.25), p1=gs*vec3(0.5,0.35,0.75)) # breaking dam
phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

sampleFlagsWithParticles(flags=flags, parts=pp, discretization=particleNumber, randomness=0.5)
copyFlagsToFlags(flags, flagsPos)
flags.initDomain(boundaryWidth=1, phiWalls=phiObs)

adaptiveTimeSteps=True

if (GUI):
	gui = Gui()
	gui.show()
	# show all particles shaded by velocity
	if(dim==3):
		gui.nextPdata()
		gui.nextPartDisplay()
		gui.nextPartDisplay()
	gui.pause()

#main loop
for t in range(2000):
	# use adaptive time stepping
	if (adaptiveTimeSteps):
		maxVel = vel.getMax()
		s.adaptTimestep( maxVel )

	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	# APIC
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=2, deleteInObstacle=False, stopInObstacle = False) 
		
	if (usePositionSolver):
		copyFlagsToFlags(flags, flagsPos)
		mapMassToGrid(flags=flagsPos, density=density, parts=pp, source=pMass, deltaX=deltaX, phiObs=phiObs, dt=s.timestep, particleMass=mass, noDensityClamping =  resampleParticles)			
		
		# This is not implemented in the most efficient way (see comment in the C++ code) and the current implementation has some computational overhead. 
		# The degenerate configurations rarely occur so that it is ok to switch this off for most scenarios.
		if (resampleParticles):
			gridParticleIndex(parts=pp, indexSys=pindex, flags=flags, index=gpi, counter=gCnt)
			apicMapPartsToMAC(flags=flags, vel=vel, parts=pp, partVel=pVel, cpx=apic_pCx, cpy=apic_pCy, cpz=apic_pCz, mass=apic_mass)
			resampeOverfullCells(vel=vel, density=density, index=gpi, indexSys=pindex, part=pp, pVel=pVel, dt=s.timestep)
	
		# position solver
		solvePressureSystem(rhs=density, vel=vel, pressure=Lambda, flags=flagsPos, cgAccuracy = 1e-3)
		computeDeltaX(deltaX=deltaX, Lambda=Lambda, flags=flagsPos)
		mapMACToPartPositions(flags=flagsPos, deltaX=deltaX, parts=pp, dt=s.timestep)
				
	apicMapPartsToMAC(flags=flags, vel=vel, parts=pp, partVel=pVel, cpx=apic_pCx, cpy=apic_pCy, cpz=apic_pCz, mass=apic_mass)
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
	markFluidCells( parts=pp, flags=flags)
	
	addGravityNoScale(flags=flags, vel=vel, gravity=(0,-0.01,0))

	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	solvePressure(flags=flags, vel=vel, pressure=pressure, cgAccuracy = 1e-3)
	setWallBcs(flags=flags, vel=vel)
	
	# we dont have any levelset, ie no extrapolation, so make sure the velocities are valid. For small time steps distance can be reduced
	extrapolateMACSimple( flags=flags, vel=vel, distance=5 )
	
	# APIC velocity update
	apicMapMACGridToParts(partVel=pVel, cpx=apic_pCx, cpy=apic_pCy, cpz=apic_pCz, parts=pp, vel=vel, flags=flags)
			
	s.step()
