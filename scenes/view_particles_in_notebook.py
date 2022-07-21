#
# NUMPY file format test
#
from manta import *

# solver params
dim = 3
# dim = 2
particleNumber = 2
res = 64
gs = vec3(res,res,res)

if (dim==2):
	gs.z=1
	particleNumber = 3      # use more particles in 2d

s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5

# prepare grids and particles
flags    = s.create(FlagGrid)
flags2   = s.create(FlagGrid)
vel      = s.create(MACGrid)
vel2      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
pressure2 = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
pp       = s.create(BasicParticleSystem)
# add velocity data to particles
pVel     = pp.create(PdataVec3)
# apic part
mass  = s.create(MACGrid)
pCx   = pp.create(PdataVec3)
pCy   = pp.create(PdataVec3)
pCz   = pp.create(PdataVec3)

# scene setup
flags.initDomain(boundaryWidth=0)
# enable one of the following
# fluidbox = Box(parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
fluidbox = Box(parent=s, p0=gs*vec3(0.4,0.72,0.4), p1=gs*vec3(0.6,0.92,0.6)) # centered falling block
phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles(flags=flags, parts=pp, discretization=particleNumber, randomness=0.2)

#main loop
for t in range(200):
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	# Affine Particle in Cell (APIC) method
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False)
	apicMapPartsToMAC(flags=flags, vel=vel, parts=pp, partVel=pVel, cpx=pCx, cpy=pCy, cpz=pCz, mass=mass)
	extrapolateMACFromWeight(vel=vel , distance=2, weight=tmpVec3)
	markFluidCells(parts=pp, flags=flags)

	addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))

	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setWallBcs(flags=flags, vel=vel)

	# we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
	extrapolateMACSimple(flags=flags, vel=vel)

	# APIC velocity update
	apicMapMACGridToParts(partVel=pVel, cpx=pCx, cpy=pCy, cpz=pCz, parts=pp, vel=vel, flags=flags)

	s.step()

	file_name = './frames/particles-frame-%d.nptxt' % s.frame
	pp.save(file_name)
	

	


