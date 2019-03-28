#
# NUMPY file format test
#
from manta import *

# solver params
dim = 2
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
fluidbox = Box(parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
#fluidbox = Box(parent=s, p0=gs*vec3(0.4,0.72,0.4), p1=gs*vec3(0.6,0.92,0.6)) # centered falling block
phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles(flags=flags, parts=pp, discretization=particleNumber, randomness=0.2)


if 0 and (GUI):
	gui = Gui()
	gui.show()
	#gui.pause()

#main loop
for t in range(2500):
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	# APIC
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

	#gui.screenshot('flipt_%04d.png' % t);
	s.step()

	pressure.save("pressure_TEST.npz".format(t))
	pressure2.load("pressure_TEST.npz".format(t))
	mantaMsg('\nTesting Real Grid numpy write/read')
	mantaMsg('Min/Max Orig: %f %f' % (pressure.getMin(), pressure.getMax()))
	mantaMsg('Min/Max New: %f %f' % (pressure2.getMin(), pressure2.getMax()))
	pressure2.sub(pressure)
	mantaMsg('Min/Max New (sub old): %f %f' % (pressure2.getMin(), pressure2.getMax()))

	vel.save("velocity_TEST.npz".format(t))
	vel2.load("velocity_TEST.npz".format(t))
	mantaMsg('\nTesting MAC Grid numpy write/read')
	mantaMsg('Min/Max Orig: %f %f' % (vel.getMin(), vel.getMax()))
	mantaMsg('Min/Max New: %f %f' % (vel2.getMin(), vel2.getMax()))
	vel2.sub(vel)
	mantaMsg('Min/Max New (sub old): %f %f' % (vel2.getMin(), vel2.getMax()))

	flags.save("flags_TEST.npz".format(t))
	flags2.load("flags_TEST.npz".format(t))
	mantaMsg('\nTesting Int Grid numpy write/read')
	mantaMsg('Min/Max Orig: %f %f' % (flags.getMin(), flags.getMax()))
	mantaMsg('Min/Max New: %f %f' % (flags2.getMin(), flags2.getMax()))
	flags2.sub(flags)
	mantaMsg('Min/Max New (sub old): %f %f' % (flags2.getMin(), flags2.getMax()))

# externally test load in python with:
#>>> import numpy as np
#>>> v = np.load("velocity_TEST.npz")
#>>> print(format(v["grid"].shape))

