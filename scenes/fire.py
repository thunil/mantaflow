#
# Simulation of a flame with smoke (and with adaptive time-stepping)
#
from manta import *

# solver params
dim = 3
res = 32
gs = vec3(res, res, res)
if dim==2:
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# set time step range
s.frameLength = 1.2   # length of one frame (in "world time")
s.timestepMin = 0.2   # time step range
s.timestepMax = 2.0
s.cfl         = 3.0   # maximal velocity per cell
s.timestep    = (s.timestepMax+s.timestepMin)*0.5
timings = Timings()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
react = s.create(RealGrid)
fuel = s.create(RealGrid)
heat = s.create(RealGrid)
red = s.create(RealGrid)
green = s.create(RealGrid)
blue = s.create(RealGrid)
flame = s.create(RealGrid)
pressure = s.create(RealGrid)

# toggle open domain
doOpen = True

# how many frames to calculate 
frames = 250

# noise field
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

# needs positive gravity because of addHeatBuoyancy2()
gravity = vec3(0,0.0981,0)

# initialize domain with boundary
bWidth=1
flags.initDomain( boundaryWidth=bWidth )
flags.fillGrid()
if doOpen:
	setOpenBound( flags, bWidth,'yY',FlagOutflow|FlagEmpty )

if (GUI):
	gui = Gui()
	gui.show(True)
	gui.pause()

# source: cube in center of domain (x, y), standing on bottom of the domain
boxSize = vec3(res/8, res/8, res/8)
boxCenter = gs*vec3(0.5, 0.1, 0.5)
sourceBox = s.create( Box, center=boxCenter, size=boxSize )

# buoyancy parameters
smokeDensity = -0.001 # alpha
smokeTempDiff = 0.1 # beta

# main loop
while s.frame < frames:
	maxvel = vel.getMaxValue()
	s.adaptTimestep( maxvel )
	
	if s.timeTotal<200:
		densityInflow( flags=flags, density=density, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=heat, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=fuel, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=react, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=red, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=green, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=blue, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )

	processBurn( fuel=fuel, density=density, react=react, red=red, green=green, blue=blue, heat=heat )

	advectSemiLagrange( flags=flags, vel=vel, grid=density, order=2 )
	advectSemiLagrange( flags=flags, vel=vel, grid=heat, order=2 )
	advectSemiLagrange( flags=flags, vel=vel, grid=fuel, order=2 )
	advectSemiLagrange( flags=flags, vel=vel, grid=react, order=2 )
	advectSemiLagrange( flags=flags, vel=vel, grid=red, order=2 )
	advectSemiLagrange( flags=flags, vel=vel, grid=green, order=2 )
	advectSemiLagrange( flags=flags, vel=vel, grid=blue, order=2 )
	advectSemiLagrange( flags=flags, vel=vel, grid=vel, order=2, openBounds=doOpen, boundaryWidth=bWidth )

	if doOpen:
		resetOutflow( flags=flags, real=density )

	vorticityConfinement( vel=vel, flags=flags, strength=0.1 )

	addBuoyancy2( flags=flags, grid=density, vel=vel, gravity=gravity, coefficient=smokeDensity )
	addBuoyancy2( flags=flags, grid=heat, vel=vel, gravity=gravity, coefficient=smokeTempDiff )

	setWallBcs( flags=flags, vel=vel )
	solvePressure( flags=flags, vel=vel, pressure=pressure )

	updateFlame( react=react, flame=flame )

	timings.display()
	s.step()

