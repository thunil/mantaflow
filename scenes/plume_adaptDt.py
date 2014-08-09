#
# Simple example scene (hello world)
# Simulation of a buoyant smoke density plume

from manta import *


# how many frames to calculate 
frames    = 200

# maximal velocity per cell, adaptDt variables
cflFac = 3.0

# solver params
dim = 2
res = 64
#res = 34
gs = vec3(res,1.5*res,res)
gs = vec3(res,res,res) # NT_DEBUG
if (dim==2):
	gs.z=1
s = FluidSolver(name='main', gridSize = gs, dim=dim)
s.cfl         = cflFac
s.frameLength = 1.2          # length of one frame (in "world time")
s.timestepMin = 0.1
s.timestepMax = 2.0
s.timestep    = (s.timestepMax+s.timestepMin)*0.5
timings = Timings()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

flags.initDomain()
flags.fillGrid()

if (GUI):
	gui = Gui()
	gui.show( dim==2 )

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))


#main loop
while s.frame < frames:
	
	maxvel = vel.getMaxValue()
	s.adaptTimestep(maxvel)
	
	if s.timeTotal<50.:
		densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=1.0)
	
	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-6e-3,0), flags=flags)
	
	solvePressure( flags=flags, vel=vel, pressure=pressure )
	setWallBcs(flags=flags, vel=vel)
	#density.save('den%04d.uni' % t)
	
	#timings.display()
	s.step()

	if 0 and (GUI) and (lastFrame!=frame) and (s.frame%1==0):
		gui.screenshot( 'addt05_%04d.png' % s.frame );
	lastFrame = s.frame

