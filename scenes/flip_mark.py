# This is a test file for the upcoming full FLIP with surface reconstruction kernels
# It will not work with the current mantaflow version yet.
#

from manta import *

assert(FULLFLIP)

# solver params
meshing = True
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
flip = s.create(FlipSystem)
mesh = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=0)
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,0.3,1))
phi = fluidbox.computeLevelset()
flags.updateFromLevelset(phi)
flip.initialize(flags=flags,discretization=2)

# moving obstacle
obsFrom = gs*vec3(0.3,0.3,0.5)
obsTo = gs*vec3(0.7,0.3,0.5)
box = s.create(Box, center=obsFrom, size=gs*vec3(0.1,0.2,0.3))
mov = s.create(MovingObstacle)
mov.add(box)

if (GUI):
	gui = Gui()
	gui.show()

#main loop
for t in range(400):
	
	# FLIP advect and writeback
	flip.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4)
	mov.projectOutside(flags=flags,flip=flip)
	flip.velocitiesToGrid(vel=vel, flags=flags)
	flip.markFluidCells(flags=flags)
	#advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
	
	addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))
	
	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	setLiquidBcs(flags=flags, vel=vel)
	mov.moveLinear(t=t,t0=0,t1=140,p0=obsFrom,p1=obsTo,flags=flags,vel=vel,smooth=True)
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setLiquidBcs(flags=flags, vel=vel)    
	setWallBcs(flags=flags, vel=vel)
	mov.moveLinear(t=t,t0=0,t1=140,p0=obsFrom,p1=obsTo,flags=flags,vel=vel,smooth=True)
	
	# recompute levelset, extrapolate velocities
	phi.initFromFlags(flags=flags)
	phi.reinitMarching(flags=flags, velTransport=vel)
	
	# FLIP load
	flip.velocitiesFromGrid(vel=vel, flags=flags, flipRatio=0.96)
	
	s.step()
