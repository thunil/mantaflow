#
# Slightly more complicated flip setup
# uses levelset surface, resampling, and ghost fluid BCs
# 
from manta import *

# solver params
dim = 3
res = 44
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
pressure = s.create(RealGrid)
phi      = s.create(LevelsetGrid)

flip     = s.create(FlipSystem) 
mesh     = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=0)
# enable one of the following
#fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
#fluidbox = s.create(Box, p0=gs*vec3(0.4,0.4,0.4), p1=gs*vec3(0.6,0.8,0.6)) # centered falling block
fluidbox = s.create(Box, p0=gs*vec3(0.3,0.4,0.3), p1=gs*vec3(0.7,0.8,0.7)) # centered falling block
phi = fluidbox.computeLevelset()
flags.updateFromLevelset(phi)

flip.initialize( flags=flags, discretization=2, randomness=0.2 )
minParticles = pow(2,dim)

if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

#main loop
for t in range(2500):
	
	# FLIP 
	flip.advectInGrid(flaggrid=flags, vel=vel, integrationMode=IntRK4)
	flip.velocitiesToGrid(vel=vel, flags=flags)
	flip.markFluidCells(flags=flags) # dont use levelset for flags!

	# create simple surface, resample particles
	unionParticleLevelset( flip, phi )
	phi.reinitMarching(flags=flags, maxTime=2 )
	flip.adjustNumber( vel=vel, flags=flags, minParticles=minParticles, maxParticles=2*minParticles, phi=phi ) 
	
	addGravity(flags=flags, vel=vel, gravity=(0.0,-0.002,0))
	
	# pressure solve
	setWallBcs(flags=flags, vel=vel) 
	# with ghost fluid
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
	setWallBcs(flags=flags, vel=vel)
	extrapolateMACSimple( flags=flags , vel=vel )
	
	# FLIP velocity update
	flip.velocitiesFromGrid(vel=vel, flags=flags, flipRatio=0.97)

	if (dim==3):
		phi.createMesh(mesh)
		#mesh.save('phi%04d.bobj.gz' % t)
	
	#s.printTimings()
	#gui.screenshot( 'flipt_%04d.png' % t );
	s.step()

