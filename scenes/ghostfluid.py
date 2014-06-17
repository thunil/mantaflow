#
# Smooth surface for ghost fluid testing
# 

from manta import *

# solver params
dim = 2
res = 64
res = 128
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.25
ghostFluid = True
accuracy = 5e-5

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=0)
basin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,0.2,1))
drop  = s.create(Sphere, center=gs*vec3(0.5,0.68,0.5), radius=res*0.15) 
#drop  = s.create(Sphere, center=gs*vec3(0.5,-0.28,0.5), radius=res*0.5) # low, gf test
phi = basin.computeLevelset()
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)
		
if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
	

#main loop
for t in range(2000):
	
	phi.reinitMarching(flags=flags, velTransport=vel) 
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
	flags.updateFromLevelset(phi)
	
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
	addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.025,0))
	
	# print current maximal velocity
	maxvel = vel.getMaxValue()
	print "Current max velocity %f " % maxvel
   
	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	if ghostFluid:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, useResNorm=False, \
			phi=phi ,  gfClamp=0.0001 ) 
	else:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, useResNorm=False)
	setWallBcs(flags=flags, vel=vel)
	
	if (dim==3):
		phi.createMesh(mesh)
	
	s.step()
	#gui.screenshot( 'out_%04d.png' % t );



