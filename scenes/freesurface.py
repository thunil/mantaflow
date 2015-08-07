#
# Simple example for free-surface simulation 
# (and optionally second order free surface boundaries, or open/outflow boundaries)
#

from manta import *

# solver params
dim = 2
res = 64
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.25
accuracy = 5e-5

ghostFluid = True
doOpen     = False

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

# scene setup
bWidth=1
flags.initDomain(boundaryWidth=bWidth)
basin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,0.2,1))
drop  = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.15)
phi = basin.computeLevelset()
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)

if doOpen:
	setOpenBound(flags,bWidth,'xXzZ',FlagOutflow|FlagEmpty) 
		
if (GUI):
	gui = Gui()
	gui.show()
	#gui.pause()
	

#main loop
for t in range(2000):
	
	# update and advect levelset
	phi.reinitMarching(flags=flags, velTransport=vel) 
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
	resetPhiInObs(flags, phi)

	if doOpen:
		resetOutflow(flags=flags,phi=phi)
	flags.updateFromLevelset(phi)
	
	# velocity self-advection
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=1, openBounds=doOpen, depth=bWidth+1)
	addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.025,0))
	
	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	if ghostFluid:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, phi=phi )
	else:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy)
	setWallBcs(flags=flags, vel=vel)
	
	if (dim==3):
		phi.createMesh(mesh)
		#mesh.save('phi%04d.bobj.gz' % t)
	
	s.step()
	#gui.pause()
	#gui.screenshot( 'screenOn_%04d.png' % t );



