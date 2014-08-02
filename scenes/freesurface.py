#
# Simple example for free-surface simulation with MacCormack advection
#

from manta import *

# solver params
dim = 3
res = 64
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
drop  = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.15)
phi = basin.computeLevelset()
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)
		
if (GUI):
	gui = Gui()
	gui.show()
	#gui.pause()
	

#main loop
for t in range(2000):
	
	# update and advect levelset
	phi.reinitMarching(flags=flags, velTransport=vel) #, ignoreWalls=False)
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
	flags.updateFromLevelset(phi)
	
	# velocity self-advection
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
	addGravity(flags=flags, vel=vel, gravity=vec3(0,-0.025,0))
	
	# pressure solve
	setWallBcs(flags=flags, vel=vel)
	if ghostFluid:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, phi=phi )
	else:
		solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy)
	setWallBcs(flags=flags, vel=vel)
	
	# note: these meshes are created by fast marching only, should smooth
	#       geometry and normals before rendering (only in 3D for now)
	if (dim==3):
		phi.createMesh(mesh)
		#mesh.save('phi%04d.bobj.gz' % t)
	
	s.step()
	#gui.pause()
	#gui.screenshot( 'screenOn_%04d.png' % t );



