#
# Simple example for free-surface simulation 
#

from manta import *

# solver params
dim = 3
res = 64
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep  = 0.25

# scene file params
ghostFluid  = True
doOpen      = False
accuracy    = 5e-4
# using fast marching is more accurate, but currently causes asymmetries
useMarching = False
# level set advection 1st/2nd order, 2nd order also introduces slight asymmetries
lsOrder     = 1

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)
phiBackup = s.create(LevelsetGrid)

# scene setup
bndWidth=1
flags.initDomain(boundaryWidth=bndWidth)
basin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,0.2,1))
drop  = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.125)
phi = basin.computeLevelset()
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)

# optionally, enable open boundaries here and below...
if doOpen:
	setOpenBound(flags,bndWidth,'xXzZ',FlagOutflow|FlagEmpty) 
		
if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
	

#main loop
for t in range(2000):
	
	# update and advect levelset
	if useMarching:
		phi.reinitMarching(flags=flags, velTransport=vel) 
	else:
		extrapolateLsSimple(phi=phi, distance=5, inside=False)
		extrapolateLsSimple(phi=phi, distance=5, inside=True )
		extrapolateMACSimple( flags=flags, vel=vel, distance=5 )

	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=lsOrder) 
	phi.setBoundNeumann(bndWidth)
	if doOpen:
		resetOutflow(flags=flags,phi=phi) # open boundaries
	flags.updateFromLevelset(phi)
	
	# velocity self-advection
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2 )
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
	#gui.screenshot( 'freesurface_%04d.png' % t );



