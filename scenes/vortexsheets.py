#
# Example scene for mesh-based smoke simulations using Vortex Sheets
# ( http://graphics.ethz.ch/publications/papers/paperPfa12.php )
#
# This example requires CUDA to run

from manta import *
assert(CUDA), 'Requires CUDA. Please compile with -DCUDA=ON'

scale = 0.2

# solver params
res = 64
gs = vec3(res,res*1.5,res)
dx = 1/(1.5*res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5

gravity = vec3(0,-0.15,0)
velInflow = vec3(0,0.52,0)

# prepare grids
flags = s.create(FlagGrid)
pressure = s.create(RealGrid, show=False)
vel = s.create(MACGrid)
density = s.create(RealGrid)

# noise field
noise = s.create(NoiseField)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

mesh = s.create(VortexSheetMesh)
vp = s.create(VortexParticleSystem)

flags.initDomain()
flags.fillGrid()

if (GUI):
	gui = Gui()
	gui.show()

# inflow mesh
source = s.create(Cylinder, center=gs*vec3(0.5,0.13,0.5), radius=res*0.14, z=gs*vec3(0, 0.03, 0))
mesh.fromShape(source)
subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=2*scale)

# fix the nodes at the bottom of mesh
fixedRegion = Box( parent=s, center=gs*vec3(0.5,0.09,0.5), size=gs*vec3(0.4,0.03,0.4))
	
#main loop
for t in range(180):
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	# seed smoke within the source region and apply inflow velocity condition
	densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
	source.applyToGrid(grid=vel, value=velInflow)
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)
	markAsFixed(mesh=mesh, shape=fixedRegion)

	# advect mesh in large-scale grid
	mesh.calcCirculation()        
	mesh.advectInGrid(vel=vel, flags=flags, integrationMode=IntRK4)  
	mesh.calcVorticity()
	
	# vortex sheet integration
	#filterVorticityCuda(mesh=mesh, sigma=1.6)
	smoothVorticity(mesh=mesh, iter=50, alpha=1.0, sigma=0.6) # fast approximation for Gauss kernel
	meshApplyBuoyancyLocalCuda(mesh=mesh, scale=0.6*1e-2, cutoffCells=5, regularization=0.5)
	mesh.calcVorticity()
	
	# mesh cosmetics
	smoothMesh(mesh=mesh, strength=1e-4, steps=1)
	subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=2*scale, cutTubes=True)
	killSmallComponents(mesh=mesh, elements=20)
	
	# base solver
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
	setWallBcs(flags=flags, vel=vel)
	addBuoyancy(density=density, vel=vel, gravity=gravity*dx, flags=flags)
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setWallBcs(flags=flags, vel=vel)
	
	# vorticity source terms
	vorticitySource(mesh=mesh, gravity=gravity, scale=0.1, maxAmount=200)
	
	#mesh.save('d%04d.bobj.gz' % t)
	s.step()
	
