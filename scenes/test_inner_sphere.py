from manta import *

new_BC     = True
dim        = 2
res        = 16
gs         = vec3(res,res,res)
if (dim==2): gs.z = 1
s          = FluidSolver(name='main', gridSize = gs, dim=dim)
s.timestep = 1
timings    = Timings()

flags     = s.create(FlagGrid)
vel       = s.create(MACGrid)
pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)
density   = s.create(RealGrid)

flags.initDomain()
flags.fillGrid()

center = gs*vec3(0.5,0.5,0.5)
radius = res*0.4
sphere = s.create(Sphere, center=center, radius=radius)
phiObs = sphere.computeLevelset()
phiObs.multConst(-1)

initVortexVelocity(phiObs=phiObs, vel=vel, center=center, radius=radius)

updateFractions(flags=flags, phiObs=phiObs, fractions=fractions)

if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

#main loop
for t in range(25000):

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, orderSpace=1)  
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=1.0)

	if(new_BC):
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
		extrapolateMACSimple( flags=flags, vel=vel, distance=1 );
	else:
		setWallBcs(flags=flags, vel=vel)

	if(new_BC):
		solvePressure( flags=flags, vel=vel, pressure=pressure, fractions=fractions)
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
		extrapolateMACSimple( flags=flags, vel=vel, distance=1 );

	else:
		solvePressure( flags=flags, vel=vel, pressure=pressure )
		setWallBcs(flags=flags, vel=vel)
	
	timings.display()
	s.step()

	if(t==5):
		if(new_BC):
			gui.screenshot( 'new_shot_%04d.png' % t);
		else:
			gui.screenshot( 'old_shot_%04d.png' % t);
