from manta import *

new_BC     = True
dim        = 2
res        = 64
gs         = vec3(2*res,res,res)
if (dim==2): gs.z = 1
s          = FluidSolver(name='main', gridSize = gs, dim=dim)
s.timestep = 1.
timings    = Timings()
bWidth     = 0	# number of ghostcells minus 1

flags     = s.create(FlagGrid)
vel       = s.create(MACGrid)
pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)
density   = s.create(RealGrid)
phiWalls  = s.create(RealGrid)

flags.initDomain(inflow="xX", phiWalls=phiWalls, boundaryWidth=bWidth)
flags.fillGrid()

sphere  = s.create(Sphere, center=gs*vec3(0.25,0.5,0.5), radius=res*0.2)
phiObs = sphere.computeLevelset()

CombineLevelsets(phiObs,phiWalls)
updateFractions(flags=flags, phiObs=phiObs, fractions=fractions)

velInflow = vec3(.5,0,0)
initVel(phiObs=phiObs, vel=vel, val=velInflow)

if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

#main loop
for t in range(25000):

	applyDensAtObstacle(phiObs=phiObs, dens=density)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2, orderSpace=1)  
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=1.0)

	if(new_BC):
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
		extrapolateMACSimple( flags=flags, vel=vel, distance=bWidth );
	else:
		setWallBcs(flags=flags, vel=vel)

	setInflowBcs(vel=vel,dir='xX',value=velInflow)

	if(new_BC):
		solvePressure( flags=flags, vel=vel, pressure=pressure, fractions=fractions)

		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
		extrapolateMACSimple( flags=flags, vel=vel, distance=bWidth );
	else:
		solvePressure( flags=flags, vel=vel, pressure=pressure )

		setWallBcs(flags=flags, vel=vel)

	setInflowBcs(vel=vel,dir='xX',value=velInflow)

	timings.display()
	s.step()

	if(t==5):
		if(new_BC):
			gui.screenshot( 'new_shot_%04d.png' % t);
		else:
			gui.screenshot( 'old_shot_%04d.png' % t);