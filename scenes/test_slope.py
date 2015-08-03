from manta import *

new_BC     = True
dim        = 3
res        = 32
gs         = vec3(2*res,res,res)
if (dim==2): gs.z = 1
s          = FluidSolver(name='main', gridSize = gs, dim=dim)
s.timestep = 1
timings    = Timings()

flags     = s.create(FlagGrid)
vel       = s.create(MACGrid)
pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)
density   = s.create(RealGrid)
phiWalls  = s.create(RealGrid)
mesh      = s.create(Mesh)

flags.initDomain(inflow="xXY", phiWalls=phiWalls)
flags.fillGrid()

slope  = s.create(Slope, anglexy=-15.,angleyz=15.,origin=0.5*res, gs=gs)
phiObs = slope.computeLevelset()
#slope.generateMesh(mesh)

CombineLevelsets(phiObs,phiWalls)
updateFractions(flags=flags, phiObs=phiObs, fractions=fractions)

velInflow = vec3(.966,-.259,0)
initVelObs(phiObs=phiObs, vel=vel, val=velInflow)

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
		extrapolateMACSimple( flags=flags, vel=vel, distance=0 );
	else:
		setWallBcs(flags=flags, vel=vel)

	setInflowBcs(vel=vel,dir='xXY',value=velInflow)

	if(new_BC):
		solvePressure( flags=flags, vel=vel, pressure=pressure, fractions=fractions)
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
		extrapolateMACSimple( flags=flags, vel=vel, distance=0 );

	else:
		solvePressure( flags=flags, vel=vel, pressure=pressure )
		setWallBcs(flags=flags, vel=vel)

	setInflowBcs(vel=vel,dir='xXY',value=velInflow)
	
	# if (dim==3):
	# 	phiObs.createMesh(mesh)

	timings.display()
	s.step()

	if(t==5):
		if(new_BC):
			gui.screenshot( 'new_shot_%04d.png' % t);
		else:
			gui.screenshot( 'old_shot_%04d.png' % t);
