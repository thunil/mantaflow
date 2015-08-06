from manta import *

new_BC     = True
dim        = 2
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
phiWalls  = s.create(LevelsetGrid)
mesh      = s.create(Mesh)

flags.initDomain(inflow="xXY", phiWalls=phiWalls)

slope  = s.create(Slope, anglexy=-15.,angleyz=15.,origin=0.5*res, gs=gs)
phiObs = slope.computeLevelset()
#slope.generateMesh(mesh)

phiObs.join(phiWalls)
updateFractions(flags=flags, phiObs=phiObs, fractions=fractions)
setObstacleFlags(flags=flags, phiObs=phiObs, fractions=fractions)
flags.fillGrid()

velInflow = vec3(.966,-.259,0)
vel.setConst(velInflow)


if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

#main loop
for t in range(25000):

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2)

	if(new_BC):
		extrapolateMACSimple( flags=flags, vel=vel, distance=2 , intoObs=True);
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)
	else:
		setWallBcs(flags=flags, vel=vel)

	setInflowBcs(vel=vel,dir='xXY',value=velInflow)

	if(new_BC):
		solvePressure( flags=flags, vel=vel, pressure=pressure, fractions=fractions)

		extrapolateMACSimple( flags=flags, vel=vel, distance=5 , intoObs=True);
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phiObs=phiObs)

	else:
		solvePressure( flags=flags, vel=vel, pressure=pressure )
		setWallBcs(flags=flags, vel=vel)

	setInflowBcs(vel=vel,dir='xXY',value=velInflow)
	
	# if (dim==3):
	# 	phiObs.createMesh(mesh)

	timings.display()
	s.step()

	if 0 and (t==5):
		if(new_BC):
			gui.screenshot( 'new_shot_%04d.png' % t);
		else:
			gui.screenshot( 'old_shot_%04d.png' % t);
			


