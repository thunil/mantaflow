#
# Obstacle test case
# 

import sys
from manta import *
from helperInclude import *

dim = 2
# solver params
res = 42
gs = vec3(res,1.5*res,res)
if (dim==2): gs.z = 1  # 2D
s = FluidSolver(name='main', gridSize = gs, dim=dim)
s.timestep = 1.0

flags     = s.create(FlagGrid)
vel       = s.create(MACGrid)
density   = s.create(RealGrid)
pressure  = s.create(RealGrid)
fractions = s.create(MACGrid)

flags.initDomain(0)
flags.fillGrid()

sphere  = s.create(Sphere, center=gs*vec3(-0.2,0.75,0.5), radius=res*0.75)
sphere.applyToGrid(grid=flags, value=FlagObstacle)
phi = sphere.computeLevelset()
updateFractions(flags=flags, phi=phi, fractions=fractions)

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.01, 0))

if 0 and (GUI):
    gui = Gui()
    gui.show(); gui.pause()

#main loop
for t in range(140):

	source.applyToGrid(grid=density, value=1)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2)
	
	addBuoyancy(flags=flags, density=density, vel=vel, gravity=vec3(0,-3e-3,0))
	
	#need to call update fractions before this for better obstacle BC
	if 1:
		# with volume fractions
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phi=phi)
		solvePressure(flags=flags, vel=vel, pressure=pressure, fractions=fractions)
		setWallBcs(flags=flags, vel=vel, fractions=fractions, phi=phi)
	else:
		# old boundaries
		setWallBcs(flags=flags, vel=vel)
		solvePressure(flags=flags, vel=vel, pressure=pressure)
		setWallBcs(flags=flags, vel=vel)

	s.step()

# check final state
doTestGrid( sys.argv[0],"dens" , s, density , threshold=0.00001 , thresholdStrict=1e-08 )
doTestGrid( sys.argv[0],"vel"  , s, vel     , threshold=0.00001 , thresholdStrict=1e-08 )
