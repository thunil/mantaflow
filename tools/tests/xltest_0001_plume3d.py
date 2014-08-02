#
# Simple example scene (hello world)
# Simulation of a buoyant smoke density plume

import sys
from manta import *
from helperInclude import *

# solver params
res = 160
res = 100
#res = 50
gs = vec3(res,1.5*res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.5
testInterval = 1

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

# noise field
noise = s.create(NoiseField, loadFromFile=True )
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 0
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

flags.initDomain()
flags.fillGrid()

if 0 and (GUI):
	gui = Gui()
	gui.show()

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))
	
#main loop
for t in range(100):
	densityInflow(flags=flags, density=density, noise=noise, shape=source, scale=1, sigma=0.5)
		
	#source.applyToGrid(grid=vel, value=velInflow)
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)	
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
	
	setWallBcs(flags=flags, vel=vel)	
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-5e-2,0), flags=flags)
	
	#solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=1e-04)
	solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=2.5, cgAccuracy=1e-08)
	setWallBcs(flags=flags, vel=vel)

	# check state in intervals
	if(t%testInterval == (testInterval-1)):
		doTestGrid( sys.argv[0], ("dens_%04d"  % t), s, density , threshold=0.001 , thresholdStrict=1e-08 )
		doTestGrid( sys.argv[0], ("vel_%04d"   % t), s, vel     , threshold=0.001 , thresholdStrict=1e-08 )
		#doTestDataLoad( sys.argv[0], ("dens_%04d"  % t), s, density )
		#doTestDataLoad( sys.argv[0], ("vel_%04d"   % t), s, vel     )
	
	s.step()


