#
# Simple flip surface mesh creation scene
# 
from manta import *
import os
setDebugLevel(10) # full output

# input file 
partfile = 'flipParts_%04d.uni' 

# output file name so that blender can directly read it...
meshfile = 'fluidsurface_final_%04d.bobj.gz' 

# solver params
dim = 3
res = 128 
res = 40 
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 1.0

# kernel radius for surface creation
radiusFactor = 2.5
# triangle scale relative to cell size
scale = 0.5

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)
pp       = s.create(BasicParticleSystem) 
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# scene setup
flags.initDomain(boundaryWidth=0)
	
if 1 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
   

#main loop
for tOut in range(2500):
	# for testing, optionally skip input frames...
	tIn = tOut * 1

	meshfileCurr = meshfile % tOut 

	# already exists?
	if (os.path.isfile( meshfileCurr )):
		mesh.load( meshfileCurr )

	else:
		# generate mesh; first read input sim particles
		pp.load( partfile % tIn );
		
		# create surface
		phi.clear()
		gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
		#unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) # faster, but not as smooth
		averagedParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor , 1, 1 ) 

		if (dim==3):
			setBoundaries(phi, 0., boundaryWidth=1)
			phi.createMesh(mesh)

			# too slow right now!
			#subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=3*scale, cutTubes=False) 
			# beautify mesh
			#for iters in range(0):
				#smoothMesh(mesh=mesh, strength=1e-3, steps=10) 
				#subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=3*scale, cutTubes=True)

			# write output file:
			#mesh.save( meshfileCurr )
		
	s.step()

