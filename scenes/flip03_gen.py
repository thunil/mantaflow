#
# Simple flip surface mesh creation scene
# 
from manta import *

# solver params
dim = 3
res = 128 
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 1.0

# kernel radius for surface creation
radiusFactor = 2.5

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
	#gui.pause()
   

#main loop
for tOut in range(2500):
	# for testing, optionally skip input frames...
	tIn = tOut * 1
	# read input sim particles
	pp.load( 'flipParts_%04d.uni' % tIn );
	
	# create surface
	phi.clear()
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	#unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) # faster, but not as smooth
	averagedParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor , 1, 1 ) 

	if (dim==3):
		setBoundaries(phi, 0., boundaryWidth=1)
		phi.createMesh(mesh)
		# output file name so that blender can directly read it...
		mesh.save( 'fluidsurface_final_%04d.bobj.gz' % tOut );
	
	s.step()

