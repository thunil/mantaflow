#
# Simple flip with level set and basic resampling
# 
from manta import *

# solver params
dim = 3
res = 64
res = 120
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5
minParticles = pow(2,dim)

# size of particles 
radiusFactor = 5.1

# prepare grids and particles
flags    = s.create(FlagGrid)
phi      = s.create(LevelsetGrid)

pp       = s.create(BasicParticleSystem) 
pPos     = pp.create(PdataVec3) 
pVel     = pp.create(PdataVec3) 
# test real value, not necessary for simulation
#pTest    = pp.create(PdataReal) 
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
for t in range(2500):
	# debug, FF
	t = t * 5
	# read input sim
	if (dim==2):
		pp.load( 'flipOut02d_%04d.uni' % t );
	if (dim==3):
		pp.load( 'flipOut05_3d_%04d.uni' % t );
	# adjust local coordinates to current grid
	pp.getPosPdata(pPos)
	# ...
	# NT_DBEUG , todo add basic ops, like grids
	pp.setPosPdata(pPos)
	testDiscardNth( pp , 25000 )
	
	# create surface
	phi.clear()
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	#unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) 
	averagedParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) 
	# NT_DEBUG , why does marching need flags?
	phi.reinitMarching(flags=flags, maxTime=int(2*radiusFactor) )


	if (dim==3):
		phi.createMesh(mesh)
	
	#s.printMemInfo()
	s.printTimings()
	s.step()

	#pp.save( 'flipOut01_%04d.uni' % t );
	if 0 and (GUI):
		gui.screenshot( 'flipt11n3d_%04d.png' % t );

