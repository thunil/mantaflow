#
# More complex flip setup, breaking dam with resampling and
# additional particle data fields
# 
import sys
from manta import *
from helperInclude import *

# solver params
dim = 3
res = 32
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1

basename1 = "t01.py"
#basename2 = basename1 # optionally, make different
basename2 = "t02.py"

s = Solver(name='main', gridSize = gs, dim=dim)
flags    = s.create(FlagGrid)
real1    = s.create(RealGrid)
real2    = s.create(RealGrid)
#ls1      = s.create(LevelsetGrid)
#ls2      = s.create(LevelsetGrid)
#mac1     = s.create(MACGrid)
#mac2     = s.create(MACGrid)
realErr   = s.create(RealGrid)

parts1    = s.create(BasicParticleSystem) 
pDens1    = parts1.create(PdataReal) 
pTest1    = parts1.create(PdataReal) 

parts2    = s.create(BasicParticleSystem) 
pDens2    = parts2.create(PdataReal) 

flags.initDomain(boundaryWidth=0)

if 1 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
    
# try to load uni file if it exists
def tryToLoad( grid, basename, suffix ):
	rfile = referenceFilename( basename, suffix ) 
	if(os.path.isfile(rfile)):
		grid.load(rfile)
	else:
		grid.clear()
	return 1


#main loop
for t in range(150):

	tryToLoad( real1, basename1, ("dens_%04d"  % t) )
	tryToLoad( real2, basename2, ("dens_%04d"  % t) )
	realErr.sub(real1,real2);
	realErrMax = gridMaxDiff(real1, real2)

	# load particles
	tryToLoad( parts1 , basename1, ("parts_%04d"  % t) )
	tryToLoad( parts2 , basename2, ("parts_%04d"  % t) )

	#tryToLoad( pDens1 , basename1, ("pdens_%04d"  % t) )
	
	realErr.print(zSlice=15) 
	print("Max difference in step " +str(t) + " = "+ str(realErrMax) )

	s.step()
    
