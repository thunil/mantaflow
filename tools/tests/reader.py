#
# Read in generated data, display & compare it
# Needs manual file paths and grid dimension settings for now...
# 
import sys
from manta import *
from helperInclude import *

# solver params
dim = 3
res = 64
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
gs = vec3(52,52,52)

# input file prefixes
basename1 = "t01a.py"
basename2 = "t01b.py"

# print info about running build, and those used to create data files?
buildInfo=1

# solver setup
s = Solver(name='main', gridSize = gs, dim=dim)
flags    = s.create(FlagGrid)
real1    = s.create(RealGrid)
real2    = s.create(RealGrid)
realErr  = s.create(RealGrid)

#ls1      = s.create(LevelsetGrid)
#ls2      = s.create(LevelsetGrid)
#lsErr    = s.create(RealGrid)

mac1     = s.create(VecGrid)
mac2     = s.create(VecGrid)
macErr   = s.create(VecGrid)

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
    
# use numbered files? or just base name?
appendNumber = False
# try to load uni file if it exists
def tryToLoad( grid, basename, suffix, number ):
	if(appendNumber==True):
		suffix = suffix+("_04d" % number )
	rfile = referenceFilename( basename, suffix ) 
	print("Trying to load " + rfile)
	if(os.path.isfile(rfile)):
		grid.load(rfile)
		if(buildInfo==1):
			printUniFileInfoString(rfile) # more detailed build info
	else:
		grid.clear()
	return 1

if(buildInfo==1):
	printBuildInfo() # more detailed build info , about what's running

# to be initialized later on...
realErrMax = 0
macErrMax = 0
lsErrMax = 0
partErrMax = 0

#main loop
for t in range(150):

	if(1):
		tryToLoad( real1, basename1, "dens"  , t )
		tryToLoad( real2, basename2, "dens"  , t )
		realErr.sub(real1,real2);
		realErrMax = gridMaxDiff(real1, real2)
	
		#realErr.print(zSlice=15) 
		print("Max difference in step " +str(t) + " = "+ str(realErrMax) )

	if(0):
		tryToLoad( mac1, basename1, "vel"  , t )
		tryToLoad( mac2, basename2, "vel"  , t )
		macErr.sub(mac1,mac2);
		macErrMax = gridMaxDiffVec3(mac1, mac2)
	
		#macErr.print(zSlice=15) 
		print("Max vec3 difference in step " +str(t) + " = "+ str(macErrMax) )

	# load particles
	if(0):
		tryToLoad( parts1 , basename1, "parts"  , t )
		tryToLoad( parts2 , basename2, "parts"  , t )

	#tryToLoad( pDens1 , basename1, ("pdens_%04d"  % t) )

	s.step()
    
