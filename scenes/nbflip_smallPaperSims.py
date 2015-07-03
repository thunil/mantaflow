#
# Flip scene with particle seeding only in a narrow band around the surface
# 
from manta import *
import sys, os, time, random
random.seed(5481232)

# Command line arguments
GUI = True
if len(sys.argv)>1: GUI = (sys.argv[1] in ['True', 'true', '1'])

simtypeno = 2
if len(sys.argv)>2: simtypeno = int(sys.argv[2])

# 0:
scene = 0
if len(sys.argv)>3: scene = int(sys.argv[3])


# Configuration 
if scene == 0: simname = "2dtest"
simtype = ["levelset", "flip0", "flip", "nbflip","nbflipd"][simtypeno]
narrowBand  = 4; # nbflip only: no. of cells around surface which contain particles
combineBand = 3; # nbflip only: no. of cells around surface which are influenced by particles
kernelType  = 1; # XflipY only: particle mapping kernel (1: d-linear, 2: smooth SPH kernel with support of 4^d grid points)
if simtype in ["levelset", "flip0", "flip"]: narrowBand = combineBand = 0
if simtype in ["levelset", "flip0"]: kernelType = 0

# solver params
dim = 2

if scene == 0: res = 64; gs = vec3(res,res,1)
s = Solver(name='main', gridSize = gs, dim=dim)
s.print('Solver grid resolution is: %i x %i x %i' % (gs.x, gs.y, gs.z))

# Print sim type info
s.print("Sim Type: %s" % simtype)
if simtype == "nbflip": s.print("NB%i, CB%i" % (narrowBand, combineBand))
if simtype in ["flip", "nbflip", "nbflipd"]:  s.print(", kernelType %i" % kernelType)
s.print("")

# adaptive time stepping
s.frameLength = 1.0   # length of one frame (in "world time")
s.timestep    = 0.5
s.timestepMin = 0.2   # time step range
s.timestepMax = 1.0
s.cfl         = 5.0   # maximal velocity per cell, 0 to use fixed timesteps

gravity = (0,-0.003,0)

minParticles = pow(2,dim)
timings = Timings()

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)

phiParts = s.create(LevelsetGrid)
phi      = s.create(LevelsetGrid)
phiMesh  = s.create(LevelsetGrid)
ttt      = s.create(RealGrid)
perCellCorr = s.create(RealGrid)
pressure = s.create(RealGrid)
phiBackup = s.create(RealGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
velParts = s.create(MACGrid)
velDiff  = s.create(MACGrid)
mapWeights = s.create(MACGrid)

pp        = s.create(BasicParticleSystem) 
pphi      = pp.create(PdataReal) 
pVel      = pp.create(PdataVec3) 
mesh      = s.create(Mesh)
meshparts = s.create(Mesh)
meshspheres = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# geometry in world units (to be converted to grid space upon init)
flags.initDomain(boundaryWidth=0)
phi.initFromFlags(flags)

phiObs = s.create(Box, p0=vec3(1,1,1), p1=gs-1).computeLevelset()
phiObs.multConst(-1)

if scene == 0: # Dropping stuff
	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.3,1.0)) # basin
	phi.join( fluidBasin.computeLevelset() )
	fluidSphere = s.create(Sphere, center=gs*vec3(0.5,0.3,0.5), radius=gs.x*0.2) # basin
	phi.join( fluidSphere.computeLevelset() )

	dropcounter = 0
			
flags.updateFromLevelset(phi)
targetVolume = calcFluidVolume(flags)


if simtype in ["flip0", "flip", "nbflip", "nbflipd"]:
	sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.4 )

if simtype in ["flip", "nbflip", "nbflipd"]:
	mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

if 1 and (GUI):
	gui = Gui()
	gui.show( dim==2 )
	gui.pause()
	  
	# show all particles shaded by velocity
	if(dim==2):
		gui.nextPdata()
		#gui.nextPartDisplay()
		#gui.nextPartDisplay()

outdir = '../sim/'+simname+'/'+simtype+'/'
for subdir in ['mesh', 'meshparts', 'meshspheres']:
	os.makedirs(outdir + subdir, exist_ok=True)


#main loop
lastframe = 0
timestep = -1
while s.frame < 2000:
	timestep = timestep + 1
	s.print('\n### Frame %i - Timestep %i ###\n' % (s.frame,timestep))
	
	maxVel = vel.getMaxValue()
	if (s.cfl < 1000): s.adaptTimestep( maxVel )

	# velocities are extrapolated at the end of each step
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=(scene in [1,2]) ) 
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
	flags.updateFromLevelset(phi)
	#if simtype != "flip":
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts , radiusFactor )

	if simtype in ["nbflip", "nbflipd"]:
		phi.addConst(1.); # shrink slightly
		phi.join( phiParts );
		extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True, flags=flags, ignoreWalls=True)
		#phi.reinitMarching(flags, maxTime=narrowBand+2, ignoreWalls=True)
	elif simtype in ["flip0", "flip"]:
		phi.copyFrom( phiParts );
		extrapolateLsSimple(phi=phi, distance=4, inside=True, flags=flags, ignoreWalls=True)
		#phi.reinitMarching(flags, maxTime=4, ignoreWalls=True)

	if simtype in ["flip0", "flip", "nbflip", "nbflipd"]:
		extrapolateLsSimple(phi=phi, distance=3, flags=flags, ignoreWalls=True)
		flags.updateFromLevelset(phi)

	# make sure we have velocities throught liquid region
	if simtype in ["nbflip"]:
		mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights, kernelType=kernelType );
		velDiff.copyFrom(vel); velDiff.multConst(vec3(-1)) 
		combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, phi=phi, narrowBand=combineBand, thresh=0.0001)
		velDiff.add(vel)
		velOld.copyFrom(vel)
	elif simtype in ["nbflipd"]:
		mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights, kernelType=kernelType );
		velDiff.copyFrom(vel); velDiff.multConst(vec3(-1)) 
		combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, thresh=1)
		velDiff.add(vel)
		velOld.copyFrom(vel)
	elif simtype == "flip":
		velDiff.copyFrom(vel); velDiff.multConst(vec3(-1)) 
		mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights, kernelType=kernelType  );
		velDiff.add(vel)
		
	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=gravity)
	setWallBcs(flags=flags, vel=vel)
	vol = calcFluidVolume(flags)
	perCellCorr.setConst(0.02 * ( targetVolume-vol ) / vol)
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi, perCellCorr=perCellCorr)
	setWallBcs(flags=flags, vel=vel)

	# make sure we have proper velocities for levelset advection
	if simtype in ["flip", "nbflip", "nbflipd"]:
		extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)) )
		#extrapolateMACFromWeight( vel=vel , distance=(int(maxVel*1.1 + 2.)), weight=tmpVec3 ) 
	elif simtype in ["levelset", "flip0"]:
		phiBackup.copyFrom(phi)
		phi.reinitMarching(flags=flags, velTransport=vel, ignoreWalls=False);	
		phi.copyFrom(phiBackup);
		phi.reinitExact(flags=flags)

	if simtype in ["flip", "nbflip"]:
		#updateParticleDistance(parts=pp, partPhi=pphi, phi=phi, t=1)
		#flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95, partPhi=pphi, pphiPIC=-combineBand, pphiFLIP=0 )
		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )
	if simtype in ["nbflipd"]:
		updateParticleDistance(parts=pp, partPhi=pphi, phi=phi, t=0.666)
		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95, partPhi=pphi, pphiPIC=-combineBand, pphiFLIP=0 )

	if scene == 2: pour() # Pour				
		
	# set source grids for resampling, used in adjustNumber!
	if simtype in ["flip", "nbflip", "nbflipd"]:
		pVel.setSource( vel, isMAC=True )
	else:
		pVel.setSource( 0, isMAC=False )

	nrg = calcTotalEnergy(flags,vel,gravity)

	if simtype in ["nbflip", "nbflipd"]:
		adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor , narrowBand=narrowBand ) 
	elif simtype in ["flip0", "flip"]:
		adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

	#timings.display()
	#s.printMemInfo()
	#s.print('%i particles, %i bytes each, %i MB total' % (pp.size(),(12+4+12),pp.size()*(12+4+12)/(1<<20)))
	s.step()
		
	# optionally particle data , or screenshot, or stats
	if 0 and s.frame!=lastframe:
		if simtype in ["nbflip", "nbflipd"]: 
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1, inc=10)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_preview_%04d.bobj.gz' % (s.frame-1) )
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )
		else:
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1, phi=phiMesh, minPhi=-1.5, inc=10)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_preview_%04d.bobj.gz' % (s.frame-1) )
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1, phi=phiMesh, minPhi=-1.5)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )

		#if GUI:
		#    gui.screenshot( '../vid/frames/flip05_st%i_nb%02i_cb%02i_kt%i_%04d.png' % (simtypeno,round(narrowBand*10),round(combineBand*10),kernelType, s.frame) )

		lastframe = s.frame


timings.saveMean(outdir + 'meantimings.txt')