#
# Flip scene with particle seeding only in a narrow band around the surface
# 
from manta import *
import sys, os, time

# Command line arguments
GUI = True
if len(sys.argv)>1: GUI = (sys.argv[1] in ['True', 'true', '1'])

simtypeno = 2
if len(sys.argv)>2: simtypeno = int(sys.argv[2])

# Configuration 
simname = "bdam3"
simtype = ["levelset", "flip0", "flip", "nbflip"][simtypeno]
narrowBand  = 4; # nbflip only: no. of cells around surface which contain particles
combineBand = 3; # nbflip only: no. of cells around surface which are influenced by particles
kernelType  = 1; # XflipY only: particle mapping kernel (1: d-linear, 2: smooth SPH kernel with support of 4^d grid points)
if simtype in ["levelset", "flip0", "flip"]: narrowBand = combineBand = 0
if simtype in ["levelset", "flip0"]: kernelType = 0


print("Sim Type:", simtype, end='')
if simtype == "nbflip": print(", NB", narrowBand, ", CB", combineBand, end='')
if simtype in ["flip", "nbflip"]:  print(", kernelType", kernelType, end='')
print("")

# solver params
dim = 3
res = 96
gs = vec3(4*res,2*res,res)
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# adaptive time stepping
s.frameLength = 2.0   # length of one frame (in "world time")
s.timestep    = 1.0
s.timestepMin = 0.2   # time step range
s.timestepMax = 2.0
s.cfl         = 3.0   # maximal velocity per cell, 0 to use fixed timesteps

gravity = (0,-0.003,0)

minParticles = pow(2,dim)
timings = Timings()

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)

phiParts = s.create(LevelsetGrid)
phi      = s.create(LevelsetGrid)
phiSpheres = s.create(LevelsetGrid)
ttt      = s.create(RealGrid)
perCellCorr = s.create(RealGrid)
pressure = s.create(RealGrid)
phiBackup = s.create(RealGrid)

vel      = s.create(MACGrid)
combineDiff = s.create(MACGrid)
velOld   = s.create(MACGrid)
velParts = s.create(MACGrid)
mapWeights = s.create(MACGrid)

pp        = s.create(BasicParticleSystem) 
pVel      = pp.create(PdataVec3) 
mesh      = s.create(Mesh)
meshparts = s.create(Mesh)
meshspheres = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem) 
gpi    = s.create(IntGrid)

# geometry in world units (to be converted to grid space upon init)
flags.initDomain(boundaryWidth=1)

#obssphere = s.create(Box, p0=gs*vec3(0.2,0.4,0), p1=gs*vec3(0.4,0.6,1))
#obssphere.applyToGrid(flags, value=2)
domainbox = s.create(Box, p0=vec3(1,1,1), p1=gs-1)
phiDomain = domainbox.computeLevelset()

# breaking dam & basin
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.2,0.95,1)) # breaking dam
fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.2,1.0)) # basin

phi.initFromFlags(flags)
phi.join( fluidbox.computeLevelset() )
phi.join( fluidBasin.computeLevelset() )
	
flags.updateFromLevelset(phi)

#if simtype in ["flip0", "flip", "nbflip"]:
sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.4 )

if simtype in ["flip", "nbflip"]:
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


#statstxt = open(outdir + 'stats.txt', 'w')

#main loop
lastframe = 0
timestep = -1
while s.frame < 200:
	timestep = timestep + 1
	print('\n### Frame %i - Timestep %i ###\n' % (s.frame,timestep))
	
	if (s.cfl < 1000): s.adaptTimestep( vel.getMaxValue() )
	
	# velocities are extrapolated at the end of each step
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
	flags.updateFromLevelset(phi)
	advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts , radiusFactor )

	if simtype == "nbflip":
		phi.addConst(1.); # shrink slightly
		phi.join( phiParts );
		extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True)
	elif simtype in ["flip0", "flip"]:
		phi.copyFrom( phiParts );
		extrapolateLsSimple(phi=phi, distance=4, inside=True)

	#extrapolateLsSimple(phi=phiParts, distance=4, inside=True)
		
	if simtype in ["flip0", "flip", "nbflip"]:
		extrapolateLsSimple(phi=phi, distance=3)
		flags.updateFromLevelset(phi)

	# make sure we have velocities throught liquid region
	if simtype == "nbflip":
		mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights, kernelType=kernelType );
		combineDiff.copyFrom(vel)
		combineDiff.multConst(vec3(-1))
		combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, phi=phi, narrowBand=combineBand, thresh=0.1)
		combineDiff.add(vel)
		velOld.copyFrom(vel)
	elif simtype == "flip":
		mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights, kernelType=kernelType  );
		
	# forces & pressure solve
	addGravity(flags=flags, vel=vel, gravity=gravity)
	setWallBcs(flags=flags, vel=vel)
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
	setWallBcs(flags=flags, vel=vel)

	# make sure we have proper velocities for levelset advection
	if simtype in ["flip", "nbflip"]:
		extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)) )
		#extrapolateMACFromWeight( vel=vel , distance=(int(maxVel*1.1 + 2.)), weight=tmpVec3 ) 
	elif simtype in ["levelset", "flip0"]:
		phiBackup.copyFrom(phi)
		phi.reinitMarching(flags=flags, velTransport=vel, ignoreWalls=False);	
		phi.copyFrom(phiBackup);
		phi.reinitExact(flags=flags)

	if simtype in ["flip", "nbflip"]:	
		#flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )
		#flipVelocityUpdateNb(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95, narrowBand=narrowBand , phi=phi, test=ttt )
		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )

	# set source grids for resampling, used in adjustNumber!
	if simtype in ["flip", "nbflip"]:
		pVel.setSource( vel, isMAC=True )
	else:
		pVel.setSource( 0, isMAC=False )

	#kin = calcKineticEnergy(flags,vel)
	vol = calcFluidVolume(flags)
	nrg = calcTotalEnergy(flags,vel,gravity)

	if simtype == "nbflip":
		adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor , narrowBand=narrowBand ) 
	else:
		adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 
		
	timings.display()
	#s.printMemInfo()
	s.step()
		
	# optionally particle data , or screenshot, or stats
	if 1 and s.frame!=lastframe:
		#pp.save( outdir + 'parts/parts_%04d.uni' % (s.frame-1) )
		
		phi.createMesh(mesh)
		mesh.save( outdir + 'mesh/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )
		
		phiParts.createMesh(meshparts)
		meshparts.save( outdir + 'meshparts/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )
		
		if simtype == "nbflip": pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1)
		else:
			phiSpheres.copyFrom(phiDomain)			
			phiSpheres.intersect(phi)
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1, phi=phiSpheres, minPhi=-2)
		meshspheres.save( outdir + 'meshspheres/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )
				
		#statstxt.write('%.6f %i\n' % (nrg,vol))
		
		#if GUI:
		#    gui.screenshot( '../vid/frames/flip05_st%i_nb%02i_cb%02i_kt%i_%04d.png' % (simtypeno,round(narrowBand*10),round(combineBand*10),kernelType, s.frame) )

		lastframe = s.frame

timings.saveMean(outdir + 'meantimings.txt')
#statstxt.close()
	

	



