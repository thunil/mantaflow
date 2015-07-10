#
# Flip scene with particle seeding only in a narrow band around the surface
# 
from manta import *
import sys, os, time, random
random.seed(5481232)

# Command line arguments
GUI = True
if len(sys.argv)>1: GUI = (sys.argv[1] in ['True', 'true', '1'])

simtypeno = 3
if len(sys.argv)>2: simtypeno = int(sys.argv[2])

# 0: Dropping stuff, 
# 1: Cylinders, 
# 2: Pour
scene = 0
if len(sys.argv)>3: scene = int(sys.argv[3])


# Configuration 
if scene == 0: simname = "abc128_test"
if scene == 1: simname = "cylinders_test"
if scene == 2: simname = "pour_test"
simtype = ["levelset", "flip0", "flip", "nbflip"][simtypeno]
narrowBand  = 3; # nbflip only: no. of cells around surface which contain particles
combineBand = 2; # nbflip only: no. of cells around surface which are influenced by particles
kernelType  = 1; # XflipY only: particle mapping kernel (1: d-linear, 2: smooth SPH kernel with support of 4^d grid points)
if simtype in ["levelset", "flip0", "flip"]: narrowBand = combineBand = 0
if simtype in ["levelset", "flip0"]: kernelType = 0


# solver params
dim = 3

if scene == 0: res = 128; gs = vec3(2*res,1.5*res,res)
if scene == 1: res = 100; gs = vec3(2*res,2*res,2*res)
if scene == 2: res = 128; gs = vec3(res,2*res,res)
if dim==2: gs.z = 1;
s = Solver(name='main', gridSize = gs, dim=dim)
s.print('Solver grid resolution is: %i x %i x %i' % (gs.x, gs.y, gs.z))

# Print sim type info
s.print("Sim Type:")
if simtype == "nbflip": s.print("NB%i, CB%i" % (narrowBand, combineBand))
if simtype in ["flip", "nbflip"]:  s.print(", kernelType %i" % kernelType)
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
mapWeights  = s.create(MACGrid)
mapWeights2 = s.create(MACGrid)

pp        = s.create(BasicParticleSystem) 
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
	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.2,1.0)) # basin
	phi.join( fluidBasin.computeLevelset() )
	
	dropcounter = 0

	def dropstuff(dropcounter_):
		s.print('DROPPING STUFF')

		dropphi = s.create(LevelsetGrid)
		dropmesh = s.create(Mesh)
		capitalletters = [chr(ord('A')+i) for i in range(0,26)]
		dropmesh.load( '../scenes/letters/' + capitalletters[dropcounter_] + '.obj')
		dropmesh.scale( vec3(gs.x/6.0) );
		dropx = [1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]
		dropz = [1, 0, 1, 2, 0, 1, 0, 2, 1, 0, 1]
		t = vec3(dropx[dropcounter_]/2, random.random(), dropz[dropcounter_]/2)
		dropmesh.offset( t*gs*vec3(0.2,0.45,0.2) + (1-t)*gs*vec3(0.8,0.53,0.8) );
		dropmesh.computeLevelset(dropphi, 2);
		dropphi.applyToGrid(vel, value=vec3(0,0,0), isoval=1)
		phi.join( dropphi ) 
		flags.updateFromLevelset(phi)
		sampleLevelsetWithParticles( phi=dropphi, flags=flags, parts=pp, discretization=2, randomness=0.4 )

elif scene == 1: # Cylinders	
	for cen in [vec3(0.5,1,1.75), 
			    vec3(0.6,1,1.50), 
			    vec3(0.7,1,1.25), 
			    vec3(0.8,1,1.00), 
			    vec3(0.9,1,0.75), 
			    vec3(1.0,1,0.50), 
				vec3(1.1,1,0.25)]:
		c = s.create(Cylinder, center=res*cen+vec3(0,1,0), radius=res*0.05, z=res*vec3(0, 1, 0))
		phiObs.join(c.computeLevelset())
		c.applyToGrid(flags, value=2)

	fluidBasin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.2,1.0)) # basin
	phi.join( fluidBasin.computeLevelset() )

	fluidBox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.30,0.6,0.30)) # basin
	phi.join( fluidBox.computeLevelset() )

	phi.difference(phiObs)

elif scene == 2: # Pour	
	flags.setConst(2)
	H = vec3(1,1,dim-2)
	
	c = s.create(Cylinder, center=gs*vec3(0.6,0.5,0.5), radius=gs.x*0.4-1, z=gs*vec3(0,0.5-1/gs.y,0))
	c.applyToGrid(flags, value=4)
	phiObs = c.computeLevelset()
	c = s.create(Box, p0=gs*vec3(0,0.7,0)+H, p1=gs*vec3(1,1,1)-H)
	c.applyToGrid(flags, value=4)
	phiObs.join( c.computeLevelset() )
	phiObs.multConst(-1)
	phiObs.addConst(0.5)

	def pour():
		if s.frame == 152:
			c = s.create(Box, p0=gs*vec3(0,0.7,0)+H, p1=gs*vec3(1,1,1)-H)
			c.applyToGrid(flags, value=2)
			c = s.create(Cylinder, center=gs*vec3(0.6,0.85,0.5), radius=gs.x*0.4-1, z=gs*vec3(0,0.15-1/gs.y,0))
			c.applyToGrid(flags, value=4)
	
		if s.frame > 140: return
		s.print('POURING')
		source  = s.create(Sphere, center=gs*vec3(0.3,0.85,0.5), radius=gs.x*0.15)
		source2 = s.create(Sphere, center=gs*vec3(0.3,0.85,0.5), radius=gs.x*0.15+1)
		phi.join( source.computeLevelset() )
		source2.applyToGrid( grid=vel , value=gs*vec3(.03,-.01,0) )
		flags.updateFromLevelset(phi)
		
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


#main loop
lastframe = 0
timestep = -1
while s.frame < [200, 250, 250][scene]:
	timestep = timestep + 1
	s.print('\n### Frame %i - Timestep %i ###\n' % (s.frame,timestep))
	
	maxVel = vel.getMaxValue()
	if (s.cfl < 1000): s.adaptTimestep( maxVel )

	# velocities are extrapolated at the end of each step
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=(scene in [1]) ) 
	advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
	flags.updateFromLevelset(phi)
	if simtype != "flip":
		advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

	# create approximate surface level set, resample particles
	gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
	unionParticleLevelset( pp, pindex, flags, gpi, phiParts , radiusFactor )

	if simtype == "nbflip":
		phi.addConst(1.); # shrink slightly
		phi.join( phiParts );
		extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True, flags=flags, ignoreWalls=True)
	elif simtype in ["flip0", "flip"]:
		phi.copyFrom( phiParts );
		extrapolateLsSimple(phi=phi, distance=4, inside=True, flags=flags, ignoreWalls=True)

	if simtype in ["flip0", "flip", "nbflip"]:
		extrapolateLsSimple(phi=phi, distance=3, flags=flags, ignoreWalls=True, copyIntoBnd=1)
		flags.updateFromLevelset(phi)

	# make sure we have velocities throught liquid region
	if simtype == "nbflip":
		mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights, kernelType=kernelType );
		mapWeights2.copyFrom(mapWeights)
		extrapolateMACFromWeight( vel=velParts , distance=2, weight=mapWeights2 ) 
		combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, phi=phi, narrowBand=combineBand, thresh=0.0001)
		velOld.copyFrom(vel)
	elif simtype == "flip":
		mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights, kernelType=kernelType  );
		extrapolateMACFromWeight( vel=vel , distance=2, weight=mapWeights ) 
		
	if scene == 2: pour() # Pour

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
		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )
		
	if scene == 0: # Dropping stuff
		timings.disable()
		if s.frame >= dropcounter*10 and dropcounter < 9:
			dropstuff(dropcounter)
			dropcounter = dropcounter + 1
		timings.enable()
		
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

	if dim==3:
		timings.disable()
		phiMesh.copyFrom(phi)
		phiMesh.difference(phiObs)
		phiMesh.createMesh(mesh)
		timings.enable()
		
	timings.display()
	s.printMemInfo()
	s.print('%i particles, %i bytes each, %i MB total' % (pp.size(),(12+4+12),pp.size()*(12+4+12)/(1<<20)))
	s.step()
		
	# optionally particle data , or screenshot, or stats
	if 0 and s.frame!=lastframe:
		timings.disable()
		mesh.save( outdir + 'mesh/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )
				
		if simtype == "nbflip": 
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1, inc=10)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_preview_%04d.bobj.gz' % (s.frame-1) )
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )
		else:
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1, phi=phiMesh, minPhi=-1.5, inc=10)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_preview_%04d.bobj.gz' % (s.frame-1) )
			pp.createSphereMesh(meshspheres, radius=0.4, sphereQual=1, phi=phiMesh, minPhi=-1.5)
			meshspheres.save( outdir + 'meshspheres/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )

		phiMesh.copyFrom(phiParts)
		phiMesh.difference(phiObs)
		phiMesh.createMesh(meshparts)
		meshparts.save( outdir + 'meshparts/fluidsurface_final_%04d.bobj.gz' % (s.frame-1) )
				
		#if GUI:
		#    gui.screenshot( '../vid/frames/flip05_st%i_nb%02i_cb%02i_kt%i_%04d.png' % (simtypeno,round(narrowBand*10),round(combineBand*10),kernelType, s.frame) )

		lastframe = s.frame
		timings.enable()


timings.saveMean(outdir + 'meantimings.txt')