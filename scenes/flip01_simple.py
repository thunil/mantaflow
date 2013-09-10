#
# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

# solver params
dim = 2
particleNumber = 2
res = 64
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
	particleNumber = 3 # use more particles in 2d
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
phiInit  = s.create(LevelsetGrid)
pressure = s.create(RealGrid)
flip     = s.create(FlipSystem) 

# scene setup
flags.initDomain(boundaryWidth=0)
# enable one of the following
fluidbox = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # breaking dam
#fluidbox = s.create(Box, p0=gs*vec3(0.4,0.4,0.4), p1=gs*vec3(0.6,0.8,0.6)) # centered falling block
phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
flip.initialize( flags=flags, discretization=particleNumber, randomness=0.2 )
    
if (GUI):
    gui = Gui()
    gui.show()
    gui.pause()
    
#main loop
for t in range(2500):
    
    # FLIP 
    flip.advectInGrid(flaggrid=flags, vel=vel, integrationMode=IntRK4)
    flip.velocitiesToGrid(vel=vel, flags=flags)
    flip.markFluidCells(flags=flags)
    
    addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)

    # we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
    extrapolateMACSimple( flags=flags , vel=vel )
    
    # FLIP velocity update
    flip.velocitiesFromGrid(vel=vel, flags=flags, flipRatio=0.97)
    
    #gui.screenshot( 'flipt_2rk4_%04d.png' % t );
    s.step()

