#
# Simple example for free-surface simulation with MacCormack advection
#

from manta import *
import sys

l = len(sys.argv)
openB = sys.argv[1] if l>=2 else 'xX'
bWidth = int(sys.argv[2]) if l>=3 else 1
factor = vec3(int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])) if l>=6 else vec3(1,1,1)

# solver params
dim = 2
res = 64
gs = vec3(factor.x*res,factor.y*res,factor.z*res)
if (dim==2):
    gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.25
ghostFluid = True
accuracy = 5e-5

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

# scene setup
flags.initDomain(boundaryWidth=bWidth)
p0=vec3(0,0,0)
p1=vec3(1,0.2,1)
a=vec3(0.5,0.5,0.5)*gs+p0-vec3(0.5,0.5,0.5)*(gs/factor)
b=a+(p1-p0)*(gs/factor)
basin = s.create(Box, p0=a, p1=b)
drop  = s.create(Sphere, center=vec3(0.5,0.5,0.5)*gs, radius=res*0.15)
phi = basin.computeLevelset()
phi.join(drop.computeLevelset())
flags.updateFromLevelset(phi)

setOpenBound(flags,bWidth,openB,FlagOutflow|FlagEmpty) 

        
if (GUI):
    gui = Gui()
    gui.show()
    #gui.pause()
    

#main loop
for t in range(500):
    
    # update and advect levelset
    phi.reinitMarching(flags=flags, velTransport=vel, ignoreWalls=True)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1) 
    resetOutflow(flags=flags,phi=phi)
    flags.updateFromLevelset(phi)
    
    # velocity self-advection
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2, openBounds=True, depth=bWidth+1)
    addGravity(flags=flags, vel=vel, gravity=vec3(0,(-0.025/max(factor.x,max(factor.y,factor.z))),0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)
    if ghostFluid:
        solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, phi=phi )
    else:
        solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy)
    setWallBcs(flags=flags, vel=vel)
    
    # note: these meshes are created by fast marching only, should smooth
    #       geometry and normals before rendering (only in 3D for now)
    if (dim==3):
        phi.createMesh(mesh)
        #mesh.save('phi%04d.bobj.gz' % t)
    
    s.step()
    gui.screenshot( 'freesurface_'+openB+'_'+str(bWidth)+'_('+str(int(factor.x))+','+str(int(factor.y))+')_%04d.png' % t)
    
    #gui.pause()
    #gui.screenshot( 'freesurface_xX_%04d.png' % t )



