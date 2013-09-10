#
# Testing shapes
# 

from manta import *
from helperInclude import *

# solver params
res = 42
gs  = vec3(res,res,res)
s   = Solver(name='main', gridSize = gs, dim=3)


# prepare grids

density = s.create(RealGrid)
vel     = s.create(MACGrid)

setConstant    ( density, 0. )
setConstantVec3( vel    , vec3(0, 0, 0) )

sph1 = s.create(Sphere, center=gs*vec3(0.3,0.4,0.5), radius=res*0.2)
sph1.applyToGrid(grid=density, value=0.302)

velVal = vec3(0.1,0.1,0.4)
sph2 = s.create(Sphere, center=gs*vec3(0.6,0.5,0.4), radius=res*0.25)
sph2.applyToGrid(grid=vel, value=velVal )

doTestGrid( __file__, "densSph" , s, density  )
doTestGrid( __file__, "velSph"  , s, vel      )


# sphere

setConstant    ( density, 0. )
setConstantVec3( vel    , vec3(0, 0, 0) )

box1 = s.create(Box, p0=gs*vec3(0.2,0.2,0.3), p1=gs*vec3(0.9,0.8,0.9) )
box1.applyToGrid(grid=density, value=0.812)

velVal = vec3(0.5,0.1,0.1)
box2 = s.create(Box, p0=gs*vec3(0.2,0.2,0.3), p1=gs*vec3(0.9,0.8,0.9) )
box2.applyToGrid(grid=vel, value=velVal )

doTestGrid( __file__, "densBox" , s, density  )
doTestGrid( __file__, "velBox"  , s, vel      )


# cylinder

setConstant    ( density, 0. )
setConstantVec3( vel    , vec3(0, 0, 0) )

cyl1 = s.create(Cylinder, center=gs*vec3(0.5,0.5,0.5), radius=res*0.2, z=gs*vec3(0, 0.3, 0))
cyl1.applyToGrid(grid=density, value=0.432)

velVal = vec3(0.4,0.3,0.2)
cyl2 = s.create(Cylinder, center=gs*vec3(0.5,0.5,0.5), radius=res*0.2, z=gs*vec3(0, 0.3, 0))
cyl2.applyToGrid(grid=vel, value=velVal )

doTestGrid( __file__, "densCyl" , s, density  )
doTestGrid( __file__, "velCyl"  , s, vel      )


