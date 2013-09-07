#
# Very basic test, set constant values, mainly testing IO
# 

from manta import *
from helperInclude import *

# solver params
gs  = vec3(10, 20, 30)
s   = Solver(name='main', gridSize = gs, dim=3)

# flags not needed for now
#flags   = s.create(FlagGrid)
#flags.initDomain()
#flags.fillGrid()


# prepare grids

density = s.create(RealGrid)
vel     = s.create(MACGrid)

# set some value != 0
setConstant    ( density, 0.123 )
setConstantVec3( vel    , vec3(0.1, 0.2, 0.3) )

# verify
doTestReal( __file__,"dens" , s, density  )
doTestVec3( __file__,"vel"  , s, vel      )


