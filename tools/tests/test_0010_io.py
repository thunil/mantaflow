#
# Very basic test, set constant values, mainly testing IO
# 

import sys
print ("Running python "+sys.version)

from manta import *
from helperInclude import *


# solver params
gs  = vec3(10, 20, 30)
s   = Solver(name='main', gridSize = gs, dim=3)

# prepare grids
density = s.create(RealGrid)
vel     = s.create(MACGrid)

# set some value != 0
setConstant    ( density, 0.123 )
setConstantVec3( vel    , vec3(0.1, 0.2, 0.3) )

# verify
doTestGrid( sys.argv[0], "dens" , s, density  )
doTestGrid( sys.argv[0], "vel"  , s, vel      )


