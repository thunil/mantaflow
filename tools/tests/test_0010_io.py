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
flags   = s.create(IntGrid)

# set some value != 0
density.setConst( 0.123 )
vel.setConst( vec3(0.1, 0.2, 0.3) )
flags.setConst( 7193 )

#density.printGrid( zSlice=15 ) # debug info

# verify
doTestGrid( sys.argv[0], "dens" , s, density  , threshold=1e-08 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "vel"  , s, vel      , threshold=1e-08 , thresholdStrict=1e-14  )
doTestGrid( sys.argv[0], "flag" , s, flags    , threshold=1e-14 , thresholdStrict=1e-14  )

doGenerateInfo( sys.argv[0] )

