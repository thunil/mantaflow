#
# Inverted test, this should "fail", ie give a large difference...
# 

from manta import *
from helperInclude import *

# solver params
gs  = vec3(17, 177, 27)
s   = Solver(name='main', gridSize = gs, dim=3)


# prepare grids
density = s.create(RealGrid)
dummy   = s.create(RealGrid)

# set some value != 0
setConstant    ( density,  25.01 )
setConstant    ( dummy  , -25.00 )

# verify , note - this should fail!
if (getGenRefFileSetting()==1):
	doTestGrid( __file__, "dens" , s, density )
else:
	doTestGrid( __file__, "dens" , s, dummy , threshold=50. , invertResult=True )


