#
# Helper functions for test runs in mantaflow
# 

from manta import *
import os
import shutil




def outputFilename( file, gridname ):
	return file +"_"+ gridname + "_out.uni" 

def referenceFilename( file, gridname ):
	return file +"_"+ gridname + "_ref.uni" 


def checkResult( name, result , thresh, invertResult=False ):
	print ("Checking '%s', result=%f , thresh=%f" % ( name , result , thresh) )

	if ( ( result > 0) and (result < 1e-04) ):
		print ("Debug, small difference: scaled output %f" % ( result * 100000.0 ) ) # debugging...

	allGood = 0
	if ( result <= thresh) :
		allGood = 1

	# for checks that should fail
	if ( invertResult == True) :
		if ( allGood == 0) :
			allGood = 1
		else:
			allGood = 0

	# now react on outcome...
	if ( allGood == 1 ):
		print("OK! Results for "+name+" match...")
		return 0
	else:
		print("FAIL! Allowed "+name+" threshold "+str(thresh)+", results differ by "+str(result))
		return 1


def getGenRefFileSetting( ):
	# check env var for generate data setting
	ret = int(os.getenv('MANTA_GEN_TEST_DATA', 0))
	# print("Gen-data-setting: " + str(ret))
	return ret


# note, typeId encode what type of grid we have - 
# 
def doTestGrid( file , name, solver , grid, threshold=1e-10, invertResult=False  ):

	# handle grid types that need conversion
	#print( "doTestGrid, incoming grid type :" + type(grid).__name__)
	if ( type(grid).__name__ == "MACGrid" ):
		gridTmpMac = solver.create(VecGrid)
		convertMacToVec3(grid , gridTmpMac )
		return doTestGrid( file, name, solver, gridTmpMac , threshold)
	if ( type(grid).__name__ == "LevelsetGrid" ):
		gridTmpLs = solver.create(RealGrid)
		convertLevelsetToReal(grid , gridTmpLs )
		return doTestGrid( file, name, solver, gridTmpLs  , threshold)
	if ( type(grid).__name__ == "IntGrid" ):
		print( "Error doTestGrid - int grids not yet supported...")
		return 1

	# now we should only have real & vec3 grids

	# create temp grid
	if ( type(grid).__name__ == "RealGrid" ):
		compareTmpGrid = solver.create(RealGrid)
	elif ( type(grid).__name__ == "VecGrid" ):
		compareTmpGrid = solver.create(VecGrid)
	else:
		print( "Error doTestGrid - unknown grid type " + type(grid).__name__ )
		return 1

	genRefFiles = getGenRefFileSetting()

	if (genRefFiles==1):
		#grid.save( outputFilename( file, name ) )
		#shutil.copyfile( outputFilename( file, name ) , referenceFilename( file, name ) )
		grid.save( referenceFilename( file, name ) )
		print( "OK! Generated reference file '" + referenceFilename( file, name ) + "'")
		return 0
	else:
		compareTmpGrid.load( referenceFilename( file, name ) )

		errVal = 1e10
		if ( type(grid).__name__ == "RealGrid" ):
			errVal = gridMaxDiff    ( grid, compareTmpGrid )
		elif ( type(grid).__name__ == "VecGrid" ):
			errVal = gridMaxDiffVec3( grid, compareTmpGrid )

		# finally, compare max error to allowed threshold, and return result
		return checkResult( name, errVal , threshold , invertResult )




