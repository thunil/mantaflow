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

def checkResult( result , thresh ):
	#print ("Checking , result=%f , thresh=%f" % (result , thresh) )
	if ( result <= thresh) :
		print "OK! Results match..."
		return 0
	else:
		print "FAIL! Allowed threshold "+str(thresh)+", results differ by "+str(result)
		return 1

def doTestReal( file , name, solver , grid, threshold=1e-10  ):
	# create temp grid
	tmp = solver.create(RealGrid)

	# check env var for generate data setting
	genRefFiles = int(os.getenv('MANTA_GEN_TEST_DATA', 0))
	#print "Gen-data-setting: " + str(genRefFiles)
	if (genRefFiles==1):
		#grid.save( outputFilename( file, name ) )
		#shutil.copyfile( outputFilename( file, name ) , referenceFilename( file, name ) )
		grid.save( referenceFilename( file, name ) )
		print "OK! Generated reference file '" + referenceFilename( file, name ) + "'"
		return 0
	else:
		tmp.load( referenceFilename( file, name ) )
		errVal = gridMaxDiff( grid, tmp )
		return checkResult( errVal , threshold )


def doTestVec3( file , name, solver , grid , threshold=1e-10 ):
	# create temp grid
	tmp = solver.create(MACGrid)

	# check env var for generate data setting
	genRefFiles = int(os.getenv('MANTA_GEN_TEST_DATA', 0))
	#print "Gen-data-setting: " + str(genRefFiles)
	if (genRefFiles==1):
		#grid.save( outputFilename( file, name ) )
		#shutil.copyfile( outputFilename( file, name ) , referenceFilename( file, name ) )
		grid.save( referenceFilename( file, name ) )
		print "OK! Generated reference file '" + referenceFilename( file, name ) + "'"
		return 0
	else:
		tmp.load( referenceFilename( file, name ) )
		errVal = gridMaxDiffVec3( grid, tmp )
		return checkResult( errVal , threshold )



