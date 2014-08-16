#!/usr/bin/python
#
# Simple script to run all base tests
# The test scripts are named test_XXXX_description.py
# They are assumed to be in ascending order or complexity.
#
# Rough ordering:
# 0xxx tests are very basic (mostly single operator calls)
# 1xxx is for 2d sims
# 2xxx for "real" 3d sims
# xl... larger 3d sims (similar to 2xxx cases)
# 

import os
import shutil
import sys
import re
from subprocess import check_output
from helperGeneric import *

# note - todo, right now this script assumes the test_X.py files are
# in the current working directory... this should be changed at some 
# point.

# debugging, print outputs from all manta calls
printAllOutpus = 0
filePrefix = "test_"

if(len(sys.argv)<2):
	print "Usage runTests.py <manta-executable>"
	exit(1)

manta = sys.argv[1]
print "Using mantaflow executable '" + manta + "' " 

# check parameters:
# xl - whether complex/large tests should be run
if(len(sys.argv)>2):
	for i in range(2, len(sys.argv)):
		if(sys.argv[i].lower() == "xl"):
			print "Running XL tests! This can take a while..."
			filePrefix = "xltest_"
		else:
			print "Error: unknown parameter '" + sys.argv[i] +"' "
			exit(1);

files = os.popen("ls "+str(filePrefix)+"????_*.py").read() 
#print "Debug - using test scene files: "+files

# NT_DEBUG , todo use , extract file name!
datadir = "./data"
if not os.path.exists( datadir ):
	os.makedirs( datadir )	

if getGenRefFileSetting():
	print "\nNote - generating test data for all tests!"
	print "Tests results will not be evaluated...\n"

currdate = os.popen("date \"+%y%m%d%H%M\"").read() 
currdate = "201408162205" # NT_DEBUG
# in visual mode, also track runtimes
if getVisualSetting():
	dirname = "./runtimes"
	if not os.path.exists( dirname ):
		os.makedirs( dirname )	

num = 0
numOks = 0
numFail = 0
failedTests = ""

files = files.split('\n')
for file in files:
	if ( len(file) < 1):
		continue

	(utime1, stime1, cutime1, cstime1, elapsed_time1) = os.times() 

	num += 1
	print "Running '" + file + "' "
	result = os.popen(manta + " " + file + " 2>&1 ").read() 

	(utime2, stime2, cutime2, cstime2, elapsed_time2) = os.times() 

	oks = re.findall(r"\nOK!", result)
	#print oks
	numOks += len(oks)

	fails = re.findall(r"\nFAIL!", result)
	# also check for "Errors"
	if (len(fails)==0) :
		# note - if there are errors, the overall count of tests won't be valid anymore... 
		fails = re.findall(r"Error", result) 
	
	if (len(fails)!=0) :
		# we had failed tests! add to list...
		failedTests = failedTests+" "+file 
	#print fails
	numFail += len(fails)

	if (len(fails)>0) or (printAllOutpus==1):
		print
		print "Full output: " + result
		print

	# store benchmarking results (if theres any output)
	if getVisualSetting() and os.path.isfile( "%s_0001.ppm"%(file) ):
	#if getVisualSetting():
		runtime = elapsed_time2-elapsed_time1 
		timefile = "%s_v%d.time" % (file, getVisualSetting()) 
		
		#print("echo %s %f >> %s " % (currdate, runtime, timefile) ) # debug
		os.popen("echo %s %f >> %s " % (currdate, runtime, timefile) ).read() 

		os.popen("./helperGnuplot.sh %s"%(timefile)) 
		#exit(1);


if getGenRefFileSetting():
	print "Test data generated"
	exit(0)

if getVisualSetting():
	print "Viusal data generated"
	# now convert & remove all ppms"
	os.popen("for i in ./test_*.ppm ; do convert $i $(basename $i .ppm).png; done")
	outpngdir = "./result_%s"%(currdate)
	if not os.path.exists( outpngdir ):
		os.makedirs( outpngdir )	
	os.popen("mv ./test_*.png ./result_%s/"
	trashdir = "./trash"
	if not os.path.exists( trashdir ):
		os.makedirs( trashdir )	
	os.popen("mv ./test_*.ppm %s"%(trashdir)
	exit(0)

print
print " ============================================= "
print
print "Test summary: "  + str(numOks) + " passed, " + str(numFail) + " failed.   (from "+str(num) + " files) "

if (numFail==0) and (numOks==0):
	print "Failure! Perhaps manta executable didnt work, or runTests not called in test directory \n"
	exit(1)
elif (numFail==0) and (numOks>0):
	print "All good :) \n"
	exit(0)
else:
	print "Oh no :( the following tests failed: %s \n" % failedTests
	exit(2)


