#!/usr/bin/python
#
# Simple script to run all base tests
# 

import os
import shutil
import sys
import re
from subprocess import check_output

# debugging, print outputs from all manta calls
printAllOutpus = 0

if(len(sys.argv)<2):
	print "Usage runTests.py <manta-executable>"
	exit(1)

manta = sys.argv[1]
print "Using mantaflow executable '" + manta + "' " 

files = os.popen("ls test_????_*.py").read() 
#print files

genRefFiles = int(os.getenv('MANTA_GEN_TEST_DATA', 0))
if (genRefFiles==1):
	print "\nNote - generating test data for all tests!"
	print "Tests results will not be evaluated...\n"

num = 0
numOks = 0
numFail = 0

files = files.split('\n')
for file in files:
	if ( len(file) < 1):
		continue
	num += 1
	print "Running '" + file + "' "
	result = os.popen(manta + " " + file + " 2>&1 ").read() 

	oks = re.findall(r"OK!", result)
	#print oks
	numOks += len(oks)

	fails = re.findall(r"FAIL!", result)
	# also check for "Errors"
	if (len(fails)==0) :
		fails = re.findall(r"Error", result) 
	#print fails
	numFail += len(fails)

	if (len(fails)>0) or (printAllOutpus==1):
		print
		print "Full output: " + result
		print


if (genRefFiles==1):
	print "Test data generated"
	exit(0)

print
print "Test summary, " +str(num) + " runs, " + str(numOks) + " passed, " + str(numFail) + " failed."

if (numFail==0) and (numOks==0):
	print "Failure! Perhaps manta executable didnt work, or runTests not called in test directory \n"
	exit(1)
elif (numFail==0) and (numOks>0):
	print "All good :) \n"
	exit(0)
else:
	print "Some tests failed :( \n"
	exit(2)


