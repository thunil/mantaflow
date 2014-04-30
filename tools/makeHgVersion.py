#!/usr/bin/python
import os
import shutil
import sys
import re
doDebug = False

# helper function to write output file
def writeHeader( filename, content ):
	try:
		outfile = open(filename, "w")
		outfile.write( content )
		outfile.close()
		if(doDebug):
			print( "Wrote '" + filename +"' " )
	except IOError:
		print("Warning, unable to write file '"+filename+"' ")
		exit(1)


# params

if(len(sys.argv)<2):
	print("Usage makeHgVersion.py <out-file> <optional: path-to-hg> ")
	print("Warning, the target file <out-file> will be overwritten! ")
	exit(1)

# target file
outname = sys.argv[1]

# path to git/hg executable, try a few options
# note /opt/local/bin/xxx is a double entry, can be overwritten by command line arg
exenames = [ "--replace--", "--replace--", "/opt/local/bin/git", "/usr/local/bin/git" ]
# check default
exenames[1] = os.popen("which git").read() 
exenames[1] = exenames[1].rstrip('\n')
# optionally, make argument
if(len(sys.argv)>2):
	exenames[0] = sys.argv[2]

exename = ""
for nameCheck in exenames:
	#print "exe entry '"+nameCheck+"' "  # debug
	if( os.path.isfile(nameCheck) ):
		exename = nameCheck

# write empty file if no exe found
if(exename == ""):
	writeHeader( outname, "\n// no executable found!\n\n" )
	print("Warning, no exe found - writing dummy header")
	exit(0); # dont throw an error for make, we can still continue...

if(doDebug):
	print("Params: outname '"+outname+"' , exename '"+exename)

# read old contents
oldContent = ""
doWrite    = True
try:
	infile = open(outname, "r")
	oldContent = infile.read()
	infile.close()
	if(doDebug):
		print("\n Old file content '"+oldContent+"' end \n")
except IOError:
	if(doDebug):
		print("Old file not found...")


# get hg version
#hgVersion = os.popen(exename+" id").read() 
# get gid id
hgVersion = os.popen(exename+" log -1 ").read() 
# remove newlines...
hgVersion = hgVersion.splitlines()[0]
hgVersion = hgVersion.rstrip('\n')
if(doDebug):
	print( "Got hg info: '" + hgVersion +"' " )

# matches old?
newContent = "\n\n#define MANTA_HG_VERSION \"" + hgVersion + "\" \n\n" 

if(newContent == oldContent):
	if(doDebug):
		print("MATCHES! No write")
	doWrite = False
else:
	if(doDebug):
		print("Old info different, writing")

# write temp file
if(doWrite):
	writeHeader( outname, newContent )
	print( "Updated hg info header , "+hgVersion )


