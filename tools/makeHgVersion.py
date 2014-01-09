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
	print "Usage makeHgVersion.py <out-file> <optional: path-to-hg> "
	print "Warning, the target file <out-file> will be overwritten! "
	exit(1)

# target file
outname = sys.argv[1]

# path to hg executable, try a few options
hgnames = [ "/opt/local/bin/hg", "/opt/local/bin/hg", "/usr/local/bin/hg" ]
# note /opt/local/bin/hg is a double entry, can be overwritten by command line arg
hgname = ""
# optionally, make argument
if(len(sys.argv)>2):
	hgnames[0] = sys.argv[2]

for hgnameCheck in hgnames:
	if( os.path.isfile(hgnameCheck) ):
		hgname = hgnameCheck

# write empty file if no hg found
if(hgname == ""):
	writeHeader( outname, "\n// no hg found!\n\n" )
	print("Warning, no hg found - writing empty header")
	exit(0); # dont throw an error for make, we can still continue...

if(doDebug):
	print("Params: outname '"+outname+"' , hgname '"+hgname)

# read old contents
oldContent = ""
doWrite    = True
try:
	infile = open(outname, "r")
	oldContent = infile.read()
	infile.close()
	if(doDebug):
		print "\n Old file content '"+oldContent+"' end \n"
except IOError:
	if(doDebug):
		print "Old file not found..."


# get hg version
hgVersion = os.popen(hgname+" id").read() 
# remove newlines...
hgVersion = hgVersion.rstrip('\n')
if(doDebug):
	print( "Got hg info: '" + hgVersion +"' " )

# matches old?
newContent = "\n\n#define MANTA_HG_VERSION \"" + hgVersion + "\" \n\n" 

if(newContent == oldContent):
	if(doDebug):
		print "MATCHES! No write"
	doWrite = False
else:
	if(doDebug):
		print "Old info different, writing"

# write temp file
if(doWrite):
	writeHeader( outname, newContent )
	print( "Updated hg info header , "+hgVersion )


