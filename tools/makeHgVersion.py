#!/usr/bin/python
import os
import shutil
import sys
import re
doDebug = False

# params

if(len(sys.argv)<2):
	print "Usage makeHgVersion.py <out-file> <optional: path-to-hg> "
	print "Warning, the target file <out-file> will be overwritten! "
	exit(1)

# target file
outname = sys.argv[1]

# optionally, make argument
hgname = "/opt/local/bin/hg"
if(len(sys.argv)>2):
	hgname = sys.argv[2]

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
	try:
		outfile = open(outname, "w")
		outfile.write( newContent )
		outfile.close()
		if(doDebug):
			print( "Wrote '" + outname +"' " )
		print( "Updated hg info header , "+hgVersion )
	except IOError:
		print("Warning, unable to write file '"+outname+"' ")
		exit(1)


