#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2019 Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# Uni to numpy converter, 3 command line parameters:
# 1) uni input
# 2) npz output
# [3) optional: no of channels, default 1, 3 for vec3]
# 
#******************************************************************************

import sys , os , re , shutil
from manta import *
import numpy as np

# test npz files with:
# import sys,numpy as np; a=np.load(filename); a=a[a.files[-1]]; print(format([a.shape, np.mean(a),np.std(a), np.min(a),np.max(a)] ));

print("Converting uni file '"+sys.argv[1]+ "' to npz '" + sys.argv[2] + "' ") 

filename = sys.argv[1] 
if(os.path.isfile(sys.argv[2])):
	print("Error: target already exists")
	exit(1)

channels = 1 # scalar
if len(sys.argv)>3:
	channels = int(sys.argv[3])
	if channels!=1 and channels!=3:
		print("Error: invalid no of channels (only 1 or 3 allowed)")
		exit(1)

def tryToGetSize( filename ):
	rfile = filename 
	size = vec3(0,0,0)
	if(os.path.isfile(rfile)):
		size = getUniFileSize(rfile) 
		#print("Tried to read " + str(rfile) + " with " + format(size) )
	return size

def tryToLoad( grid, filename):
	rfile = filename 
	print("Trying to load " + rfile)
	if(os.path.isfile(rfile)):
		grid.load(rfile)
		#printUniFileInfoString(rfile) # more detailed build info

# ------------------------------------------------------------------------------------------
# setup

gs = vec3(0,0,0)
gs = tryToGetSize( filename )

if(gs.x==0):
	print("Couldnt get file size!")
	exit(1)
	
dim = 3
if (gs.z==1):
	dim=2
print("Using grid size " + str(gs) + " , dim "+str(dim) )


# solver setup
s = Solver(name='main', gridSize = gs, dim=dim)
realGrid   = s.create(RealGrid)
vecGrid    = s.create(VecGrid)


data = np.zeros([ int(gs.z), int(gs.y), int(gs.x), channels ])
if channels==1:
	tryToLoad( realGrid, filename )
	copyGridToArrayReal(source=realGrid, target=data )
else:
	tryToLoad( vecGrid, filename )
	copyGridToArrayVec3(source=vecGrid, target=data )
	if dim==2: 
		data = data[...,0:2] # remove z components of vec3s

if dim==2: 
	data = np.reshape(data, data.shape[1:] )

print("Copied to np array "+format(data.shape))

np.savez_compressed( sys.argv[2] , data )	
#np.save( sys.argv[2] , data ) # for testing, write uncompressed .npy files
print("Wrote '" + sys.argv[2] +"' ")
