#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2019 Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# numpty to uni, 2 command line parameters:
# 1) uni input
# 2) npz output
# (distinction between scalar vs vec3 content is determined by npz format, in contrast to uni2numpy.py)
# 
#******************************************************************************

import sys , os , re , shutil
from manta import *
import numpy as np

print("Converting npz file '"+sys.argv[1]+ "' to '" + sys.argv[2] + "' ") 

filename = sys.argv[1] 
if(os.path.isfile(sys.argv[2])):
	print("Error: target already exists")
	exit(1)

a=np.load(filename); 
a=a[a.files[-1]] # for npz files , note: not necessary for .npy
print("Read npz file with shape & stats: "+format([a.shape, np.mean(a), np.std(a), np.min(a),np.max(a)] ))
channels = a.shape[-1] # scalar

if (len(a.shape)==4):
	dim=3
	if channels!=1 and channels!=3:
		print("Error: invalid no of channels in numpy array (only 1 or 3 allowed)")
		exit(1)
elif (len(a.shape)==3):
	dim = 2
	if channels!=1 and channels!=2:
		print("Error: invalid no of channels in numpy array (only 1 or 2 allowed for 2D)")
		exit(1)
	if channels==2:
		# append third coord for 2d vectors
		print("Current shape & stats: "+format([a.shape, np.mean(a),np.std(a), np.min(a),np.max(a)] ))
		app = np.zeros([ a.shape[0], a.shape[1] , 1 ])
		a = np.concatenate( [a,app], axis=-1)
		print("New shape & stats: "+format([a.shape, np.mean(a),np.std(a), np.min(a),np.max(a)] ))
else:
	print("Error: invalid shape")
	exit(1)

gs = vec3(0,0,0)
if dim==3:
	gs = vec3( int(a.shape[2]), int(a.shape[1]), int(a.shape[0]) )  # mantaflow uses [x,y,z]
if dim==2:
	gs = vec3( int(a.shape[1]), int(a.shape[0]) , 1 )

print("Using grid size " + str(gs) + " , dim "+str(dim) )

# solver setup
s = Solver(name='main', gridSize = gs, dim=dim)
realGrid   = s.create(RealGrid)
vecGrid    = s.create(VecGrid)

if channels==1:
	copyArrayToGridReal(source=a, target=realGrid )
	realGrid.save(sys.argv[2])
else:
	copyArrayToGridVec3(source=a, target=vecGrid )
	vecGrid.save(sys.argv[2])

print("Wrote '" + sys.argv[2] +"' ")
