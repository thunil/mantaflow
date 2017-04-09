#!/bin/bash
# 
# note, uses & deletes test indices 7xx, and 8xx
# to really clean all tests, remove test_07?? and test_08?? dirs
# expects training data in sim 1000, and applies to 2007 and 2008 (shortened data sets)

BPATH=../data/
SAFEDEL='rm -r'
SCENE_TRAIN=1000
SCENE_NP=2008
SCENE_UNI=2007

echo Using data path ${BPATH}

# --- uni file tests ---

# optionally, clear all tiles
#find ${BPATH}/sim_${SCENE_TRAIN}/ -iname "tiles*" -exec rm -fr \{\} \;
#find ${BPATH}/sim_${SCENE_NP}/    -iname "tiles*" -exec rm -fr \{\} \;

# remove old output dir
echo
echo "************** Test 1 **************"
${SAFEDEL} ${BPATH}test_0700
${SAFEDEL} ${BPATH}test_0800

# train a model
python tf_train.py out 0 basePath ${BPATH}  useVelocities 0  trainingEpochs 1000   alwaysSave 1  testPathStartNo 700 fromSim ${SCENE_TRAIN}        simSizeLow 64

# and apply to small data set
python tf_train.py out 1 basePath ${BPATH}  useVelocities 0  testPathStartNo 800  fromSim ${SCENE_NP} toSim -1  loadModelTest 700 loadModelNo 0    simSizeLow 128


# same for a model using velocities
echo
echo "************** Test 2 **************"
${SAFEDEL} ${BPATH}test_0710
${SAFEDEL} ${BPATH}test_0810

python tf_train.py out 0 basePath ${BPATH}  useVelocities 1  trainingEpochs 1000  alwaysSave 1  testPathStartNo 710 fromSim ${SCENE_TRAIN}         simSizeLow 64 
python tf_train.py out 1 basePath ${BPATH}  useVelocities 1  testPathStartNo 810  fromSim ${SCENE_NP} toSim -1  loadModelTest 710 loadModelNo 0    simSizeLow 128 

# --- uni file tests ---

# clear all tiles
#find ${BPATH}/sim_${SCENE_TRAIN}/ -iname "tiles*" -exec rm -fr \{\} \;
#find ${BPATH}/sim_${SCENE_UNI}/   -iname "tiles*" -exec rm -fr \{\} \;

# train a model without vels
echo
echo "************** Test 3 **************"
${SAFEDEL} ${BPATH}test_0720
${SAFEDEL} ${BPATH}test_0820

python tf_train.py out 0 basePath ${BPATH}   useVelocities 0  trainingEpochs 1000  alwaysSave 1 testPathStartNo 720 fromSim ${SCENE_TRAIN}         simSizeLow 64  fileFormat uni
python tf_train.py out 1 basePath ${BPATH}   useVelocities 0  testPathStartNo 820  fromSim ${SCENE_UNI} toSim -1  loadModelTest 720 loadModelNo 0  simSizeLow 128 fileFormat uni

# train a model with vels
echo
echo "************** Test 4 **************"
${SAFEDEL} ${BPATH}test_0730
${SAFEDEL} ${BPATH}test_0830

python tf_train.py out 0 basePath ${BPATH}   useVelocities 1  trainingEpochs 1000  alwaysSave 1 testPathStartNo 730 fromSim ${SCENE_TRAIN}         simSizeLow 64  fileFormat uni
python tf_train.py out 1 basePath ${BPATH}   useVelocities 1  testPathStartNo 830  fromSim ${SCENE_UNI} toSim -1  loadModelTest 730 loadModelNo 0  simSizeLow 128 fileFormat uni


