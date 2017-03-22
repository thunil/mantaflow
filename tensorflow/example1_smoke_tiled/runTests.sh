#!/bin/bash
# 
# note, uses & deletes test indices 7xx, and 8xx

BPATH=../data/
SAFEDEL=rm
SCENE1=2008
SCENE2=2007
echo Using data path ${BPATH}

# --- uni file tests ---

# optionally, clear all tiles
#find ${BPATH}/sim_1000/ -iname "tiles*" -exec rm -fr \{\} \;
#find ${BPATH}/sim_2008/ -iname "tiles*" -exec rm -fr \{\} \;

# remove old output dir
${SAFEDEL} ${BPATH}/test_0700
${SAFEDEL} ${BPATH}/test_0800

# train a model
python tf_model_np.py out 0 basePath ${BPATH}  useVelocities 0  trainingEpochs 1000   alwaysSave 1  testPathStartNo 700 fromSim ${SCENE1}

# and apply to small data set
python tf_model_np.py out 1 basePath ${BPATH}  useVelocities 0  testPathStartNo 800  fromSim ${SCENE1} toSim -1  load_model_test 700 load_model_no 0

# same for a model using velocities

${SAFEDEL} ${BPATH}/test_0710
${SAFEDEL} ${BPATH}/test_0810

python tf_model_np.py out 0 basePath ${BPATH}  useVelocities 1  trainingEpochs 1000  alwaysSave 1  testPathStartNo 710 fromSim ${SCENE1}

python tf_model_np.py out 1 basePath ${BPATH}  useVelocities 1  testPathStartNo 810  fromSim ${SCENE1} toSim -1  load_model_test 710 load_model_no 0

# --- uni file tests ---

# clear all tiles
#find ${BPATH}/sim_1000/ -iname "tiles*" -exec rm -fr \{\} \;
#find ${BPATH}/sim_2007/ -iname "tiles*" -exec rm -fr \{\} \;

# train a model without vels

${SAFEDEL} ${BPATH}/test_0720
${SAFEDEL} ${BPATH}/test_0820

python tf_model_uni.py out 0 basePath ${BPATH}   useVelocities 0  trainingEpochs 1000  alwaysSave 1 testPathStartNo 720 fromSim ${SCENE2}
python tf_model_uni.py out 1 basePath ${BPATH}   useVelocities 0  testPathStartNo 820  fromSim ${SCENE2} toSim -1  load_model_test 720 load_model_no 0

# train a model with vels

${SAFEDEL} ${BPATH}/test_0730
${SAFEDEL} ${BPATH}/test_0830

python tf_model_uni.py out 0 basePath ${BPATH}   useVelocities 1  trainingEpochs 1000  alwaysSave 1 testPathStartNo 730 fromSim ${SCENE2}
python tf_model_uni.py out 1 basePath ${BPATH}   useVelocities 1  testPathStartNo 830  fromSim ${SCENE2} toSim -1  load_model_test 730 load_model_no 0


