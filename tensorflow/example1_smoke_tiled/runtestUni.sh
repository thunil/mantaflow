#!/bin/bash

# add path!
BPATH=/Users/name/data/
# optionally use sth. safer
SAFEDEL=rm

# specify pre-trained model to use for testing
MODEL=101
MODELNO=18

echo "Warning - Input your own paths and model numbers in script before using..."
exit 1

# train a model
python modelUni.py out 0 basePath ${BPATH} trainingEpochs 400

# train a model with vels
python modelUni.py out 0 basePath ${BPATH} trainingEpochs 400  useVelocities 1

# short application of old one, clear tiles
${SAFEDEL} ~/temp/flow_tf_tiled_data_/sim_2007/frame_0000/tiles_16x16
python modelUni.py out 1 basePath ${BPATH}  load_model_test ${MODEL} load_model_no ${MODELNO}

