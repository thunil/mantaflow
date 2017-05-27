#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2017 Daniel Hook, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL) 
# http://www.gnu.org/licenses
#
# Manta & tensor flow example with tiles
# loads and writes uni files
#
#******************************************************************************

import time
import os
import shutil
import sys
import math

import tensorflow as tf
import numpy as np

# load manta tools
sys.path.append("../tools")
import tilecreator as tiCr
import uniio
import paramhelpers as ph
from convautoenc import ConvolutionalAutoEncoder

# path to sim data, trained models and output are also saved here
basePath = '../data/'

# main mode switch:
outputOnly = True  # apply model, or run full training?
outputInputs = False

simSizeLow   = 64
tileSizeLow  = 16
upRes        = 1

# dont use for training! for applying model, add overlap here if necessary (i.e., cropOverlap>0) 
# note:  cropTileSizeLow + (cropOverlap * 2) = tileSizeLow
cropOverlap     = 0
cropTileSizeLow = tileSizeLow - 2*cropOverlap

emptyTileValue  = 0.01
learningRate    = 0.00005
trainingEpochs  = 10000 # for large values, stop manualy with ctrl-c...
dropout         = 0.9   # slight...
batchSize       = 96
testInterval    = 200
saveInterval    = 1000
fromSim = toSim = -1      # range of sim data directories to use
keepAll         = False
numTests        = 10      # evaluate on 10 data points from test data
randSeed        = 1
fileFormat      = "npz"
brightenOutput  = -1   # multiplied with output to brighten it up
outputDataName  = ''   # name of data to be regressed; by default, does nothing (density), e.g. if output data is pressure set to "pressure"
bWidth          = -1   # boundaryWidth to be cut away. 0 means 1 cell, 1 means two cells. In line with "bWidth" in manta scene files

# optional, add velocity as additional channels to input?
useVelocities   = 0


# ---------------------------------------------

# load an existing model when load_ values > -1
# when training , manually abort when it's good enough
# then enter test_XXXX id and model checkpoint ID below to load

loadModelTest = -1
loadModelNo   = -1
testPathStartNo = 1

# command line params, explanations mostly above with variables
outputOnly      = int(ph.getParam( "out",             outputOnly ))>0
trainingEpochs  = int(ph.getParam( "trainingEpochs",  trainingEpochs ))
loadModelTest   = int(ph.getParam( "loadModelTest",   loadModelTest))
loadModelNo     = int(ph.getParam( "loadModelNo",     loadModelNo))
basePath        =     ph.getParam( "basePath",        basePath        )
useVelocities   = int(ph.getParam( "useVelocities",   useVelocities  ))
testPathStartNo = int(ph.getParam( "testPathStartNo", testPathStartNo  ))
fromSim         = int(ph.getParam( "fromSim",         fromSim  ))
toSim           = int(ph.getParam( "toSim",           toSim  ))
alwaysSave      = int(ph.getParam( "alwaysSave",      False  ))      # by default, only save checkpoint when cost is lower, can be turned off here
randSeed        = int(ph.getParam( "randSeed",        randSeed )) 
simSizeLow      = int(ph.getParam( "simSizeLow",      simSizeLow )) 
upRes           = int(ph.getParam( "upRes",           upRes ))
fileFormat      =     ph.getParam( "fileFormat",      fileFormat) # either npz or uni
outputInputs    = int(ph.getParam( "outInputs",       outputInputs)) # create pngs for inputs
brightenOutput  = int(ph.getParam( "brightenOutput",  brightenOutput))
outputDataName  =    (ph.getParam( "outName",         outputDataName))
bWidth			= int(ph.getParam( "bWidth",          bWidth))
ph.checkUnusedParams()

# initialize
simSizeHigh  = simSizeLow   * upRes
tileSizeHigh = tileSizeLow  * upRes
if outputOnly: # dont discard
	emptyTileValue = -1.

if toSim==-1:
	toSim = fromSim
tiCr.setBasePath(basePath)

np.random.seed(randSeed)
tf.set_random_seed(randSeed)

#tiCr.copySimData( fromSim, toSim ); exit(1);  # debug, copy sim data to different ID

if not outputOnly:
	# run train!
	loadModelTest = -1
	# simSizeLow = 64
	if fromSim==-1:
		fromSim = toSim   = 1000 # short, use single sim

	if cropOverlap>0:
		print("ERROR: dont use cropOverlap != 0 for training...")
		exit(1)

else:
	keepAll = True
	# dont train, just apply to input seq, by default use plume (2004)
	if fromSim==-1:
		fromSim = toSim = 2007

# ---------------------------------------------

n_input  = tileSizeLow  ** 2 
n_output = tileSizeHigh ** 2
n_inputChannels = 1

if useVelocities:
	n_inputChannels = 2
n_input *= n_inputChannels

# create output dir
def next_test_path(folder_no = 1):
	test_path_addition = 'test_%04d/' % folder_no
	while os.path.exists(basePath + test_path_addition):
		folder_no += 1
		test_path_addition = 'test_%04d/' % folder_no 
	test_path = basePath + test_path_addition
	print("Using test dir '%s'" % test_path)
	os.makedirs(test_path)
	return test_path

# create model loading path
if not loadModelTest == -1:
	if not os.path.exists(basePath + 'test_%04d/' % loadModelTest):
		print('ERROR: Test to load does not exist.')
	# search for newest model if no loadModelNo is given
	if loadModelNo == -1:
		for currModel in range(0, 200):
			if os.path.isfile(basePath + 'test_%04d/model_%04d.ckpt.index' % (loadModelTest, currModel)):
				loadModelNo = currModel
		if loadModelNo == -1:
			print('ERROR: No checkpoint found. Either wrong model directory, or no checkpoint with an id below 200. Please specify correct "loadModelTest" directory number, or a manual checkpoint number with "loadModelNo".')
			exit()
		# print('Latest model: %d.' % loadModelNo)

	load_path = basePath + 'test_%04d/model_%04d.ckpt' % (loadModelTest, loadModelNo)

test_path = next_test_path(testPathStartNo)
uniio.backupFile(__file__, test_path)

# custom Logger to write Log to file
class Logger(object):
	def __init__(self):
		self.terminal = sys.stdout
		self.log = open(test_path + "logfile.log", "a")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

	def flush(self): 
		# to avoid errormsg, " AttributeError: 'Logger' object has no attribute 'flush' "
		pass
sys.stdout = Logger()

# print Variables to log
def print_variables():
	print('\nUsing variables:')
	print('fromSim: {}'.format(fromSim))
	print('toSim: {}'.format(toSim))
	print('simSizeLow: {}'.format(simSizeLow))
	print('tileSizeLow: {}'.format(tileSizeLow))
	print('cropOverlap: {}'.format(cropOverlap))
	print('cropTileSizeLow: {}'.format(cropTileSizeLow))
	print('upRes: {}'.format(upRes))
	print('emptyTileValue: {}'.format(emptyTileValue))
	print('learningRate: {}'.format(learningRate))
	print('trainingEpochs: {}'.format(trainingEpochs))
	print('dropout: {}'.format(dropout))
	print('batchSize: {}'.format(batchSize))
	print('\n')

print_variables()

# ---------------------------------------------
# TENSORFLOW SETUP
x = tf.placeholder(tf.float32, shape=[None, n_input])
y_true = tf.placeholder(tf.float32, shape=[None, n_output])
keep_prob = tf.placeholder(tf.float32)

xIn = tf.reshape(x, shape=[-1, tileSizeLow, tileSizeLow, n_inputChannels]) 
cae = ConvolutionalAutoEncoder(xIn)

# --- main graph setup ---
# for conv_trans nets, the output tiles have to be created in the same batch size
is_convolution_transpose_network = False

################## BEGIN 1FC MODEL #################################
# layer_1_size = 512
#
# fc_1_weight = tf.Variable(tf.random_normal([tileSizeLow*tileSizeLow * 2, layer_1_size], stddev=0.01))
# fc_1_bias   = tf.Variable(tf.random_normal([layer_1_size], stddev=0.01))
#
# fc1 = tf.add(tf.matmul(x, fc_1_weight), fc_1_bias)
# fc1 = tf.nn.tanh(fc1)
# fc1 = tf.nn.dropout(fc1, dropout)
#
# fc_out_weight = tf.Variable(tf.random_normal([tileSizeLow*tileSizeLow * 2, tileSizeHigh * tileSizeHigh], stddev=0.01))
# fc_out_bias   = tf.Variable(tf.random_normal([tileSizeHigh * tileSizeHigh], stddev=0.01))
#
# y_pred = tf.add(tf.matmul(fc1, fc_out_weight), fc_out_bias)
################## END 1FC MODEL #################################

################## BEGIN CONV DECONV MODEL #################################
is_convolution_transpose_network = True

# Create some wrappers for simplicity
def conv2d(x, W, b, strides=1):
	# Conv2D wrapper, with bias and relu activation
	x = tf.nn.conv2d(x, W, strides=[1, strides, strides, 1], padding='SAME')
	x = tf.nn.bias_add(x, b)
	return tf.nn.tanh(x)

def conv2d_trans(x, W, b, output_shape, strides=1):
	x = tf.nn.conv2d_transpose(x, W, output_shape, [1, strides, strides, 1], padding='SAME', name=None)
	# x = tf.nn.bias_add(x, b)
	return tf.nn.tanh(x)


def maxpool2d(x, k=2):
	# MaxPool2D wrapper
	return tf.nn.max_pool(x, ksize=[1, k, k, 1], strides=[1, k, k, 1],
						  padding='SAME')


# Create model
def conv_net(x, weights, biases, dropout):
	# Reshape input picture
	x = tf.reshape(x, shape=[-1, 16, 16, 2])

	# Convolution Layer
	conv1 = conv2d(x, weights['wc1'], biases['bc1'])
	# Max Pooling (down-sampling)
	conv1 = maxpool2d(conv1, k=2)

	# Convolution Layer
	conv2 = conv2d(conv1, weights['wc2'], biases['bc2'])
	# Max Pooling (down-sampling)
	conv2 = maxpool2d(conv2, k=2)

	# Fully connected layer
	# Reshape conv2 output to fit fully connected layer input (flatten)
	fc1 = tf.reshape(conv2, [-1, weights['wd1'].get_shape().as_list()[0]])
	fc1 = tf.add(tf.matmul(fc1, weights['wd1']), biases['bd1'])
	fc1 = tf.nn.tanh(fc1)
	# Apply Dropout
	fc1 = tf.nn.dropout(fc1, dropout)
	fc1 = tf.reshape(fc1, shape=[-1, 4, 4, 8])

	# Deconv Layer
	deconv1 = conv2d_trans(fc1, weights['wc3'], biases['bc3'], [batchSize, 8, 8, 4], strides=2)
	deconv2 = conv2d_trans(deconv1, weights['wc4'], biases['bc4'], [batchSize, 16, 16, 1], strides=2)
	out = tf.reshape(deconv2, [-1, n_output])

	# Additional, optional fc layer
	# out = tf.add(tf.matmul(out, weights['out']), biases['out'])
	return out

# Store layers weight & bias
weights = {
	# 5x5 conv, 2 features in, 4 out
	'wc1': tf.Variable(tf.random_normal([5, 5, 2, 4], stddev=0.01)),
	# 3x3 conv, 4 features in, 8 out
	'wc2': tf.Variable(tf.random_normal([3, 3, 4, 8], stddev=0.01)),
	# fully connected, 4*4*8 nodes
	'wd1': tf.Variable(tf.random_normal([4*4*8, 4*4*8], stddev=0.01)),
	'wd2': tf.Variable(tf.random_normal([4*4*8, 4*4*8], stddev=0.01)),
	# 5x5 conv, 2 features in, 4 out
	'wc3': tf.Variable(tf.random_normal([5, 5, 4, 8], stddev=0.01)),
	# 3x3 conv, 4 features in, 8 out
	'wc4': tf.Variable(tf.random_normal([3, 3, 1, 4], stddev=0.01)),
	# output size 256 , ie 16x16
	'out': tf.Variable(tf.random_normal([256, n_output], stddev=0.01))
}

biases = {
	'bc1': tf.Variable(tf.zeros([4])),
	'bc2': tf.Variable(tf.zeros([8])),
	'bd1': tf.Variable(tf.random_normal([128], stddev=0.01)),
	'bc3': tf.Variable(tf.zeros([8])),
	'bc4': tf.Variable(tf.zeros([4])),
	'out': tf.Variable(tf.random_normal([n_output], stddev=0.01))
}

y_pred = conv_net(x, weights, biases, dropout)
################## END CONV DECONV MODEL #################################

################## BEGIN CONV MODEL #################################
#
# # Create some wrappers for simplicity
# def conv2d(x, W, b, strides=1):
#     # Conv2D wrapper, with bias and relu activation
#     x = tf.nn.conv2d(x, W, strides=[1, strides, strides, 1], padding='SAME')
#     x = tf.nn.bias_add(x, b)
#     return tf.nn.tanh(x)
#
#
# def maxpool2d(x, k=2):
#     # MaxPool2D wrapper
#     return tf.nn.max_pool(x, ksize=[1, k, k, 1], strides=[1, k, k, 1],
#                           padding='SAME')
#
#
# # Create model
# def conv_net(x, weights, biases, dropout):
#     # Reshape input picture
#     x = tf.reshape(x, shape=[-1, 16, 16, 2])
#
#     # Convolution Layer
#     conv1 = conv2d(x, weights['wc1'], biases['bc1'])
#     # Max Pooling (down-sampling)
#     conv1 = maxpool2d(conv1, k=2)
#
#     # Convolution Layer
#     conv2 = conv2d(conv1, weights['wc2'], biases['bc2'])
#     # Max Pooling (down-sampling)
#     conv2 = maxpool2d(conv2, k=2)
#
#     # Fully connected layer
#     # Reshape conv2 output to fit fully connected layer input (flatten)
#     fc1 = tf.reshape(conv2, [-1, weights['wd1'].get_shape().as_list()[0]])
#     fc1 = tf.add(tf.matmul(fc1, weights['wd1']), biases['bd1'])
#     fc1 = tf.nn.tanh(fc1)
#     # Apply Dropout
#     fc1 = tf.nn.dropout(fc1, dropout)
#
#     # Output
#     out = tf.add(tf.matmul(fc1, weights['out']), biases['out'])
#     return out
#
# # Store layers weight & bias
# weights = {
#     # 5x5 conv, 2 inputs, 4 outputs
#     'wc1': tf.Variable(tf.random_normal([5, 5, 2, 4], stddev=0.01)),
#     # 3x3 conv, 4 inputs, 8 outputs
#     'wc2': tf.Variable(tf.random_normal([3, 3, 4, 8], stddev=0.01)),
#     # fully connected, 4*4*8 inputs, 256 outputs
#     'wd1': tf.Variable(tf.random_normal([4*4*8, 256], stddev=0.01)),
#     # 256 inputs, 256 outputs
#     'out': tf.Variable(tf.random_normal([256, n_output], stddev=0.01))
# }
#
# biases = {
#     'bc1': tf.Variable(tf.zeros([4])),
#     'bc2': tf.Variable(tf.zeros([8])),
#     'bd1': tf.Variable(tf.random_normal([256], stddev=0.01)),
#     'out': tf.Variable(tf.random_normal([n_output], stddev=0.01))
# }
#
# y_pred = conv_net(x, weights, biases, dropout)
################## END CONV MODEL #################################

################## BEGIN CONV CAE MODEL #################################
# pool = 2
# clFMs = int(8 / n_inputChannels)
# cae.convolutional_layer(clFMs, [3, 3], tf.nn.tanh)
# cae.max_pool([pool,pool], [pool,pool])
#
# flat_size = cae.flatten()
# cae.fully_connected_layer(flat_size, tf.nn.tanh)
# cae.fully_connected_layer(flat_size, tf.nn.tanh)
# cae.unflatten()
#
# cae.max_depool([pool,pool], [pool,pool])
# cae.deconvolutional_layer(8, [3, 3], tf.nn.relu)
#
# # cae.max_depool([pool,pool], [pool,pool])
# # cae.deconvolutional_layer(2, [5, 5], tf.nn.relu)
#
# y_pred = tf.reshape( cae.y(), shape=[-1, (tileSizeHigh) *(tileSizeHigh)* 1])
# print ("DOFs: %d " % cae.getDOFs())
################## END CONV CAE MODEL #################################


costFunc = tf.nn.l2_loss(y_true - y_pred)
optimizer = tf.train.AdamOptimizer(learningRate).minimize(costFunc)

# create session and saver
sess = tf.InteractiveSession()
saver = tf.train.Saver()

# init vars or load model
if loadModelTest == -1:
	sess.run(tf.global_variables_initializer())
else:
	saver.restore(sess, load_path)
	print("Model restored from %s." % load_path)


# load test data
if (fileFormat == "npz"):
	tiCr.loadTestDataNpz(fromSim, toSim, emptyTileValue, cropTileSizeLow, cropOverlap, 0.95, 0.05, load_vel=useVelocities, low_res_size=simSizeLow, upres=upRes, keepAll=keepAll, special_output_type=outputDataName, bWidth=bWidth)
elif (fileFormat == "uni"):
	print("\n ERROR: only npz file format supported for pressure tests."); exit()
else:
	print("\n ERROR: Unknown file format \"" + fileFormat + "\". Use \"npz\" or \"uni\"."); exit()


print('Reducing data to 2D velocity...')
tiCr.reduceInputsTo2DVelocity()
# TODO: Remove this from final code
# print('Normalizing tile values...')
# tiCr.normalizeInputTestData()
tiCr.splitTileData(0.95, 0.05)

# create a summary to monitor cost tensor
lossTrain  = tf.summary.scalar("loss",     costFunc)
lossTest   = tf.summary.scalar("lossTest", costFunc)
merged_summary_op = tf.summary.merge_all() 
summary_writer    = tf.summary.FileWriter(test_path, sess.graph)

# ---------------------------------------------
# START TRAINING
training_duration = 0.0
cost = 0.0
save_no = 0

if not outputOnly:
	try:
		print('\n*****TRAINING STARTED*****\n')
		print('(stop with ctrl-c)')
		avgCost = 0
		startTime = time.time()
		epochTime = startTime
		lastSave = 1
		lastCost = 1e10
		for epoch in range(trainingEpochs):
			batch_xs, batch_ys = tiCr.selectRandomTiles(batchSize)
			_, cost, summary = sess.run([optimizer, costFunc, lossTrain], feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: dropout})

			# save model
			if ((cost < lastCost) or alwaysSave) and (lastSave >= saveInterval):
				saver.save(sess, test_path + 'model_%04d.ckpt' % save_no)
				save_no += 1
				lastSave = 1
				lastCost = cost
				print('Saved Model with cost %f.' % cost)
			else:
				lastSave += 1

			# display error
			avgCost += cost
			if (epoch + 1) % testInterval == 0:
				accumulatedCost = 0.0
				for curr in range(numTests):
					batch_xs, batch_ys = tiCr.selectRandomTiles(batchSize, isTraining=False)
					cost_test, summary_test = sess.run([costFunc, lossTest], feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: 1.})
					accumulatedCost += cost_test
				accumulatedCost /= numTests

				avgCost /= testInterval
				print('\nEpoch {:04d}/{:04d} - Cost= {:.9f} - Cost_test= {:.9f}'.format((epoch + 1), trainingEpochs, avgCost, accumulatedCost))
				print('%d epochs took %.02f seconds.' % (testInterval, (time.time() - epochTime)))
				#print('Estimated time: %.02f minutes.' % ((trainingEpochs - epoch) / testInterval * (time.time() - epochTime) / 60.0))
				epochTime = time.time()
				summary_writer.add_summary(summary, epoch)
				summary_writer.add_summary(summary_test, epoch)
				avgCost = 0

	except KeyboardInterrupt:
		pass

	print('\n*****TRAINING FINISHED*****')
	training_duration = (time.time() - startTime) / 60.0
	print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, set "outputOnly" to True, and insert numbers for "loadModelTest", and "loadModelNo" (optional). E.g. "out 1 loadModelTest 42".')

else: 

	# ---------------------------------------------
	# outputOnly: apply to a full data set, and re-create full outputs from tiles
	print('Creating outputs...')
	batch_xs, batch_ys = tiCr.tile_inputs_all_complete, tiCr.tile_outputs_all_complete
	tileSizeHiCrop = upRes * cropTileSizeLow
	tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2

	img_count = 0
	outrange = len(tiCr.tile_inputs_all_complete) / tilesPerImg

	# use int to avoid TypeError: 'float' object cannot be interpreted as an integer
	for currOut in range( int(outrange) ): 
		batch_xs = []
		batch_ys = []
		# to output velocity inputs
		batch_velocity_x = []
		batch_velocity_y = []

		combine_tiles_amount = tilesPerImg
		if is_convolution_transpose_network:
			combine_tiles_amount = batchSize
		for curr_tile in range(combine_tiles_amount):
		# for curr_tile in range(tilesPerImg):
			idx = currOut * tilesPerImg + curr_tile
			if is_convolution_transpose_network and idx > len(tiCr.tile_inputs_all_complete) - 1:
				exit()
			batch_xs.append(tiCr.tile_inputs_all_complete[idx])
			# batch_ys.append(np.zeros((tileSizeHigh * tileSizeHigh), dtype='f'))
			batch_ys.append(tiCr.tile_outputs_all_complete[idx])

			# to output velocity inputs
			curr_tile = np.reshape(tiCr.tile_inputs_all_complete[idx], (2, tileSizeLow*tileSizeLow), order='C')
			batch_velocity_x.append(curr_tile[0])
			batch_velocity_y.append(curr_tile[1])

		resultTiles = y_pred.eval(feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: 1.})

		if brightenOutput > 0:
			for curr_value in range(len(resultTiles)):
				resultTiles[curr_value] *= brightenOutput
			for curr_value in range(len(batch_ys)):
				batch_ys[curr_value] *= brightenOutput
		tiCr.debugOutputPngsCrop(resultTiles, tileSizeHigh, simSizeHigh, test_path, imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg)
		tiCr.debugOutputPngsCrop(batch_ys, tileSizeHigh, simSizeHigh, test_path, imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg, name='expected_out')
		# tiCr.debugOutputPngsSingle(batch_ys, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='expected_out')

		if outputInputs:
			if not useVelocities:
				tiCr.debugOutputPngsSingle(batch_xs, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='input')
			else:
				tiCr.debugOutputPngsSingle(batch_velocity_x, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='vel_x')
				tiCr.debugOutputPngsSingle(batch_velocity_y, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='vel_y')

		# TODO: Remove this from final code
		# Debug: Outputs uni files that can be loaded into mantaflow (currently buggy)
		# tiCr.debugOutputPressureVelocityUni(batch_xs, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='pressure')

		# optionally, output references
		#tiCr.debugOutputPngsCrop(batch_ys, tileSizeHigh, simSizeHigh, test_path+"_ref", imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg)
		img_count += 1

	print('Test finished, %d pngs written to %s.' % (img_count, test_path) )

# write summary to test overview
loaded_model = ''
if not loadModelTest == -1:
	loaded_model = ', Loaded %04d, %d' % (loadModelTest , loadModelNo)
with open(basePath + 'test_overview.txt', "a") as text_file:
	text_file.write(test_path[-10:-1] + ': {:.2f} min, {} Epochs, cost {:.4f}, {}'.format(training_duration, trainingEpochs, cost, " ") + loaded_model + '\n')


