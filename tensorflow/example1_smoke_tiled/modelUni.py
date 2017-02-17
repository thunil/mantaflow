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

simSizeLow  = 128
tileSizeLow = 16
upRes       = 4
simSizeHigh = simSizeLow * upRes
tileSizeHigh= tileSizeLow  * upRes

# dont use for training! for applying model, add overlap here if necessary (i.e., cropOverlap>0) 
# note:  cropTileSizeLow + (cropOverlap * 2) = tileSizeLow
cropOverlap     = 0
cropTileSizeLow = tileSizeLow - 2*cropOverlap

emptyTileValue  = 0.01
learning_rate   = 0.00005
trainingEpochs  = 100000 # by default, stop only manualy with ctrl-c...
dropout         = 0.9 # slight...
batch_size      = 100
testInterval    = 200
saveInterval    = 1000
fromSim = toSim = -1
randSeed        = 1

# optional, add velocity as additional channels to input?
useVelocities   = 0


# ---------------------------------------------

# load an existing model when load_ values > -1
# when training , manually abort when it's good enough
# then enter test_XXXX id and model checkpoint ID below to load

load_model_test = -1
load_model_no   = -1
testPathStartNo = 1

# command line params
outputOnly      = int(ph.getParam( "out",             outputOnly ))>0
trainingEpochs  = int(ph.getParam( "trainingEpochs",  trainingEpochs ))
load_model_test = int(ph.getParam( "load_model_test", load_model_test ))
load_model_no   = int(ph.getParam( "load_model_no",   load_model_no   ))
basePath        =     ph.getParam( "basePath",        basePath        )
useVelocities   = int(ph.getParam( "useVelocities",   useVelocities  ))
testPathStartNo = int(ph.getParam( "testPathStartNo", testPathStartNo  ))
fromSim         = int(ph.getParam( "fromSim",         fromSim  ))
toSim           = int(ph.getParam( "toSim",           toSim  ))
useLegacyNet    = int(ph.getParam( "useLegacyNet",    False ))>0
randSeed        = int(ph.getParam( "randSeed",        randSeed )) 
ph.checkUnusedParams()

# initialize

if toSim==-1:
	toSim = fromSim
tiCr.setBasePath(basePath)

np.random.seed(randSeed)
tf.set_random_seed(randSeed)

if not outputOnly:
	# run train!
	load_model_test = -1 
	simSizeLow = 64
	if fromSim==-1:
		fromSim = toSim   = 1000 # short, use single sim

	if cropOverlap>0:
		print("Error - dont use cropOverlap != 0 for training...")
		exit(1)

else:
	# dont train, just apply to input seq, by default use plume (2004)
	if fromSim==-1:
		fromSim = toSim = 2007

# ---------------------------------------------

#n_input = ((tileSizeLow + overlapping * 2) * 2) ** 2
n_input  = tileSizeLow  ** 2 
n_output = tileSizeHigh ** 2
n_inputChannels = 1

if useVelocities:
	n_inputChannels = 4
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

if not load_model_test == -1:
	if not os.path.exists(basePath + 'test_%04d/' % load_model_test):
		print('ERROR: Test to load does not exist.')
	load_path = basePath + 'test_%04d/model_%04d.ckpt' % (load_model_test, load_model_no)

test_path = next_test_path(testPathStartNo)

# custom Logger to write Log to file
class Logger(object):
	def __init__(self):
		self.terminal = sys.stdout
		self.log = open(test_path + "logfile.log", "a")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

sys.stdout = Logger()

# print Variables to log
def print_variables():
	print('\nUsing variables:')
	print('FromSim: {}'.format(fromSim))
	print('ToSim: {}'.format(toSim))
	print('simSizeLow: {}'.format(simSizeLow))
	print('tileSizeLow: {}'.format(tileSizeLow))
	print('cropOverlap: {}'.format(cropOverlap))
	print('cropTileSizeLow: {}'.format(cropTileSizeLow))
	print('upRes: {}'.format(upRes))
	print('emptyTileValue: {}'.format(emptyTileValue))
	print('learning_rate: {}'.format(learning_rate))
	print('trainingEpochs: {}'.format(trainingEpochs))
	print('dropout: {}'.format(dropout))
	print('batch_size: {}'.format(batch_size))
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
pool = 4
if not useLegacyNet:
	# new, w stride
	clFMs = 8 / n_inputChannels
	cae.convolutional_layer(clFMs, [3, 3], tf.nn.relu)
	cae.max_pool([pool,pool], [pool,pool])

	flat_size = cae.flatten()
	cae.fully_connected_layer(flat_size, tf.nn.relu)
	cae.unflatten()

	cae.max_depool([pool,pool], [pool,pool])
	cae.deconvolutional_layer(4, [3, 3], tf.nn.relu)

	cae.max_depool([pool,pool], [pool,pool])
	cae.deconvolutional_layer(2, [5, 5], tf.nn.relu)
else:
	# old org, without stride, just for loading 101:18, no vel
	print("Warning - use legacy network, todo remove...")
	cae.convolutional_layer(4, [3, 3], tf.nn.relu)
	cae.max_pool([pool,pool])

	flat_size = cae.flatten()
	cae.fully_connected_layer(flat_size, tf.nn.relu)
	cae.unflatten()

	cae.max_depool([pool,pool])
	cae.deconvolutional_layer(2, [3, 3], tf.nn.relu)

	cae.max_depool([2,2])
	cae.deconvolutional_layer(2, [3, 3], tf.nn.relu)

	cae.max_depool([2,2])
	cae.deconvolutional_layer(1, [3, 3], tf.nn.relu)


y_pred = tf.reshape( cae.y(), shape=[-1, (tileSizeHigh) *(tileSizeHigh)* 1])
print "DOFs: %d " % cae.getDOFs()

costFunc = tf.nn.l2_loss(y_true - y_pred) 
optimizer = tf.train.AdamOptimizer(learning_rate).minimize(costFunc)

# create session and saver
sess = tf.InteractiveSession()
saver = tf.train.Saver()

# init vars or load model
if load_model_test == -1:
	sess.run(tf.global_variables_initializer())
else:
	saver.restore(sess, load_path)
	print("Model restored from %s." % load_path)


# load test data
tiCr.loadTestDataUni(fromSim, toSim, emptyTileValue, cropTileSizeLow, cropOverlap, 0.95, 0.05, load_vel=useVelocities, low_res_size=simSizeLow, upres=upRes)

# create a summary to monitor cost tensor
lossTrain  = tf.summary.scalar("loss", costFunc)
lossTest   = tf.summary.scalar("loss_test", costFunc)
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
		error_per_display = 0
		startTime = time.time()
		epochTime = startTime
		last_save_since = saveInterval
		last_save_cost = 1e10
		for epoch in range(trainingEpochs):
			batch_xs, batch_ys = tiCr.selectRandomTiles(batch_size)
			_, cost, summary = sess.run([optimizer, costFunc, lossTrain], feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: dropout})

			# save model
			if (cost < last_save_cost) & (last_save_since > saveInterval):
				saver.save(sess, test_path + 'model_%04d.ckpt' % save_no)
				save_no += 1
				last_save_since = 0
				last_save_cost = cost
				print('Saved Model with cost %f.' % cost)
			else:
				last_save_since += 1

			# display error
			error_per_display += cost
			if (epoch + 1) % testInterval == 0:
				cumulated_cost_test = 0.0
				test_count = 10
				for curr_test in range(test_count):
					batch_xs, batch_ys = tiCr.selectRandomTiles(batch_size, isTraining=False)
					cost_test, summary_test = sess.run([costFunc, lossTest], feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: 1.})
					cumulated_cost_test += cost_test
				cumulated_cost_test /= test_count

				error_per_display /= testInterval
				print('\nEpoch {:04d} - Cost= {:.9f} - Cost_test= {:.9f}'.format((epoch + 1), error_per_display, cumulated_cost_test))
				print('%d epoches took %.02f seconds.' % (testInterval, (time.time() - epochTime)))
				#print('Estimated time: %.02f minutes.' % ((trainingEpochs - epoch) / testInterval * (time.time() - epochTime) / 60.0))
				epochTime = time.time()
				summary_writer.add_summary(summary, epoch)
				summary_writer.add_summary(summary_test, epoch)
				error_per_display = 0

	except KeyboardInterrupt:
		pass

	print('\n*****TRAINING FINISHED*****')
	training_duration = (time.time() - startTime) / 60.0
	print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, set "outputOnly" to True, and insert numbers for "load_model_test", and "load_model_no" ')

else:


	# ---------------------------------------------
	# Test against all data

	batch_xs, batch_ys = tiCr.tile_inputs_all_complete, tiCr.tile_outputs_all_complete
	tileSizeHiCrop = upRes * cropTileSizeLow
	tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2

	img_count = 0
	for currOut in range(len(tiCr.tile_inputs_all_complete) / tilesPerImg):
		batch_xs = []
		batch_ys = []
		for curr_tile in range(tilesPerImg):
			idx = currOut * tilesPerImg + curr_tile
			batch_xs.append(tiCr.tile_inputs_all_complete[idx])
			batch_ys.append(np.zeros((tileSizeHigh * tileSizeHigh), dtype='f'))

		resultTiles = y_pred.eval(feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: 1.})

		tiCr.debugOutputPngsCrop(resultTiles, tileSizeHigh, simSizeHigh, test_path, \
			imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg)
		img_count += 1

	print('Test finished, %d pngs written to %s.' % (img_count, test_path) )


# write summary to test overview
loaded_model = ''
if not load_model_test == -1:
	loaded_model = ', Loaded %04d, %d' % (load_model_test , load_model_no)
with open(basePath + 'test_overview.txt', "a") as text_file:
	text_file.write(test_path[-10:-1] + ': {:.2f} min, {} Epochs, cost {:.4f}, {}'.format(training_duration, trainingEpochs, cost, " ") + loaded_model + '\n')

