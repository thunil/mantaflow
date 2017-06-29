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

# path to sim data, trained models and output are also saved here
basePath = '../data/'

# main mode switch:
outputOnly = True  # apply model, or run full training?
outputInputs = False

simSizeLow   = 64
tileSizeLow  = 16
upRes        = 4

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
fromSim = toSim = -1 # range of sim data directories to use
keepAll         = False
numTests        = 10      # evaluate on 10 data points from test data
randSeed        = 1
brightenOutput  = -1 # multiplied with output to brighten it up
outputDataName  = '' # name of data to be regressed; by default, does nothing (density), e.g. if output data is pressure set to "pressure"
bWidth          = -1 # boundaryWidth to be cut away. 0 means 1 cell, 1 means two cells. In line with "bWidth" in manta scene files
# optional, add velocity as additional channels to input
useDensity      = 1  # default, only density
useVelocities   = 0
# load pressure as output
outputPressure  = 0

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
basePath        =     ph.getParam( "basePath",        basePath)
useDensity		= int(ph.getParam( "useDensity",      useDensity  ))
useVelocities   = int(ph.getParam( "useVelocities",   useVelocities  ))
outputPressure	= int(ph.getParam( "outputPressure",  outputPressure  ))
testPathStartNo = int(ph.getParam( "testPathStartNo", testPathStartNo  ))
fromSim         = int(ph.getParam( "fromSim",         fromSim  )) 
toSim           = int(ph.getParam( "toSim",           toSim  ))
alwaysSave      = int(ph.getParam( "alwaysSave",      False  )) # by default, only save checkpoint when cost is lower, can be turned off here
randSeed        = int(ph.getParam( "randSeed",        randSeed )) 
simSizeLow      = int(ph.getParam( "simSizeLow",      simSizeLow )) 
upRes           = int(ph.getParam( "upRes",           upRes ))
outputInputs    = int(ph.getParam( "outInputs",       outputInputs)) # for debugging, write images for input data
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
	if fromSim==-1:
		fromSim = toSim = 1000 # short, use single sim

	if cropOverlap>0:
		print("ERROR: dont use cropOverlap != 0 for training...")
		exit(1)

else:
	keepAll = True
	# dont train, just apply to input seq, by default use plume (2004)
	if fromSim==-1:
		fromSim = toSim = 3000

# check if batchsize is multiple of tilesInImage
tileSizeHiCrop = upRes * cropTileSizeLow
tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2
if not (batchSize % tilesPerImg == 0):
	print("WARNING: #batchSize (%d) is no multiple of #tilesPerImage (%d). " % (batchSize, tilesPerImg) +
	"For convolutional transpose networks, this might create problems while generating the output (out=1).")

# ---------------------------------------------

n_input  = tileSizeLow  ** 2 
n_output = tileSizeHigh ** 2
n_inputChannels = 0

if useDensity:
	n_inputChannels += 1
if useVelocities:
	n_inputChannels += 3
if n_inputChannels == 0:
	print("ERROR: No inputs set. Use useDensity, useVelocities, ...")
	exit()
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
			print('ERROR: Model with id below 200 does not exist. Please specify model id as "loadModelNo".')
			exit()
		# print('Latest model: %d.' % loadModelNo)

	load_path = basePath + 'test_%04d/model_%04d.ckpt' % (loadModelTest, loadModelNo)

test_path = next_test_path(testPathStartNo)

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
x      = tf.placeholder(tf.float32, shape=[None, 1, tileSizeLow, tileSizeLow, n_inputChannels])  # fixed to 2D for now
y_true = tf.placeholder(tf.float32, shape=[None, 1, tileSizeHigh,tileSizeHigh, 1             ])  # fixed to 2D for now , always output 1 channel
keep_prob = tf.placeholder(tf.float32)

layer_1_size = 512

xIn = tf.reshape(x, shape=[-1, n_input ]) # flatten
fc_1_weight = tf.Variable(tf.random_normal([n_input, layer_1_size], stddev=0.01))
fc_1_bias   = tf.Variable(tf.random_normal([layer_1_size], stddev=0.01))

fc1 = tf.add(tf.matmul(xIn, fc_1_weight), fc_1_bias)
fc1 = tf.nn.tanh(fc1)
fc1 = tf.nn.dropout(fc1, dropout)

fc_out_weight = tf.Variable(tf.random_normal([tileSizeLow*tileSizeLow * 2, tileSizeHigh * tileSizeHigh], stddev=0.01))
fc_out_bias   = tf.Variable(tf.random_normal([tileSizeHigh * tileSizeHigh], stddev=0.01))

y_pred = tf.add(tf.matmul(fc1, fc_out_weight), fc_out_bias)
y_pred = tf.reshape( y_pred, shape=[-1, 1, tileSizeHigh, tileSizeHigh, 1])


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

if outputPressure:
	outputDataName = "pressure"

# load test data
tiCr.loadTestDataNpz(fromSim, toSim, emptyTileValue, cropTileSizeLow, cropOverlap, 0.95, 0.05, load_vel=useVelocities, low_res_size=simSizeLow, upres=upRes, keepAll=keepAll, special_output_type=outputDataName, bWidth=bWidth)

if useVelocities and not useDensity:
	tiCr.reduceInputsToVelocity(dimensions=3)
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
				avgCostTest = 0.0
				for curr in range(numTests):
					batch_xs, batch_ys = tiCr.selectRandomTiles(batchSize, isTraining=False)
					cost_test, summary_test = sess.run([costFunc, lossTest], feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: 1.})
					avgCostTest += cost_test
				avgCostTest /= numTests

				avgCost /= testInterval
				print('\nEpoch {:04d}/{:04d} - Avg. cost= {:.9f} - Avg. validation cost= {:.9f}'.format((epoch + 1), trainingEpochs, avgCost, avgCostTest))
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
	print('To apply the trained model, set "outputOnly" to True, and insert numbers for "load_model_test", and "load_model_no" ')

else: 

	# ---------------------------------------------
	# outputOnly: apply to a full data set, and re-create full outputs from tiles

	batch_xs, batch_ys = tiCr.tile_inputs_all_complete, tiCr.tile_outputs_all_complete
	print('Creating %d tile outputs...' % (len(batch_xs)) )

	tileSizeHiCrop = upRes * cropTileSizeLow
	tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2
	outrange = int( len(tiCr.tile_inputs_all_complete) / tilesPerImg )

	# use int to avoid TypeError: 'float' object cannot be interpreted as an integer
	for currOut in range(outrange): 
		batch_xs = []
		batch_ys = []
		batch_velocity = [] # for optional output of velocity input pngs

		combine_tiles_amount = tilesPerImg
		# use batchSize if its a multiple of tilesPerImage. Important for conv_trans networks.
		stopOutput = False
		if batchSize % tilesPerImg == 0:
			combine_tiles_amount = batchSize
		for curr_tile in range(combine_tiles_amount):
			idx = currOut * tilesPerImg + curr_tile
			if combine_tiles_amount == batchSize and idx > len(tiCr.tile_inputs_all_complete) - 1:
				# Stop output, if there are not enough tiles left to evaluate a full batch
				# This will result in missing output pngs but prevents a crashes for conv_trans networks
				stopOutput = True
				break
			batch_xs.append(tiCr.tile_inputs_all_complete[idx])
			# batch_ys.append(np.zeros((tileSizeHigh * tileSizeHigh), dtype='f'))
			batch_ys.append(tiCr.tile_outputs_all_complete[idx])

			# to output velocity inputs
			if useVelocities and outputInputs:
				batch_velocity.append(tiCr.tile_inputs_all_complete[idx])

		if stopOutput:
			break
		resultTiles = y_pred.eval(feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: 1.})

		if brightenOutput > 0:
			for i in range(len(resultTiles)):
				resultTiles[i] *= brightenOutput
			for i in range(len(batch_ys)):
				batch_ys[i] *= brightenOutput
		tiCr.debugOutputPngsCrop(resultTiles, tileSizeHigh, simSizeHigh, test_path, imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg)
		tiCr.debugOutputPngsCrop(batch_ys,    tileSizeHigh, simSizeHigh, test_path, imageCounter=currOut, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg, name='expected_out')

		if outputInputs:
			if not useVelocities:
				tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='input')
			else: 
				tiCr.debugOutputPngsSingle(batch_velocity, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='in_vel_x', channel=1)
				tiCr.debugOutputPngsSingle(batch_velocity, tileSizeLow, simSizeLow, test_path, imageCounter=currOut, name='in_vel_y', channel=2)


# write summary to test overview
loaded_model = ''
if not loadModelTest == -1:
	loaded_model = ', Loaded %04d, %d' % (loadModelTest , loadModelNo)
with open(basePath + 'test_overview.txt', "a") as text_file:
	text_file.write(test_path[-10:-1] + ': {:.2f} min, {} Epochs, cost {:.4f}, {}'.format(training_duration, trainingEpochs, cost, " ") + loaded_model + '\n')


