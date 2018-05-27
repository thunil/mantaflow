#******************************************************************************
#
# simplified L2 conv net training examples
# Copyright 2018 Nils Thuerey, You Xie, Erik Franz, Mengyu Chu
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
#******************************************************************************

import time
import os
import shutil
import sys
import math

import tensorflow as tf
from tensorflow.python.client import timeline
import numpy as np

# load manta tools
sys.path.append("../tools")
import tilecreator as tc
import uniio
import paramhelpers as ph
from GAN import GAN, lrelu
import fluiddataloader as FDL


# ---------------------------------------------

# initialize parameters / command line params
outputOnly	  = int(ph.getParam( "out",			 False ))>0 		# output/generation mode, main mode switch

basePath		=	 ph.getParam( "basePath",		'../data/' )
randSeed		= int(ph.getParam( "randSeed",		1 )) 				# seed for np and tf initialization
load_model_test = int(ph.getParam( "load_model_test", -1 )) 			# the number of the test to load a model from. can be used in training and output mode. -1 to not load a model
load_model_no   = int(ph.getParam( "load_model_no",   -1 )) 			# nubmber of the model to load

simSizeLow  	= int(ph.getParam( "simSize", 		  64 )) 			# tiles of low res sim
tileSizeLow 	= int(ph.getParam( "tileSize", 		  16 )) 			# size of low res tiles
upRes	  		= int(ph.getParam( "upRes", 		  4 )) 				# scaling factor

#Data and Output
loadPath		 =	 ph.getParam( "loadPath",		 '../data/' ) 	# path to training data
fromSim		   = int(ph.getParam( "fromSim",		 1000 )) 			# range of sim data to use, start index
toSim		   = int(ph.getParam( "toSim",		   -1   )) 			# end index
dataDimension   = int(ph.getParam( "dataDim",		 2 )) 				# dimension of dataset, can be 2 or 3. in case of 3D any scaling will only be applied to H and W (DHW)
numOut			= int(ph.getParam( "numOut",		  200 )) 			# number ouf images to output (from start of sim)
saveOut	  	    = int(ph.getParam( "saveOut",		 False ))>0 		# save output of output mode as .npz in addition to images
loadOut			= int(ph.getParam( "loadOut",		 -1 )) 			# load output from npz to use in output mode instead of tiles. number or output dir, -1 for not use output data
outputImages	=int(ph.getParam( "img",  			  True ))>0			# output images
#models
genModel		=	 ph.getParam( "genModel",		 'gen_test' ) 	# choose generator model
#Training
learning_rate   = float(ph.getParam( "learningRate",  0.0002 ))
decayLR		    = int(ph.getParam( "decayLR",			 False ))>0 		# decay learning rate?
dropout   		= float(ph.getParam( "dropout",  	  1.0 )) 			# keep prop for all dropout layers during training
dropoutOutput   = float(ph.getParam( "dropoutOutput", dropout )) 		# affects testing, full sim output and progressive output during training
beta			= float(ph.getParam( "adam_beta1",	 0.5 ))			#1. momentum of adam optimizer

k				= float(ph.getParam( "lambda",		  1.0)) 			# influence/weight of l1 term on generator loss
batch_size	    = int(ph.getParam( "batchSize",  	  128 ))			# batch size for pretrainig and output, default for batchSizeDisc and batchSizeGen
trainingEpochs  = int(ph.getParam( "trainingEpochs",  10000 )) 		# for GAN training
batch_norm		= int(ph.getParam( "batchNorm",	   True ))>0			# apply batch normalization to conv and deconv layers
use_spatialdisc = int(ph.getParam( "use_spatialdisc",		   True )) #use spatial discriminator or not

useVelocities   = int(ph.getParam( "useVelocities",   0  )) 			# use velocities or not

useDataAugmentation = int(ph.getParam( "dataAugmentation", 0 ))		 # use dataAugmentation or not
minScale = float(ph.getParam( "minScale",	  0.85 ))				 # augmentation params...
maxScale = float(ph.getParam( "maxScale",	  1.15 ))
rot	     = int(ph.getParam( "rot",		  2	 ))		#rot: 1: 90 degree rotations; 2: full rotation; else: nop rotation 
flip	 =   int(ph.getParam( "flip",		  1	 ))

#Test and Save
testPathStartNo = int(ph.getParam( "testPathStartNo", 0  ))
testInterval	= int(ph.getParam( "testInterval", 	  20  )) 			# interval in epochs to run tests should be lower or equal outputInterval
numTests		= int(ph.getParam( "numTests", 		  10  )) 			# number of tests to run from test data each test interval, run as batch
outputInterval	= int(ph.getParam( "outputInterval",  100  ))			# interval in epochs to output statistics
saveInterval	= int(ph.getParam( "saveInterval",	  200  ))	 		# interval in epochs to save model
alwaysSave	    = int(ph.getParam( "alwaysSave",	  True  )) 			#
maxToKeep		= int(ph.getParam( "keepMax",		 3  )) 			# maximum number of model saves to keep in each test-run
genTestImg		= int(ph.getParam( "genTestImg",	  -1 )) 			# if > -1 generate test image every output interval
note			= ph.getParam( "note",		   "" )					# optional info about the current test run, printed in log and overview
data_fraction	= float(ph.getParam( "data_fraction",		   0.3 ))
frame_max		= int(ph.getParam( "frame_max",		   200 ))
frame_min		= int(ph.getParam( "frame_min",		   0 ))

ph.checkUnusedParams()

# initialize
simSizeHigh 	= simSizeLow * upRes
tileSizeHigh	= tileSizeLow  * upRes

if not (dataDimension == 2 or dataDimension == 3):
	print('Unsupported data dimension {}. Only 2 and 3 are supported'.format(dataDimension))
	exit(1)

if toSim==-1:
	toSim = fromSim

channelLayout_low = 'd'
lowfilename = "density_low_%04d.uni"
highfilename = "density_high_%04d.uni"
lowfilename = "density_low_%04d.npz"
highfilename = "density_high_%04d.npz" # NT_DEBUG
mfl = ["density"]
mfh = ["density"]
if outputOnly: 
	highfilename = None
	mfh = None
if useVelocities:
	channelLayout_low += ',vx,vy,vz'
	mfl= np.append(mfl, "velocity")

dirIDs = np.linspace(fromSim, toSim, (toSim-fromSim+1),dtype='int16')

if (outputOnly): 
	data_fraction = 1.0
	useTempoD = False
	useTempoL2 = False
	useDataAugmentation = 0

#if ((not useTempoD) and (not useTempoL2)): # should use the full sequence, not use multi_files
tiCr = tc.TileCreator(tileSizeLow=tileSizeLow, simSizeLow=simSizeLow , dim =dataDimension, dim_t = 1, channelLayout_low = channelLayout_low, upres=upRes)
floader = FDL.FluidDataLoader( print_info=1, base_path=loadPath, filename=lowfilename, oldNamingScheme=False, filename_y=highfilename, filename_index_min=frame_min, filename_index_max=frame_max, indices=dirIDs, data_fraction=data_fraction, multi_file_list=mfl, multi_file_list_y=mfh)

if useDataAugmentation:
	tiCr.initDataAugmentation(rot=rot, minScale=minScale, maxScale=maxScale ,flip=flip)
inputx, y, xFilenames  = floader.get()
if (not outputOnly): 
	tiCr.addData(inputx,y)
elif dataDimension == 3:
	simLowLength = inputx.shape[1]
	simLowWidth = inputx.shape[2]
	simLowHeight = inputx.shape[3]

print("Random seed: {}".format(randSeed))
np.random.seed(randSeed)
tf.set_random_seed(randSeed)

# ---------------------------------------------

# 2D: tileSize x tileSize tiles; 3D: tileSize x tileSize x tileSize chunks
n_input  = tileSizeLow  ** 2
n_output = tileSizeHigh ** 2
if dataDimension == 3:
	n_input  *= tileSizeLow
	n_output *= (tileSizeLow*upRes)
n_inputChannels = 1

if useVelocities:
	n_inputChannels += 3
n_input *= n_inputChannels

# init paths
if not load_model_test == -1:
	if not os.path.exists(basePath + 'test_%04d/' % load_model_test):
		print('ERROR: Test to load does not exist.')
	load_path = basePath + 'test_%04d/model_%04d.ckpt' % (load_model_test, load_model_no)
	if outputOnly:
		out_path_prefix = 'out_%04d-%04d' % (load_model_test,load_model_no)
		test_path,_ = ph.getNextGenericPath(out_path_prefix, 0, basePath + 'test_%04d/' % load_model_test)

	else:
		test_path,_ = ph.getNextTestPath(testPathStartNo, basePath)

else:
	test_path,load_model_test_new = ph.getNextTestPath(testPathStartNo, basePath)

# logging & info
sys.stdout = ph.Logger(test_path)
print('Note: {}'.format(note))
print("\nCalled on machine '"+ os.uname()[1] +"' with: " + str(" ".join(sys.argv) ) )
print("\nUsing parameters:\n"+ph.paramsToString())
ph.writeParams(test_path+"params.json") # export parameters in human readable format

if outputOnly:
	print('*****OUTPUT ONLY*****')

if 0 and not outputOnly:
	os.makedirs(test_path+"/zbu_src")
	uniio.backupFile(__file__, test_path+"/zbu_src/")
	uniio.backupFile("../tools/tilecreator.py", test_path+"/zbu_src/")
	uniio.backupFile("../tools/GAN.py", test_path+"/zbu_src/") 
	uniio.backupFile("../tools/fluiddataloader.py", test_path+"/zbu_src/")

# ---------------------------------------------
# TENSORFLOW SETUP

import scipy.misc

def save_img(out_path, img):
	img = np.clip(img * 255.0, 0, 255).astype(np.uint8)
	scipy.misc.imsave(out_path, img)

def save_img_3d(out_path, img): 
	data = np.concatenate([np.sum(img, axis=0), np.sum(img, axis=1), np.sum(img, axis=2)], axis=0)
	save_img(out_path, data)
	


#input for gen
x = tf.placeholder(tf.float32, shape=[None, n_input])
#reference input for disc
x_disc = tf.placeholder(tf.float32, shape=[None, n_input])
#real input for disc
y = tf.placeholder(tf.float32, shape=[None, n_output])
kk = tf.placeholder(tf.float32)
#keep probablity for dropout
keep_prob = tf.placeholder(tf.float32)

print("x: {}".format(x.get_shape()))
# --- main graph setup ---

rbId = 0
def resBlock(gan, inp, s1,s2, reuse, use_batch_norm, filter_size=3):
	global rbId
	# note - leaky relu (lrelu) not too useful here

	# convolutions of resnet block
	if dataDimension == 2:
		filter = [filter_size,filter_size]
		filter1 = [1,1]
	elif dataDimension == 3:
		filter = [filter_size,filter_size,filter_size]
		filter1 = [1,1,1]

	gc1,_ = gan.convolutional_layer(  s1, filter, tf.nn.relu, stride=[1], name="g_cA%d"%rbId, in_layer=inp, reuse=reuse, batch_norm=use_batch_norm, train=train) #->16,64
	#gc1,_ = gan.convolutional_layer(  s1, filter, lrelu, stride=[1], name="g_cA%d"%rbId, in_layer=inp, reuse=reuse, batch_norm=use_batch_norm, train=train) #->16,64
	gc2,_ = gan.convolutional_layer(  s2, filter, None      , stride=[1], name="g_cB%d"%rbId,               reuse=reuse, batch_norm=use_batch_norm, train=train) #->8,128

	# shortcut connection
	gs1,_ = gan.convolutional_layer(s2, filter1 , None       , stride=[1], name="g_s%d"%rbId, in_layer=inp, reuse=reuse, batch_norm=use_batch_norm, train=train) #->16,64
	resUnit1 = tf.nn.relu( tf.add( gc2, gs1 )  )
	#resUnit1 = lrelu( tf.add( gc2, gs1 )  )
	rbId += 1
	return resUnit1

def gen_resnet(_in, reuse=False, use_batch_norm=False, train=None):
	global rbId
	print("\n\tGenerator (resize-resnett3-deep)")
	with tf.variable_scope("generator", reuse=reuse) as scope:

		if dataDimension == 2:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, n_inputChannels]) #NHWC
			patchShape = [2,2]
		elif dataDimension == 3:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, n_inputChannels]) #NDHWC
			patchShape = [2,2,2]
		rbId = 0
		gan = GAN(_in)
	
		gan.max_depool()
		inp = gan.max_depool()
		ru1 = resBlock(gan, inp, n_inputChannels*2,n_inputChannels*8,  reuse, use_batch_norm,5)# with tf.device('/device:GPU:1'):  after here, gpu1 crash

		ru2 = resBlock(gan, ru1, 128, 128,  reuse, use_batch_norm,5)# with tf.device('/device:GPU:1'):  after here, gpu0 crash
		inRu3 = ru2
		ru3 = resBlock(gan, inRu3, 32, 8,  reuse, use_batch_norm,5)
		ru4 = resBlock(gan, ru3, 2, 1,  reuse, False,5)
		resF = tf.reshape( ru4, shape=[-1, n_output] )
		print("\tDOFs: %d , %f m " % ( gan.getDOFs() , gan.getDOFs()/1000000.) ) 
		return resF


def gen_resnetSm(_in, reuse=False, use_batch_norm=False, train=None):
	global rbId
	print("\n\tGenerator (resize-resnett3-deep)")
	with tf.variable_scope("generator", reuse=reuse) as scope:

		if dataDimension == 2:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, n_inputChannels]) #NHWC
			patchShape = [2,2]
		elif dataDimension == 3:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, n_inputChannels]) #NDHWC
			patchShape = [2,2,2]
		rbId = 0
		gan = GAN(_in)
		gan.max_depool()
		inp = gan.max_depool()
		ru1 = resBlock(gan, inp, n_inputChannels*2,n_inputChannels*8,  reuse, use_batch_norm,3)
		ru2 = resBlock(gan, ru1, 16, 16,  reuse, use_batch_norm,3)
		inRu3 = ru2
		ru3 = resBlock(gan, inRu3, 8, 4,  reuse, use_batch_norm,3)
		ru4 = resBlock(gan, ru3, 2, 1,  reuse, False,5)
		resF = tf.reshape( ru4, shape=[-1, n_output] )
		print("\tDOFs: %d , %f m " % ( gan.getDOFs() , gan.getDOFs()/1000000.) ) 
		return resF



############################################gen_test###############################################################
def gen_test(_in, reuse=False, use_batch_norm=False, train=None):
	global rbId
	print("\n\tGenerator-test")
	with tf.variable_scope("generator-test", reuse=reuse) as scope:
		if dataDimension == 2:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, n_inputChannels]) #NHWC
			patchShape = [2,2]
		elif dataDimension == 3:
			_in = tf.reshape(_in, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, n_inputChannels]) #NDHWC
			patchShape = [2,2,2]
		rbId = 0
		gan = GAN(_in)

		gan.max_depool()
		i2np,_ = gan.deconvolutional_layer(32, patchShape, None, stride=[1,1], name="g_D1", reuse=reuse, batch_norm=False, train=train, init_mean=0.99) #, strideOverride=[1,1] )
		gan.max_depool()
		inp,_  = gan.deconvolutional_layer(1                   , patchShape, None, stride=[1,1], name="g_D2", reuse=reuse, batch_norm=False, train=train, init_mean=0.99) #, strideOverride=[1,1] )
		return 	tf.reshape( inp, shape=[-1, n_output] )


#change used models for gen and disc here #other models in NNmodels.py
gen_model = locals()[genModel]
#training or testing for batch norm
train = tf.placeholder(tf.bool)

if not outputOnly: #setup for training
	gen_part = gen_model(x, use_batch_norm=batch_norm, train=train)
	if genTestImg > -1: sampler = gen_part
else: #setup for generating output with trained model
	sampler = gen_model(x, use_batch_norm=batch_norm, train=False)

sys.stdout.flush()

if not outputOnly:

	#additional generator losses
	gen_l2_loss = tf.nn.l2_loss(y - gen_part)
	gen_l1_loss = tf.reduce_mean(tf.abs(y - gen_part)) #use mean to normalize w.r.t. output dims. tf.reduce_sum(tf.abs(y - gen_part))

	#uses sigmoid cross entropy and l1 - see cGAN paper
	gen_loss_complete = gen_l1_loss*kk 

	# set up decaying learning rate, if enabled
	lr_global_step = tf.Variable(0, trainable=False)
	learning_rate_scalar = learning_rate
	if decayLR:
		learning_rate = tf.train.polynomial_decay(learning_rate, lr_global_step, trainingEpochs//2, learning_rate_scalar*0.05, power=1.1)

	update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
	gen_update_ops = update_ops[:]
	ori_gen_update_ops = update_ops[:]
	pre_update_ops = update_ops[:]

	#variables to be used in the different otimization steps
	vars = tf.trainable_variables()
	g_var = [var for var in vars if "g_" in var.name]
	if use_spatialdisc:
		dis_update_ops = update_ops[:]
		d_var = [var for var in vars if "d_" in var.name]

	with tf.control_dependencies(pre_update_ops): 
		pretrain_optimizer = tf.train.AdamOptimizer(learning_rate).minimize(gen_l2_loss, var_list=g_var)


# create session and saver
config = tf.ConfigProto(allow_soft_placement=True)
sess = tf.InteractiveSession(config = config)
saver = tf.train.Saver(max_to_keep=maxToKeep)


# init vars or load model
if load_model_test == -1:
	sess.run(tf.global_variables_initializer())
else:
	saver.restore(sess, load_path)
	print("Model restored from %s." % load_path)



if not outputOnly:
	lossPretrain_gen  = tf.summary.scalar("generator_L2_loss",     gen_l2_loss)	
	merged_summary_op = tf.summary.merge_all()
	summary_writer    = tf.summary.FileWriter(test_path, sess.graph)

save_no = 0
tileSizeHi = upRes * tileSizeLow
if dataDimension == 2:
	tilesPerImg = (simSizeHigh // tileSizeHi) ** 2
else:
	tilesPerImg = (simSizeHigh // tileSizeHi) ** 3
image_no = 0
if not outputOnly:
	os.makedirs(test_path+'test_img/')
	#if pretrain>0: # or pretrain_gen > 0 or pretrain_disc>0:
	#	os.makedirs(test_path+'pretrain_test_img/') # NT_DEBUG ?

def modifyVel(Dens,Vel):
	return velout # not active right now...

def getInput(index = 1, randomtile = True, isTraining = True, batch_size = 1, useDataAugmentation = False):
	if randomtile == False:
		batch_xs, batch_ys = tiCr.getFrameTiles(index) 
	else:
		batch_xs, batch_ys = tiCr.selectRandomTiles(selectionSize = batch_size, augment=useDataAugmentation)	

	batch_xs = np.reshape(batch_xs, (-1, n_input))
	batch_ys = np.reshape(batch_ys, (-1, n_output))
	return batch_xs, batch_ys

#evaluate the generator (sampler) on the first step of the first simulation and output result
def generateTestImage(sim_no = fromSim, frame_no = 1, outPath = test_path,imageindex = 0):
	if (not outputOnly):
		batch_xs, _ = getInput(randomtile = False, index = (sim_no-fromSim)*frame_max + frame_no)
	else:
		batch_xs = inputx[frame_no]
	resultTiles = []
	for tileno in range(batch_xs.shape[0]):
		batch_xs_in = np.reshape(batch_xs[tileno],[-1, n_input])
		results = sess.run(sampler, feed_dict={x: batch_xs_in, keep_prob: dropoutOutput, train: False})
		resultTiles.extend(results)
	resultTiles = np.array(resultTiles)
	if dataDimension == 2: # resultTiles may have a different size
		imgSz = int(resultTiles.shape[1]**(1.0/2) + 0.5)
		resultTiles = np.reshape(resultTiles,[resultTiles.shape[0],imgSz,imgSz, 1])
	else:
		imgSz = int(resultTiles.shape[1]**(1.0/3) + 0.5)
		resultTiles = np.reshape(resultTiles,[resultTiles.shape[0],imgSz,imgSz,imgSz])
	tiles_in_image=[int(simSizeHigh/tileSizeHigh),int(simSizeHigh/tileSizeHigh)]
	tc.savePngsGrayscale(resultTiles,outPath, imageCounter=imageindex, tiles_in_image=tiles_in_image)

def saveModel(cost, exampleOut=-1, imgPath = test_path):
	global save_no
	saver.save(sess, test_path + 'model_%04d.ckpt' % save_no)
	msg = 'Saved Model %04d with cost %f.' % (save_no, cost)
	if exampleOut > -1:
		generateTestImage(imageindex = save_no, outPath = imgPath)
	save_no += 1
	return msg

# write summary to test overview
loaded_model = ''
if not load_model_test == -1:
	loaded_model = ', Loaded %04d, %04d' % (load_model_test , load_model_no)
with open(basePath + 'test_overview.log', "a") as text_file:
	if not outputOnly:
		text_file.write(test_path[-10:-1] + ': {}D, \"{}\"\n'.format(dataDimension, note))
		text_file.write('\t{} Epochs, gen: {} '.format(trainingEpochs, gen_model.__name__) + loaded_model + '\n')
		text_file.write('\tlambda: {}, dropout: {:.4f}({:.4f})'.format(k, dropout, dropoutOutput) + '\n')
	else:
		text_file.write('Output:' + loaded_model + ' (' + test_path[-28:-1] + ')\n')
		text_file.write('\ttile size: {}, seed: {}, dropout-out: {:.4f}'.format(tileSizeLow, randSeed, dropoutOutput) + '\n')

	
	
# ---------------------------------------------
# ---------------------------------------------
training_duration = 0.0

#train generator using L2 loss
if (not outputOnly): # and pretrain>0:
	try:
		print('Generator using L2' + '{} epochs\n'.format(trainingEpochs))
		print('\n*****TRAINING STARTED***** (stop with ctrl-c)\n')

		startTime = time.time()
		epochTime = startTime
		avgCost = 0
		for epoch in range(trainingEpochs):
			batch_xs, batch_ys = getInput(batch_size = batch_size, useDataAugmentation = useDataAugmentation)
			saved=False
			_, gen_cost, summary = sess.run([pretrain_optimizer, gen_l2_loss, lossPretrain_gen], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: True})
			summary_writer.add_summary(summary, epoch)
			avgCost += gen_cost

			if (epoch + 1) % saveInterval == 0:
				print('%05d / %d: last interval: %.02f seconds, %.02f min remaining. avg cost: %.02f' % (epoch+1, trainingEpochs, (time.time() - epochTime), ((trainingEpochs - epoch) * (time.time() - startTime) / epoch / 60.0), (avgCost / outputInterval)))
				epochTime = time.time()
				avgCost = 0
				saved = True
				print(saveModel(gen_cost, genTestImg, test_path+"test_img/")) 

		#training_duration = (time.time() - startTime) / 60.0
		#print('Training needed %.02f minutes.' % (training_duration))
		#sys.stdout.flush()

	except KeyboardInterrupt:
		print("training interrupted")
		sys.stdout.flush()
		with open(basePath + 'test_overview.log', "a") as text_file:
			text_file.write('\ttraining interrupted after %d epochs' % (epoch + 1) + '\n')

	if not saved:
		print(saveModel(gen_cost, genTestImg, test_path+"test_img/"))

	print('\n*****TRAINING FINISHED*****')
	training_duration = (time.time() - startTime) / 60.0
	print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, call the script with command line parameters "out 1  load_model_test %d  load_model_no %d " ' % (load_model_test_new, (save_no-1)) )
	sys.stdout.flush()
	with open(basePath + 'test_overview.log', "a") as text_file:
		text_file.write('\ttraining duration: %.02f minutes' % training_duration + '\n')

### OUTPUT MODE ###

elif outputOnly: #may not work if using tiles smaller than full sim size
	print('*****OUTPUT ONLY*****')

	for layerno in range(0,frame_max-frame_min):
		print('Generating %d' % (layerno))
		if dataDimension == 2:
			generateTestImage(fromSim,layerno, outPath = test_path, imageindex = layerno)
		else:
			print('Not supported at the moment...')	
			#generate3DUni(fromSim,layerno,outPath = test_path, imageindex = layerno)

	print('Test finished, %d outputs written to %s.' % (frame_max-frame_min, test_path) )

