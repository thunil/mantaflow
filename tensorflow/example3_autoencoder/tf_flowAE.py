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
dt			= float(ph.getParam( "dt", 		  0.5 )) 				# step time of training data
#Data and Output
loadPath		 =	 ph.getParam( "loadPath",		 '../data/' ) 	# path to training data
fromSim		   = int(ph.getParam( "fromSim",		 1000 )) 			# range of sim data to use, start index
toSim		   = int(ph.getParam( "toSim",		   -1   )) 			# end index
dataDimension   = int(ph.getParam( "dataDim",		 2 )) 				# dimension of dataset, can be 2 or 3. in case of 3D any scaling will only be applied to H and W (DHW)
numOut			= int(ph.getParam( "numOut",		  200 )) 			# number ouf images to output (from start of sim)
saveOut	  	    = int(ph.getParam( "saveOut",		 False ))>0 		# save output of output mode as .npz in addition to images
loadOut			= int(ph.getParam( "loadOut",		 -1 )) 			# load output from npz to use in output mode instead of tiles. number or output dir, -1 for not use output data
outputImages	=int(ph.getParam( "img",  			  True ))>0			# output images
outputGif		= int(ph.getParam( "gif",  			  False ))>0		# output gif
outputRef		= int(ph.getParam( "ref",			 False ))>0 		# output "real" data for reference in output mode (may not work with 3D)
#models
genModel		=	 ph.getParam( "genModel",		 'gen_test' ) 	# choose generator model
discModel		=	 ph.getParam( "discModel",		 'disc_test' ) 	# choose discriminator model
#Training
learning_rate   = float(ph.getParam( "learningRate",  0.0002 ))
decayLR		    = int(ph.getParam( "decayLR",			 False ))>0 		# decay learning rate?
dropout   		= float(ph.getParam( "dropout",  	  1.0 )) 			# keep prop for all dropout layers during training
dropoutOutput   = float(ph.getParam( "dropoutOutput", dropout )) 		# affects testing, full sim output and progressive output during training
beta			= float(ph.getParam( "adam_beta1",	 0.5 ))			#1. momentum of adam optimizer

weight_dld		= float(ph.getParam( "weight_dld",	1.0)) 			# ? discriminator loss factor ?
k				= float(ph.getParam( "lambda",		  1.0)) 			# influence/weight of l1 term on generator loss
#k2				= float(ph.getParam( "lambda2",		  0.0)) 			# influence/weight of d_loss term on generator loss
#k_f				= float(ph.getParam( "lambda_f",		  1.0)) 			# changing factor of k
batch_size	    = int(ph.getParam( "batchSize",  	  128 ))			# batch size for pretrainig and output, default for batchSizeDisc and batchSizeGen
batch_size_disc = int(ph.getParam( "batchSizeDisc",   batch_size )) 	# batch size for disc runs when training gan
batch_size_gen  = int(ph.getParam( "batchSizeGen",	batch_size )) 	# batch size for gen runs when training gan
trainGAN		= int(ph.getParam( "trainGAN",   	  True ))>0 		# GAN trainng can be switched off to use pretrainig only
trainingEpochs  = int(ph.getParam( "trainingEpochs",  100000 )) 		# for GAN training
discRuns 		= int(ph.getParam( "discRuns",  	  1 )) 				# number of discrimiinator optimizer runs per epoch
genRuns  		= int(ph.getParam( "genRuns",  		  1 )) 				# number or generator optimizer runs per epoch
batch_norm		= int(ph.getParam( "batchNorm",	   True ))>0			# apply batch normalization to conv and deconv layers
bn_decay		= float(ph.getParam( "bnDecay",	   0.999 ))			# decay of batch norm EMA
use_spatialdisc = int(ph.getParam( "use_spatialdisc",		   True )) #use spatial discriminator or not

useVelocities   = int(ph.getParam( "useVelocities",   0  )) 			# use velocities or not

useDataAugmentation = int(ph.getParam( "dataAugmentation", 0 ))		 # use dataAugmentation or not
minScale = float(ph.getParam( "minScale",	  0.85 ))				 # augmentation params...
maxScale = float(ph.getParam( "maxScale",	  1.15 ))
rot	     = int(ph.getParam( "rot",		  2	 ))		#rot: 1: 90 degree rotations; 2: full rotation; else: nop rotation 
flip	 =   int(ph.getParam( "flip",		  1	 ))

#Pretraining
pretrain		= int(ph.getParam( "pretrain",		0 )) 				# train generator with L2 loss before alternating training, number of epochs
#pretrain_disc	= int(ph.getParam( "pretrainDisc",   0 )) 				# train discriminator before alternating training
#pretrain_gen	= int(ph.getParam( "pretrainGen",	0 ))				# train generator using pretrained discriminator before alternating training

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
ADV_flag		= int(ph.getParam( "adv_flag",		   True )) # Tempo parameter, add( or not) advection to pre/back frame to align
change_velocity = int(ph.getParam( "change_velocity",		   False )) 
saveMD          = int(ph.getParam( "saveMetaData", 0 ))      # profiling, add metadata to summary object? warning - only main training for now
overlap         = int(ph.getParam( "overlap",		   3 )) # parameter for 3d unifile output, overlap of voxels

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


############################################discriminator network###############################################################
def disc_binclass(in_low, in_high, reuse=False, use_batch_norm=False, train=None):
	#in_low: low res reference input, same as generator input (condition)
	#in_high: real or generated high res input to classify
	#reuse: variable reuse
	#use_batch_norm: bool, if true batch norm is used in all but the first con layers
	#train: if use_batch_norm, tf bool placeholder
	print("\n\tDiscriminator (conditional binary classifier)")
	with tf.variable_scope("discriminator", reuse=reuse):
		if dataDimension == 2:
			#in_low,_,_ = tf.split(in_low,n_inputChannels,1)
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, 1])).max_depool(height_factor = upRes,width_factor=upRes) #NHWC
			print(in_low)
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, 1])
			filter=[4,4]
			stride = [2]
			stride2 = [2]
		elif dataDimension == 3:
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, 1])).max_depool(depth_factor = upRes,height_factor = upRes,width_factor = upRes) #NDHWC
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, 1]) # dim D is not upscaled
			filter=[4,4,4]
			stride = [2,2]
			stride2 = [2]

		#merge in_low and in_high to [-1, tileSizeHigh, tileSizeHigh, 2]
		gan = GAN(tf.concat([in_low, in_high], axis=-1), bn_decay=bn_decay) #64
		d1,_ = gan.convolutional_layer(32, filter, lrelu, stride=stride2, name="d_c1", reuse=reuse) #32

		d2,_ = gan.convolutional_layer(64, filter, lrelu, stride=stride2, name="d_c2", reuse=reuse, batch_norm=use_batch_norm, train=train) #64

		d3,_ = gan.convolutional_layer(128, filter, lrelu, stride=stride, name="d_c3", reuse=reuse, batch_norm=use_batch_norm, train=train) #128

		d4,_ = gan.convolutional_layer(256, filter, lrelu, stride=[1], name="d_c4", reuse=reuse, batch_norm=use_batch_norm, train=train) #256

		shape=gan.flatten()
		gan.fully_connected_layer(1, None, name="d_l5")

		print("\tDOFs: %d " % gan.getDOFs())
		return gan.y()
		

############################################ Tempo discriminator network ############################################################
def disc_binclass_cond_tempo(in_high, n_t_channels=3, reuse=False, use_batch_norm=False, train=None):
	# NO in_low: low res reference input, same as generator input (no condition)
	# in_high: real or generated high res input to classify, shape should be batch, dim_z, dim_y, dim_x, channels
	# reuse: variable reuse
	# use_batch_norm: bool, if true batch norm is used in all but the first con layers
	# train: if use_batch_norm, tf bool placeholder
	print("\n\tDiscriminator for Tempo (conditional binary classifier)")
	print("\n\tTempo, nearby frames packed as channels, number %d" % n_t_channels)
	with tf.variable_scope("discriminatorTempo", reuse=reuse):
		if dataDimension == 2:
			# in_low,_,_ = tf.split(in_low,n_inputChannels,1)
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, n_t_channels])
			filter=[4,4]
			stride = [2]
			stride2 = [2]
		elif dataDimension == 3:
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, n_t_channels]) # dim D is not upscaled
			filter=[4,4,4]
			stride = [2,2]
			stride2 = [2]

		# merge in_low and in_high to [-1, tileSizeHigh, tileSizeHigh, 2]
		gan = GAN(in_high, bn_decay=bn_decay)  # 64
		t1, _ = gan.convolutional_layer(32, filter, lrelu, stride=stride2, name="t_c1", reuse=reuse)  # 32
		t2, _ = gan.convolutional_layer(64, filter, lrelu, stride=stride2, name="t_c2", reuse=reuse,
										batch_norm=use_batch_norm, train=train)  # 64
		t3, _ = gan.convolutional_layer(128, filter, lrelu, stride=stride, name="t_c3", reuse=reuse,
										batch_norm=use_batch_norm, train=train)  # 128
		t4, _ = gan.convolutional_layer(256, filter, lrelu, stride=[1], name="t_c4", reuse=reuse,
										batch_norm=use_batch_norm, train=train)  # 256
		shape = gan.flatten()
		gan.fully_connected_layer(1, None, name="t_l5")

		print("\tDOFs: %d " % gan.getDOFs())
		return gan.y()

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

############################################disc_test###############################################################
def disc_test(in_low, in_high, reuse=False, use_batch_norm=False, train=None):
	print("\n\tDiscriminator-test")
	with tf.variable_scope("discriminator_test", reuse=reuse):
		if dataDimension == 2:
			#in_low,_,_ = tf.split(in_low,n_inputChannels,1)
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, 1])).max_depool(height_factor = upRes,width_factor = upRes) #NHWC
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, 1])
			filter=[4,4]
			stride2 = [2]
		elif dataDimension == 3:
			shape = tf.shape(in_low)
			in_low = tf.slice(in_low,[0,0],[shape[0],int(n_input/n_inputChannels)])
			in_low = GAN(tf.reshape(in_low, shape=[-1, tileSizeLow, tileSizeLow, tileSizeLow, 1])).max_depool(depth_factor = upRes,height_factor = upRes,width_factor = upRes) #NDHWC
			in_high = tf.reshape(in_high, shape=[-1, tileSizeHigh, tileSizeHigh, tileSizeHigh, 1]) # dim D is not upscaled
			filter=[4,4,4]
			stride2 = [2]

		#merge in_low and in_high to [-1, tileSizeHigh, tileSizeHigh, 2]
		gan = GAN(tf.concat([in_low, in_high], axis=-1), bn_decay=bn_decay) #64
		d1,_ = gan.convolutional_layer(32, filter, lrelu, stride=stride2, name="d_c1", reuse=reuse) #32
		shape=gan.flatten()
		gan.fully_connected_layer(1, None, name="d_l5")
		if dataDimension == 2:
			d2 = tf.constant(1., shape = [batch_size, tileSizeLow,tileSizeLow,64])
			d3 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),128])	
			d4 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),256])
		elif dataDimension == 3:
			d2 = tf.constant(1., shape = [batch_size, tileSizeLow,tileSizeLow,tileSizeLow,64])
			d3 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),int(tileSizeLow/2),128])	
			d4 = tf.constant(1., shape = [batch_size, int(tileSizeLow/2),int(tileSizeLow/2),int(tileSizeLow/2),256])
		print("\tDOFs: %d " % gan.getDOFs())
		return gan.y()

#change used models for gen and disc here #other models in NNmodels.py
gen_model = locals()[genModel]
disc_model = locals()[discModel]
disc_time_model = disc_binclass_cond_tempo # tempo dis currently fixed

#set up GAN structure
bn=batch_norm
#training or testing for batch norm
train = tf.placeholder(tf.bool)

if not outputOnly: #setup for training
	gen_part = gen_model(x, use_batch_norm=bn, train=train)
	if use_spatialdisc:
		disc  = disc_model(x_disc, y, use_batch_norm=bn, train=train)
		gen   = disc_model(x_disc, gen_part, reuse=True, use_batch_norm=bn, train=train)
	if genTestImg > -1: sampler = gen_part
else: #setup for generating output with trained model
	sampler = gen_model(x, use_batch_norm=bn, train=False)


sys.stdout.flush()
#exit(1) #used for net layout check

if not outputOnly:
	#for discriminator [0,1] output
	if use_spatialdisc:
		disc_sigmoid = tf.reduce_mean(tf.nn.sigmoid(disc))
		gen_sigmoid = tf.reduce_mean(tf.nn.sigmoid(gen))

		# loss of the discriminator with real input 
		disc_loss_disc = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=disc, labels=tf.ones_like(disc)))
		#loss of the discriminator with input from generator
		disc_loss_gen = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=gen, labels=tf.zeros_like(gen)))
		disc_loss = disc_loss_disc * weight_dld + disc_loss_gen
		#loss of the generator
		gen_loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=gen, labels=tf.ones_like(gen)))
	
	else:
		gen_loss = tf.zeros([1])
		#disc_loss_layer = tf.zeros([1])

	#additional generator losses
	gen_l2_loss = tf.nn.l2_loss(y - gen_part)
	gen_l1_loss = tf.reduce_mean(tf.abs(y - gen_part)) #use mean to normalize w.r.t. output dims. tf.reduce_sum(tf.abs(y - gen_part))

	#uses sigmoid cross entropy and l1 - see cGAN paper
	gen_loss_complete = gen_loss + gen_l1_loss*kk 

	# set up decaying learning rate, if enabled
	lr_global_step = tf.Variable(0, trainable=False)
	learning_rate_scalar = learning_rate
	if decayLR:
		learning_rate = tf.train.polynomial_decay(learning_rate, lr_global_step, trainingEpochs//2, learning_rate_scalar*0.05, power=1.1)

	update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
	gen_update_ops = update_ops[:]
	ori_gen_update_ops = update_ops[:]
	pre_update_ops = update_ops[:]
	#print(update_ops)

	#variables to be used in the different otimization steps
	vars = tf.trainable_variables()
	g_var = [var for var in vars if "g_" in var.name]
	if use_spatialdisc:
		dis_update_ops = update_ops[:]
		d_var = [var for var in vars if "d_" in var.name]
				
	if use_spatialdisc:
		with tf.control_dependencies(dis_update_ops):
			#optimizer for discriminator, uses combined loss, can only change variables of the disriminator
			disc_optimizer_adam = tf.train.AdamOptimizer(learning_rate, beta1=beta)
			disc_optimizer = disc_optimizer_adam.minimize(disc_loss, var_list=d_var)

	with tf.control_dependencies(gen_update_ops): 
		# optimizer for generator, can only change variables of the generator,
		gen_optimizer = tf.train.AdamOptimizer(learning_rate, beta1=beta).minimize(gen_loss_complete, var_list=g_var)

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
	# create a summary to monitor cost tensor
	#training losses
	if use_spatialdisc:
		lossTrain_disc  = tf.summary.scalar("discriminator-loss train",     disc_loss)
		lossTrain_gen  = tf.summary.scalar("generator-loss train",     gen_loss)

	#testing losses
	if use_spatialdisc:
		lossTest_disc_disc   = tf.summary.scalar("discriminator-loss test real", disc_loss_disc)
		lossTest_disc_gen   = tf.summary.scalar("discriminator-loss test generated", disc_loss_gen)
		lossTest_disc = tf.summary.scalar("discriminator-loss test", disc_loss)
		lossTest_gen   = tf.summary.scalar("generator-loss test", gen_loss)

	#discriminator output [0,1] for real input
	if use_spatialdisc:
		outTrain_disc_real = tf.summary.scalar("discriminator-out train", disc_sigmoid)
		outTrain_disc_gen = tf.summary.scalar("generator-out train", gen_sigmoid)

	#discriminator output [0,1] for generated input
	if use_spatialdisc:
		outTest_disc_real = tf.summary.scalar("discriminator-out test", disc_sigmoid)
		outTest_disc_gen = tf.summary.scalar("generator-out test", gen_sigmoid)

	#pretrain losses
	#lossPretrain_disc  = tf.summary.scalar("discriminator_pretrain_loss",     disc_loss)
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
	if pretrain>0: # or pretrain_gen > 0 or pretrain_disc>0:
		os.makedirs(test_path+'pretrain_test_img/') # NT_DEBUG ?

def modifyVel(Dens,Vel):
	return velout # not active right now...

tCnt=0
def getInput(index = 1, randomtile = True, isTraining = True, batch_size = 1, useDataAugmentation = False, modifyvelocity = False, useVelocities = False):
	global tCnt
	if randomtile == False:
		batch_xs, batch_ys = tiCr.getFrameTiles(index) 
	else:
		batch_xs, batch_ys = tiCr.selectRandomTiles(selectionSize = batch_size, augment=useDataAugmentation)	
		
	if useVelocities and modifyvelocity:
		Densinput = batch_xs[:,:,:,:,0:1]
		Velinput = batch_xs[:,:,:,:,1:4]
		Veloutput = modifyVel(Densinput, Velinput)
		batch_xs = np.concatenate((Densinput, Veloutput), axis = 4)

	batch_xs = np.reshape(batch_xs, (-1, n_input))
	batch_ys = np.reshape(batch_ys, (-1, n_output))
	return batch_xs, batch_ys

#evaluate the generator (sampler) on the first step of the first simulation and output result
def generateTestImage(sim_no = fromSim, frame_no = 1, outPath = test_path,imageindex = 0, modifyvelocity = False):
	if (not outputOnly):
		batch_xs, _ = getInput(randomtile = False, index = (sim_no-fromSim)*frame_max + frame_no, modifyvelocity = modifyvelocity, useVelocities = useVelocities)
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

def generate3DUni(sim_no = fromSim, frame_no = 1, outPath = test_path,imageindex = 0):
	if dataDimension == 2:
		print("ERROR: only for 3D Uni files output!")	
		exit(1)

	if (overlap*2 > tileSizeLow) or (tileSizeLow > simLowLength):
		print("Wrong parameters for 3d output!")	
		exit(1)
	batch_xs = inputx[frame_no]
	if useVelocities and change_velocity:
		batch_xs = np.reshape(batch_xs,[1,simLowLength,simLowWidth,simLowHeight,-1])
		Densinput = batch_xs[:,:,:,:,0:1]
		Velinput = batch_xs[:,:,:,:,1:4]
		Veloutput = modifyVel(Densinput, Velinput)
		batch_xs = np.concatenate((Densinput, Veloutput), axis = 4)
		batch_xs = np.reshape(batch_xs,[simLowLength,simLowWidth,simLowHeight,4])
	tiles = []
	batch_xs=np.reshape(batch_xs,[simLowLength,simLowWidth,simLowHeight,-1])

	lengthnum = ((simLowLength-overlap*2+tileSizeLow-overlap*2-1)//(tileSizeLow-overlap*2))
	widthnum = ((simLowWidth-overlap*2+tileSizeLow-overlap*2-1)//(tileSizeLow-overlap*2))
	heightnum = ((simLowHeight-overlap*2+tileSizeLow-overlap*2-1)//(tileSizeLow-overlap*2))

	for i in range(lengthnum):
		for j in range(widthnum):
			for k in range(heightnum):
				ifrom = (tileSizeLow-overlap*2)*i
				ito = (tileSizeLow-overlap*2)*i+tileSizeLow
				jfrom = (tileSizeLow-overlap*2)*j
				jto = (tileSizeLow-overlap*2)*j+tileSizeLow
				kfrom = (tileSizeLow-overlap*2)*k
				kto = (tileSizeLow-overlap*2)*k+tileSizeLow
				if ito >simLowLength:
					ifrom = simLowLength-tileSizeLow
					ito = simLowLength
				if jto >simLowWidth:
					jfrom = simLowWidth-tileSizeLow
					jto = simLowWidth
				if kto >simLowHeight:
					kfrom = simLowHeight-tileSizeLow
					kto = simLowHeight
				low = batch_xs[ifrom:ito, jfrom:jto, kfrom:kto, :]
				tiles.append(low)
	batch_xs = np.array(tiles)
	resultTiles = []
	for tileno in range(batch_xs.shape[0]):
		batch_xs_in = np.reshape(batch_xs[tileno],[-1, n_input])
		results = sess.run(sampler, feed_dict={x: batch_xs_in, keep_prob: dropoutOutput, train : False})
		results = np.array(results)
		resultTiles.extend(results)
	resultTiles = np.array(resultTiles)
	resulttiles = np.reshape(resultTiles,[resultTiles.shape[0],tileSizeHigh,tileSizeHigh,tileSizeHigh])
	high = np.zeros([simLowLength*upRes,simLowWidth*upRes,simLowHeight*upRes])
	for i in range(lengthnum):
		for j in range(widthnum):
			for k in range(heightnum):
				ihighfrom = (tileSizeLow-overlap*2)*upRes*(i-1)+(tileSizeLow-overlap)*upRes
				ihighto = ihighfrom + (tileSizeLow-overlap*2)*upRes
				jhighfrom = (tileSizeLow-overlap*2)*upRes*(j-1)+(tileSizeLow-overlap)*upRes
				jhighto = jhighfrom+(tileSizeLow-overlap*2)*upRes
				khighfrom = (tileSizeLow-overlap*2)*upRes*(k-1)+(tileSizeLow-overlap)*upRes
				khighto = khighfrom+(tileSizeLow-overlap*2)*upRes
				ifrom = overlap*upRes
				ito = (tileSizeLow-overlap)*upRes
				jfrom = overlap*upRes
				jto = (tileSizeLow-overlap)*upRes
				kfrom = overlap*upRes
				kto = (tileSizeLow-overlap)*upRes
				if i == 0:
					ifrom = 0
					ito = (tileSizeLow-overlap)*upRes
					ihighfrom = 0
					ihighto = (tileSizeLow-overlap)*upRes
				if j == 0:
					jfrom = 0
					jto = (tileSizeLow-overlap)*upRes
					jhighfrom = 0
					jhighto = (tileSizeLow-overlap)*upRes
				if k == 0:
					kfrom = 0
					kto = (tileSizeLow-overlap)*upRes
					khighfrom = 0
					khighto = (tileSizeLow-overlap)*upRes
				if i == lengthnum-1:
					ifrom = overlap*upRes
					ito = tileSizeLow*upRes
					ihighfrom = simLowLength*upRes-tileSizeLow*upRes+overlap*upRes
					ihighto = simLowLength*upRes
				if j == widthnum-1:
					jfrom = overlap*upRes
					jto = tileSizeLow*upRes
					jhighfrom = simLowWidth*upRes-tileSizeLow*upRes+overlap*upRes
					jhighto = simLowWidth*upRes
				if k == heightnum-1:
					kfrom = overlap*upRes
					kto = tileSizeLow*upRes
					khighfrom = simLowHeight*upRes-tileSizeLow*upRes+overlap*upRes
					khighto = simLowHeight*upRes
				high[ihighfrom: ihighto, jhighfrom:jhighto, khighfrom:khighto] = resulttiles[i*widthnum*heightnum+j*heightnum+k][ifrom:ito,jfrom:jto,kfrom:kto]

		high = np.reshape(high,[simLowLength*upRes,simLowWidth*upRes,simLowHeight*upRes])
		
		head, _ = uniio.readUni(loadPath + "sim_%04d/density_high_%04d.uni"%(sim_no,frame_no+frame_min))
		head['dimX'] = simLowHeight*upRes
		head['dimY'] = simLowWidth*upRes
		head['dimZ'] = simLowLength*upRes
		uniio.writeUni(outPath+'source_%04d.uni'%(frame_no+frame_min), head, high)
		if(0): # save png
			save_img_3d( outPath + 'source_{:04d}.png'.format(imageCounter),high )

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
		text_file.write('\t{} Epochs, gen: {}, disc: {}'.format(trainingEpochs, gen_model.__name__, disc_model.__name__) + loaded_model + '\n')
		text_file.write('\tgen-runs: {}, disc-runs: {}, lambda: {}, dropout: {:.4f}({:.4f})'.format(genRuns, discRuns, k, dropout, dropoutOutput) + '\n')
	else:
		text_file.write('Output:' + loaded_model + ' (' + test_path[-28:-1] + ')\n')
		text_file.write('\ttile size: {}, seed: {}, dropout-out: {:.4f}'.format(tileSizeLow, randSeed, dropoutOutput) + '\n')

	
	
#train generator using L2 loss
if (not outputOnly) and pretrain>0:
	print('\t Generator using L2')
	print('{} epochs\n'.format(pretrain))
	startTime = time.time()
	epochTime = startTime
	avgCost = 0
	for epoch in range(pretrain):
		batch_xs, batch_ys = getInput(batch_size = batch_size, useDataAugmentation = useDataAugmentation, useVelocities = useVelocities)
		_, gen_cost, summary = sess.run([pretrain_optimizer, gen_l2_loss, lossPretrain_gen], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: True})
		summary_writer.add_summary(summary, epoch)
		avgCost += gen_cost

		if (epoch + 1) % saveInterval == 0:
			print('%05d / %d: last interval: %.02f seconds, %.02f min remaining. avg cost: %.02f' % (epoch+1, pretrain, (time.time() - epochTime), ((pretrain - epoch) * (time.time() - startTime) / epoch / 60.0), (avgCost / outputInterval)))
			epochTime = time.time()
			avgCost = 0
			print(saveModel(gen_cost, genTestImg, test_path+"pretrain_test_img/")) 

	print(saveModel(gen_cost, genTestImg, test_path+"pretrain_test_img/"))
	training_duration = (time.time() - startTime) / 60.0
	print('Training needed %.02f minutes.' % (training_duration))
	sys.stdout.flush()
	
# ---------------------------------------------
# ---------------------------------------------
# START TRAINING
training_duration = 0.0
cost = 0.0

if not outputOnly and trainGAN:
	try:
		print('\n*****TRAINING STARTED*****\n')
		print('(stop with ctrl-c)')
		avgCost_disc = 0
		avgCost_gen = 0
		avgL1Cost_gen = 0
		avgOut_disc = 0
		avgOut_gen = 0

		avgTestCost_disc_real = 0
		avgTestCost_disc_gen = 0
		avgTestCost_gen = 0
		avgTestOut_disc_real = 0
		avgTestOut_disc_gen = 0
		tests = 0
		startTime = time.time()
		intervalTime = startTime
		lastOut = 1
		lastSave = 1
		lastCost = 1e10
		saved = False
		saveMsg = ''
		kkin = k

		disc_cost = 0
		gen_cost = 0
		
		avgTemCost_gen = 0
		avgTemCost_gen_l = 0
		avgTemCost_disc = 0

		avgOut_disc_t = 0
		avgOut_gen_t = 0
		avgTestCost_disc_real_t = 0
		avgTestOut_disc_real_t = 0
		avgTestCost_disc_gen_t = 0
		avgTestOut_disc_gen_t = 0
		avgTestCost_gen_t = 0
		avgTestCost_gen_t_l = 0
		
		for epoch in range(trainingEpochs):
			lrgs = max(0, epoch-(trainingEpochs//2)) # LR counter, start decay at half time... (if enabled) 

			run_options = None; run_metadata = None
			if saveMD:
				run_options = tf.RunOptions(trace_level=tf.RunOptions.FULL_TRACE)
				#run_options = tf.RunOptions(trace_level=tf.RunOptions.FULL_TRACE, output_partition_graphs=True)
				run_metadata = tf.RunMetadata()


			# TRAIN MODEL
			# discriminator variables; with real and generated input
			if use_spatialdisc:
				for runs in range(discRuns):
					batch_xs, batch_ys = getInput(batch_size = batch_size_disc, useDataAugmentation = useDataAugmentation, useVelocities = useVelocities)
					_, disc_cost, summary,disc_sig,gen_sig = sess.run([disc_optimizer, disc_loss, lossTrain_disc,disc_sigmoid,gen_sigmoid], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: True, lr_global_step: lrgs}     , options=run_options, run_metadata=run_metadata )
					avgCost_disc += disc_cost
					summary_writer.add_summary(summary, epoch)
					if saveMD: summary_writer.add_run_metadata(run_metadata, 'dstep%d' % epoch)

			# generator variables
			for runs in range(genRuns):
				batch_xs, batch_ys = getInput(batch_size = batch_size_disc, useDataAugmentation = useDataAugmentation, useVelocities = useVelocities)
				
				train_dict = {x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: True, kk: kkin,
							lr_global_step: lrgs}
				if use_spatialdisc:
					getlist = [gen_optimizer, gen_loss, gen_l1_loss, lossTrain_gen, gen_l2_loss]
				else:
					getlist = [gen_optimizer, gen_l1_loss, gen_l2_loss]

				result_list = sess.run(getlist, feed_dict=train_dict, options=run_options, run_metadata=run_metadata)

				if use_spatialdisc:
					_, gen_cost, gen_l1_cost, summary, gen_l2_cost = result_list
				else:
					_, gen_l1_cost, gen_l2_cost = result_list
				gen_tem_cost = 0
				gen_tem_cost_l = 0

				avgL1Cost_gen += gen_l1_cost
				avgTemCost_gen += gen_tem_cost
				avgTemCost_gen_l += gen_tem_cost_l
				if use_spatialdisc:
					avgCost_gen += gen_cost
					summary_writer.add_summary(summary, epoch)
				if saveMD: summary_writer.add_run_metadata(run_metadata, 'gstep%d' % epoch)


			# save model
			if ((disc_cost+gen_cost < lastCost) or alwaysSave) and (lastSave >= saveInterval):
				lastSave = 1
				lastCost = disc_cost+gen_cost
				saveMsg = saveModel(lastCost)
				saved = True
			else:
				lastSave += 1
				saved = False

			# test model
			if (epoch + 1) % testInterval == 0:
				if use_spatialdisc:
					# gather statistics from training
					# not yet part of testing!
					batch_xs, batch_ys = getInput(batch_size = numTests, useVelocities = useVelocities)
					disc_out, summary_disc_out, gen_out, summary_gen_out = sess.run([disc_sigmoid, outTrain_disc_real, gen_sigmoid, outTrain_disc_gen], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropout, train: False})
					summary_writer.add_summary(summary_disc_out, epoch)
					summary_writer.add_summary(summary_gen_out, epoch)
					avgOut_disc += disc_out
					avgOut_gen += gen_out

					# testing starts here...
					# get test data
					batch_xs, batch_ys = getInput(batch_size = numTests, isTraining=False, useVelocities = useVelocities)
					#disc with real imput
					disc_out_real, summary_test_out, disc_test_cost_real, summary_test = sess.run([disc_sigmoid, outTest_disc_real, disc_loss_disc, lossTest_disc_disc], feed_dict={x: batch_xs, x_disc: batch_xs, y: batch_ys, keep_prob: dropoutOutput, train: False})
					summary_writer.add_summary(summary_test, epoch)
					summary_writer.add_summary(summary_test_out, epoch)
					avgTestCost_disc_real += disc_test_cost_real
					avgTestOut_disc_real += disc_out_real
					#disc with generated input
					disc_out_gen, summary_test_out, disc_test_cost_gen, summary_test = sess.run([gen_sigmoid, outTest_disc_gen, disc_loss_gen, lossTest_disc_gen], feed_dict={x: batch_xs, x_disc: batch_xs, keep_prob: dropoutOutput, train: False})
					summary_writer.add_summary(summary_test, epoch)
					summary_writer.add_summary(summary_test_out, epoch)
					avgTestCost_disc_gen += disc_test_cost_gen
					avgTestOut_disc_gen += disc_out_gen
									
				#gen
				train_dict = {x: batch_xs, x_disc: batch_xs, keep_prob: dropoutOutput, train: False}
				gen_test_cost, summary_test = sess.run([gen_loss, lossTest_gen], feed_dict=train_dict)
				summary_writer.add_summary(summary_test, epoch)
				avgTestCost_gen += gen_test_cost

				tests += 1

			# output statistics
			if (epoch + 1) % outputInterval == 0:
				#training average costs
				avgCost_disc /= (outputInterval * discRuns)
				avgCost_gen /= (outputInterval * genRuns)
				avgL1Cost_gen /= (outputInterval * genRuns)
				#test average costs
				if not (tests == 0):
					avgOut_disc /= tests
					avgOut_gen /= tests
					avgTestCost_disc_real /= tests
					avgTestCost_disc_gen /= tests
					avgTestCost_gen /= tests
					avgTestOut_disc_real /= tests
					avgTestOut_disc_gen /= tests
						
				print('\nEpoch {:05d}/{}, Cost:'.format((epoch + 1), trainingEpochs))
				print('\tdisc: loss: train_loss={:.6f} - test-real={:.6f} - test-generated={:.6f}, out: train={:.6f} - test={:.6f}'.
					format(avgCost_disc, avgTestCost_disc_real, avgTestCost_disc_gen, avgOut_disc, avgTestOut_disc_real))
				print('\tT D : loss[ -train (total={:.6f}), -test (real&1={:.6f}) (generated&0={:.6f})]'.
					format(avgTemCost_disc, avgTestCost_disc_real_t, avgTestCost_disc_gen_t))
				print('\t	sigmoidout[ -test (real&1={:.6f}) (generated&0={:.6f})'.
					format(avgTestOut_disc_real_t, avgTestOut_disc_gen_t))
				print('\t gen: loss: train={:.6f} - L1(*k)={:.3f} - test={:.6f}, DS out: train={:.6f} - test={:.6f}'
					.format(avgCost_gen, avgL1Cost_gen * k, avgTestCost_gen, avgOut_gen, avgTestOut_disc_gen))
				if use_spatialdisc:
					print('\tdisc: loss: disc=%f'%(disc_sig))
					print('\tgen: loss: gen=%f'%(gen_sig))
								
				print('\t l1_cost: %f'%(gen_l1_cost))
				
				epochTime = (time.time() - startTime) / (epoch + 1)
				print('\t{} epochs took {:.2f} seconds. (Est. next: {})'.format(outputInterval, (time.time() - intervalTime), time.ctime(time.time() + outputInterval * epochTime)))
				remainingTime = (trainingEpochs - epoch) * epochTime
				print('\tEstimated remaining time: {:.2f} minutes. (Est. end: {})'.format(remainingTime / 60.0, time.ctime(time.time() + remainingTime)))
				if saved:
					print('\t' + saveMsg) # print save massage here for clarity
				if genTestImg > -1:
					generateTestImage(outPath = test_path+'test_img/', imageindex = image_no)
					image_no +=1
				sys.stdout.flush()
				intervalTime = time.time()
				avgCost_disc = 0
				avgCost_gen = 0
				avgL1Cost_gen = 0
				avgOut_disc = 0
				avgOut_gen = 0
				avgTestCost_disc_real = 0
				avgTestCost_disc_gen = 0
				avgTestCost_gen = 0
				avgTestOut_disc_real = 0
				avgTestOut_disc_gen = 0
				tests = 0
				lastOut = 0

			lastOut +=1

	except KeyboardInterrupt:
		print("training interrupted")
		sys.stdout.flush()
		with open(basePath + 'test_overview.log', "a") as text_file:
			text_file.write('\ttraining interrupted after %d epochs' % (epoch + 1) + '\n')

	print('\n*****TRAINING FINISHED*****')
	training_duration = (time.time() - startTime) / 60.0
	print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, call the script with command line parameters "out=1 load_model_test=%d  load_model_no=%d " ' % (load_model_test_new, save_no) )
	sys.stdout.flush()
	with open(basePath + 'test_overview.log', "a") as text_file:
		text_file.write('\ttraining duration: %.02f minutes' % training_duration + '\n')


### OUTPUT MODE ###

elif outputOnly: #may not work if using tiles smaller than full sim size
	print('*****OUTPUT ONLY*****')

	for layerno in range(0,frame_max-frame_min):
		print('Generating %d' % (layerno))
		if dataDimension == 2:
			generateTestImage(fromSim,layerno,outPath = test_path, imageindex = layerno, modifyvelocity = change_velocity)
		else:
			generate3DUni(fromSim,layerno,outPath = test_path, imageindex = layerno)

	print('Test finished, %d outputs written to %s.' % (frame_max-frame_min, test_path) )

