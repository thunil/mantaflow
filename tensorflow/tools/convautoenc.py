import numpy as np
import tensorflow as tf
import math

class ConvolutionalAutoEncoder(object):
	#---------------------------------------------------------------------------------
	def __init__(self, _image):
		self.layer = _image
		self.batch_size = tf.shape(_image)[0]
		self.DOFs = 0
		# stack 
		self.preFlatShapes = []
		self.weight_stack = []
		self.layer_num = 0
		print ("CAE Input shape: {}".format(self.layer.get_shape())) 
	
	#---------------------------------------------------------------------------------
	# thanks to http://robromijnders.github.io/tensorflow_basic/
	def weight_image(self):
		W = self.weight_stack[-1]
		# compute size of the image
		s = W.get_shape()
		out_channels = 1
		if int(s[3]) % 3 == 0:
			out_channels = 3
		print ("Shape {}".format(s))
		weight_patches = int( int(s[2]) * int(s[3]) / out_channels ) # e.g. the number of [3,3] patches in a CNN
		side_length = int(math.ceil(math.sqrt(weight_patches))) # image side length (in patches)
		image_patches = side_length * side_length # max number of patches that fit in the image
		# split into per filter weights
		ws = []
		ws_dim3 = tf.split(3, int(s[3] / out_channels), W) # e.g. [ [3,3,3,1], [3,3,3,1], ... ]
		for w in ws_dim3:
			# split these further
			ws.extend(tf.split(2, s[2], w))  # e.g. [ [3,3,1,1], [3,3,1,1], ... ]
		# pad image
		padding = image_patches - weight_patches
		for i in range(padding):
			ws.append(tf.zeros([s[0], s[1], 1, out_channels]))
		# build rows of image
		rows = []
		for i in range(side_length):
			start = i * side_length
			end = start + side_length
			rows.append(tf.concat(0, ws[start:end]))
		# combine rows to image
		image = tf.concat(1, rows) # [sidelength * ]
		s = [int(image.get_shape()[0]), int(image.get_shape()[1])]
		image = tf.reshape(image, [1, s[0], s[1], out_channels])
		image = tf.image.resize_images(image, [int(s[1] * 50), int(s[0] * 50)], 1)
		image_tag = "l" + str(self.layer_num) + "_weight_image"
		tf.image_summary(image_tag, image)
		print ("Image Summary: save weights as image")
		
	#---------------------------------------------------------------------------------
	def convolutional_layer(self, _filterSpread, _patchShape, activation_function=tf.nn.tanh):
		self.layer_num += 1
		# set the input and output dimension
		inChannels = int(self.layer.get_shape()[3])
		outChannels = inChannels * _filterSpread
		# create a weight matrix
		W = self.weight_variable([_patchShape[0], _patchShape[1], inChannels, outChannels])
		self.DOFs += _patchShape[0]* _patchShape[1]* inChannels* outChannels
		self.weight_stack.append(W)
		# create a bias vector
		b = self.bias_variable([outChannels])
		self.DOFs += outChannels
		# feed forward step
		self.layer = activation_function(self.conv2d(self.layer, W) + b)
		# user output
		print ("Convolutional Layer {} ({}) : {}".format(W.get_shape(), activation_function.__name__,self.layer.get_shape()))
		return self.layer

	#---------------------------------------------------------------------------------
	# 2 x 2 max pool operation
	def max_pool(self, window_size=[2,2], window_stride=[2,2]):
		self.layer = tf.nn.max_pool(self.layer, ksize=[1, window_size[0], window_size[1], 1], strides=[1, window_stride[0], window_stride[1], 1], padding="VALID")
		# user output
		print ("Max Pool {}: {}".format(window_size, self.layer.get_shape())) 
		return self.layer
	
	#---------------------------------------------------------------------------------
	def avg_pool(self, window_size=[2,2], window_stride=[2,2]):
		self.layer = tf.nn.avg_pool(self.layer, ksize=[1, window_size[0], window_size[1], 1], strides=[1, window_stride[0], window_stride[1], 1], padding="VALID")
		# user output
		print ("Avg Pool {}: {}".format(window_size, self.layer.get_shape()))
		return self.layer
	
	#---------------------------------------------------------------------------------
	# make layer flat
	# e.G. [1, 4, 4, 2] -> [1, 32]
	def flatten(self):
		# get unflat shape
		layerShape = self.layer.get_shape()
		self.preFlatShapes.append(layerShape)
		# compute flat size
		flatSize = int(layerShape[1]) * int(layerShape[2]) * int(layerShape[3])
		# make flat
		self.layer = tf.reshape(self.layer, [-1, flatSize])
		# user output
		print ("Flatten: {}".format(self.layer.get_shape())) 
		return flatSize
	
	#---------------------------------------------------------------------------------
	def fully_connected_layer(self, _numHidden, _act):
		self.layer_num += 1
		# get previous layer size
		numInput = int(self.layer.get_shape()[1])
		# build layer variables
		W = self.weight_variable([numInput, _numHidden])
		b = self.bias_variable([_numHidden])
		self.DOFs += numInput*_numHidden + _numHidden
		# activate
		self.layer = _act(tf.matmul(self.layer, W) + b)  
		# user output
		print ("Fully Connected Layer: {}".format(self.layer.get_shape()))
		return self.layer
	
	#---------------------------------------------------------------------------------
	# make layer 3D (from previously stored
	# e.G. [1, 32] -> [1, 4, 4, 2]
	def unflatten(self):
		unflatShape = self.preFlatShapes.pop()
		unflatShape = [-1, int(unflatShape[1]), int(unflatShape[2]), int(unflatShape[3])]
		self.layer = tf.reshape(self.layer, unflatShape)
		print ("Unflatten: {}".format(self.layer.get_shape()) )
		return self.layer

	#---------------------------------------------------------------------------------
	# inverse of 2 x 2 max pool
	def max_depool(self, window_size=[2, 2], window_stride=[2,2]):
		outWidth = self.layer.get_shape()[2] * window_stride[0] + window_size[0] - window_stride[0]
		outHeight = self.layer.get_shape()[1] * window_stride[1] + window_size[1] -  window_stride[1]
		self.layer = tf.image.resize_images(self.layer, [int(outHeight), int(outWidth)], 1)
		print ("Max Depool {}: {}".format(window_size, self.layer.get_shape()))
		return self.layer

	#---------------------------------------------------------------------------------
	def avg_depool(self, window_size=[2, 2], window_stride=[2,2]):
		outWidth = self.layer.get_shape()[2] * window_stride[0] + window_size[0] - window_stride[0]
		outHeight = self.layer.get_shape()[1] * window_stride[1] + window_size[1] -  window_stride[1]
		self.layer = tf.image.resize_images(self.layer, [int(outHeight), int(outWidth)], 0)
		print ("Avg Depool {}: {}".format(window_size, self.layer.get_shape()))
		return self.layer

	#---------------------------------------------------------------------------------
	def deconvolutional_layer(self, _filterSpread, _patchShape, activation_function=tf.nn.tanh):
		self.layer_num += 1
		# if mode "VALID" is enabled the output shape is a bit larger than the input, because of the patch not beeing padded in convolution
		xDimension = int(self.layer.get_shape()[1]) #+ (int(_patchShape[0]) - 1)
		yDimension = int(self.layer.get_shape()[2]) #+ (int(_patchShape[1]) - 1)
		# spread channels
		inChannels = int(self.layer.get_shape()[3])
		outChannels = int(inChannels / _filterSpread) # must always come out even
		# create a weight matrix
		W = self.weight_variable([_patchShape[0], _patchShape[1], outChannels, inChannels]) 
		self.DOFs += _patchShape[0]* _patchShape[1]* outChannels* inChannels
		# create a bias vector
		b = self.bias_variable([outChannels])
		self.DOFs += outChannels
		# feed forward step
		self.layer = activation_function(self.deconv2d(self.layer, W, [self.batch_size, xDimension, yDimension, outChannels]) + b)
		self.layer = tf.reshape(self.layer, [-1, xDimension, yDimension, outChannels])
		# user output
		print ("Deconvolutional Layer ({}): {}".format(activation_function.__name__, self.layer.get_shape())) 
		return self.layer

	#---------------------------------------------------------------------------------
	def y(self):
		return self.layer

	#---------------------------------------------------------------------------------
	def getDOFs(self):
		return self.DOFs

	#---------------------------------------------------------------------------------
	# generate random valued weight field
	def weight_variable(self, shape):
		initial = tf.truncated_normal(shape=shape, mean=0.0, stddev=0.1)
		return tf.Variable(initial)

	#---------------------------------------------------------------------------------
	# gemerate biases for the nodes
	def bias_variable(self, shape):
		initial = tf.constant(0.1, shape=shape)
		return tf.Variable(initial)

	#---------------------------------------------------------------------------------
	def conv2d(self, x, W):
		return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding="SAME")

	#---------------------------------------------------------------------------------
	def deconv2d(self, x, W, output_shape):
		return tf.nn.conv2d_transpose(x, W, output_shape=output_shape, strides=[1,1,1,1], padding="SAME")

	def variable_summaries(self, var, name):
		"""Attach a lot of summaries to a Tensor."""
		with tf.name_scope('summaries'):
			mean = tf.reduce_mean(var)
			tf.scalar_summary('mean/' + name, mean)
			with tf.name_scope('stddev'):
				stddev = tf.sqrt(tf.reduce_sum(tf.square(var - mean)))
			tf.scalar_summary('sttdev/' + name, stddev)
			tf.scalar_summary('max/' + name, tf.reduce_max(var))
			tf.scalar_summary('min/' + name, tf.reduce_min(var))
			tf.histogram_summary(name, var)



