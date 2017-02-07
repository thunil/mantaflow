#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2017 Nils Thuerey, Boris Bonev
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL) 
# http://www.gnu.org/licenses
#
# Read mantaflow uni files into numpy arrays
# note - only supports 3D grids for now
# (python2 , switch to python3 below)
#
#******************************************************************************

import gzip
import struct
import sys
from datetime import date
from collections import namedtuple
import numpy as np


# read content of grid
def RU_read_content(bytestream, head):
	assert (head['bytesPerElement'] == 12 and head['elementType'] == 2) or (head['bytesPerElement'] == 4 and (head['elementType'] == 0 or head['elementType'] == 1))

	if (head['elementType'] == 0):
		data = np.frombuffer(bytestream.read(), dtype="int32") # int grid
	else:
		data = np.frombuffer(bytestream.read(), dtype="float32") # float grid , scalar or vec3
	
	if (head['elementType'] == 2):
		return data.reshape(head['dimX'], head['dimY'], head['dimZ'], 3, order='C')
	else:
		return data.reshape(head['dimX'], head['dimY'], head['dimZ'], order='C')

# read uni file header (v3)
def RU_read_header(bytestream):
	ID = bytestream.read(4)

	if ID=="MNT2":
		# unpack header struct object
		head = namedtuple('HeaderV3', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, timestamp')
		# convert to namedtuple and then directly to a dict
		head = head._asdict(head._make(struct.unpack('iiiiii256sQ', bytestream.read(288))))

		# when writing, we'll need a v4 header field, re-pack...
		head['dimT'] = 0 
		head['info'] = head['info'][0:252]
		head4 = namedtuple('HeaderV4', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, dimT, timestamp')(**head)
		head = head4._asdict()

	elif ID=="MNT3":
		# unpack header struct object
		head = namedtuple('HeaderV4', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, dimT, timestamp')
		# convert to namedtuple and then directly to a dict
		# header is shorter for v3!
		head = head._asdict(head._make(struct.unpack('iiiiii252siQ', bytestream.read(288))))

	elif ID=="M4T2" or ID=="M4T3":
		print("read_header error - 4D grids not yet supported")
		exit(1)

	else:
		print("read_header error - unknown header '%s' " % ID)
		exit(1)
	
	return head

# use this to read the .uni file. It will return the header as dictionary and the content as np-array
def readUni(filename):
	#print("Reading '%s'" % filename) # debug
	with gzip.open(filename, 'rb') as bytestream:
		head = RU_read_header(bytestream)
		content = RU_read_content(bytestream, head)

		return head, content

# use this to write a .uni file. The head has to be supplied in the same dictionary format as the output of readuni
def writeUni(filename, head, content):
	#print("Writing '%s'" % filename) # debug
	with gzip.open(filename, 'wb') as bytestream:
		# write the header of the uni file
		#bytestream.write(b'MNT2') # v3
		#head_tuple = namedtuple('GenericDict', head.keys())(**head)
		#head_buffer = struct.pack('iiiiii256sQ', *head_tuple)
		bytestream.write(b'MNT3') # new, v4
		head_tuple = namedtuple('HeaderV4', head.keys())(**head)
		head_buffer = struct.pack('iiiiii252siQ', *head_tuple)
		bytestream.write(head_buffer)
		# write grid content
		if (head['elementType'] == 2):
			# vec3 grid
			content = content.reshape(head['dimX']*head['dimY']*head['dimZ']*3, order='C')
		else:
			# int or scalar grid
			content = content.reshape(head['dimX']*head['dimY']*head['dimZ'], order='C')

		if sys.version_info >= (3,0):
			# changed for Python3
			bytestream.write(memoryview(content))
		else:
			bytestream.write(np.getbuffer(content))


