import os,gzip
import struct
from collections import namedtuple
import glob,re
import numpy as np

# this is for reading the field
def _read_content(bytestream, head):
    #dt = np.dtype(np.float64).newbyteorder('>')
    assert (head['bytesPerElement'] == 12 and head['elementType'] == 2) or head['bytesPerElement'] == 4
    if (head['elementType'] == 0):
        data = np.frombuffer(bytestream.read(), dtype="int32")
    else:
        data = np.frombuffer(bytestream.read(), dtype="float32")
    if (head['elementType'] == 2):
        return data.reshape(head['dimX'], head['dimY'], head['dimZ'], 3, order='C')
    else:
        return data.reshape(head['dimX'], head['dimY'], head['dimZ'], order='C')

# read important information in the header such as dimensions and grid type
def _read_head(bytestream):
    ID = bytestream.read(4)
    #dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, timestamp = struct.unpack('iiiiii256sQ', bytestream.read(288))

    # unpack header struct object
    head = namedtuple('UniHeader', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, timestamp')
    # convert to namedtuple and then directly to a dict
    head = head._asdict(head._make(struct.unpack('iiiiii256sQ', bytestream.read(288))))
    #head['ID'] = ID
    
    return head

# use this to read the .uni file. It will return the header as dictionary and the content as np-array
def readuni(filename):
    with gzip.open(filename, 'rb') as bytestream:
        head = _read_head(bytestream)
        content = _read_content(bytestream, head)

        return head, content

# use this to write a .uni file. The head has to be supplied in the same dictionary format as the output of readuni
def writeuni(filename, head, content):
    with gzip.open(filename, 'wb') as bytestream:
        # write the head of the uni file
        bytestream.write(b'MNT2')
        head_tuple = namedtuple('GenericDict', head.keys())(**head)
        head_buffer = struct.pack('iiiiii256sQ', *head_tuple)
        bytestream.write(head_buffer)
        if (head['elementType'] == 2):
            content = content.reshape(head['dimX']*head['dimY']*head['dimZ']*3, order='C')
        else:
            content = content.reshape(head['dimX']*head['dimY']*head['dimZ'], order='C')

        bytestream.write(content.tobytes())

def convert_to_npz(dir, file_name, header, quantization=0):
    #TODO reactivate assert path.suffix == '.npz', ("Serialization must be to a .npz file")

    fields = []
    # gather file list of all files with ending .uni in dir
    dirlist = glob.glob(dir+"*.uni")

    # natural sort the array
    convert = lambda text: int(text) if text.isdigit() else text
    dirlist.sort(key=lambda k: [ convert(c) for c in re.split('([0-9]+)', k) ] )

    # load every file remaining in list
    for index, filename in enumerate(dirlist):
        # extract data from uni file (headers are ignored)
        field = readuni(filename)[1]
        # if 3d add 4th dimension
        if field.ndim == 3:
            field = field[...,np.newaxis]
        fields.append(field)
    fields = np.array(fields)
    if quantization != 0:
        fields = np.around(fields, decimals=quantization)

    np.savez_compressed(dir + file_name, data=fields, header=header)
    for f in glob.glob(dir + "*.uni"):
        os.remove(f)
