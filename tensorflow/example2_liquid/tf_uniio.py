# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL)
# http://www.gnu.org/licenses
#
# .uni files IO
#
# ----------------------------------------------------------------------------

from __future__ import print_function # for calling print(..., end='') in python2

import gzip, struct
import numpy as np
from collections import namedtuple

def _read_particle_data(bytestream, head, data_type=None): # data_type = {None: BasicParticleSystem; "float32": Real; "int32": Int}
    assert(head['bytesPerElement']==16 or head['bytesPerElement']==12 or head['bytesPerElement']==4)

    if(head['elementType']==0): # BasicParticleSystem
        print('(BasicParticleSystem) ', end='')
        data = np.frombuffer(bytestream.read(), dtype=np.dtype([('f1',(np.float32,3)),('f2',(np.int32,1))]))['f1']
    else:                       # head['elementType']==1: ParticleDataImpl<T>, where T = {float32: Real(4) or Vec3(12); int32: Int(4)}
        print('(ParticleDataImpl<T={}{}>) '.format(data_type, 'x3' if (head['bytesPerElement']==12) else ''), end='')
        data = np.reshape(np.frombuffer(bytestream.read(), dtype=data_type), (-1, 3 if (head['bytesPerElement']==12) else 1))

    return data

def _read_grid_data(bytestream, head, data_type=None):
    assert(head['bytesPerElement']==12 or head['bytesPerElement']==4)
    print('(Grid<T={}{}>) '.format(data_type, 'x3' if (head['bytesPerElement']==12) else ''), end='')
    data = np.frombuffer(bytestream.read(), dtype=data_type)
    if head['bytesPerElement']==12:
        return data.reshape((head['dimX'], head['dimY'], head['dimZ'], 3))
    else:
        return data.reshape((head['dimX'], head['dimY'], head['dimZ']))

def _read_particle_head(bytestream):
    ID = bytestream.read(4)     # NOTE: useless
    # unpack header struct object
    head = namedtuple('UniPartHeader', 'dim, dimX, dimY, dimZ, elementType, bytesPerElement, info, timestamp')
    # convert to namedtuple and then directly to a dict
    head = head._asdict(head._make(struct.unpack('iiiiii256sQ', bytestream.read(288))))

    return head

def _read_grid_head(bytestream):
    ID = bytestream.read(4)
    # unpack header struct object
    head = namedtuple('UniHeader', 'dimX, dimY, dimZ, gridType, elementType, bytesPerElement, info, dimT, timestamp')
    # convert to namedtuple and then directly to a dict
    head = head._asdict(head._make(struct.unpack('iiiiii252siQ', bytestream.read(288))))

    return head

# use this to read the .uni file. It will return the header as dictionary and the content as a numpy array
def read_particles(filename, data_type=None):
    print('Reading {} ... '.format(filename), end='')
    with gzip.open(filename, 'rb') as bytestream:
        head = _read_particle_head(bytestream)
        data = _read_particle_data(bytestream, head, data_type)

        print('Done.')
        return head, data

def read_grid(filename, data_type=None):
    print('Reading {} ... '.format(filename), end='')
    with gzip.open(filename, 'rb') as bytestream:
        head = _read_grid_head(bytestream)
        data = _read_grid_data(bytestream, head, data_type)

        print('Done.')
        return head, data

def drop_zdim(data):
    return np.delete(data, -1, 1)
