#******************************************************************************
#
# MantaGen
# Copyright 2018 Steffen Wiewel, Moritz Becher, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
#******************************************************************************

from manta import *
import numpy
from scenes.volumes import *

def random_vec3(vmin=[-1,-1,-1], vmax=[1,1,1]):
    return vec3(
        numpy.random.uniform(low=vmin[0], high=vmax[0]),
        numpy.random.uniform(low=vmin[1], high=vmax[1]),
        numpy.random.uniform(low=vmin[2], high=vmax[2])
    )

def random_vec3s(vmin=-1, vmax=1): # scalar params
    return vec3(
        numpy.random.uniform(low=vmin, high=vmax),
        numpy.random.uniform(low=vmin, high=vmax),
        numpy.random.uniform(low=vmin, high=vmax)
    )

def random_box(center_min=[0,0,0], center_max=[1,1,1], size_min=[0,0,0], size_max=[1,1,1], is3d=True):
    size = random_vec3(size_min, size_max)
    center = random_vec3(center_min, center_max)
    if not is3d:	
         size.z   = 1.0
         center.z = 0.5
    return BoxVolume(center=center, size=size)

def random_velocity(vel_min=[-1,-1,-1], vel_max=[1,0,1]):
    return random_vec3(vel_min, vel_max)
