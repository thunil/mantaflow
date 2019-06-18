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

from __future__ import print_function

import os
import numpy as np
import json
from random import randint, getrandbits
import datetime
from math import pi

import scenes.volumes as v
from util.path import find_dir, make_dir, get_unique_path
from util import git
from util import arguments
from enum import Enum
from manta import *


#--------------------------------
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

#--------------------------------
def plot_field(x, out_path, title="", vmin=None, vmax=None, plot_colorbar=True):
    #print("plot_field", x.shape)
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
    cmap = plt.cm.get_cmap('viridis')
    colors = cmap(np.arange(cmap.N))

    max_color = 10
    min_color = [1,1,1,1]
    final_color = colors[max_color]
    for i in range(max_color):
        for j in range(4):
            interpolant = (max_color - i) / max_color
            colors[i][j] =  min(final_color[j] * (1.0 - interpolant) + min_color[j] * interpolant, 1.0)
    cm = LinearSegmentedColormap.from_list("custom_vir", colors, N=256)

    im1 = ax1.imshow(x[:,:], vmin=vmin, vmax=vmax, cmap=cm, interpolation='nearest')

    # contours
    cmap = plt.cm.get_cmap('viridis')
    colors = cmap(np.arange(cmap.N))
    for i in range(len(colors)):
        for j in range(4):
            colors[i][j] =  min(colors[i][j] * 1.4, 1.0)
    cm = LinearSegmentedColormap.from_list("custom_vir", colors, N=256)
    #cm = plt.get_cmap('autumn')

    cs = ax1.contour(x[:,:], levels=np.arange(vmin+(vmax-vmin)*0.025,vmax,(vmax-vmin)/10.0), cmap=cm)

    #if plot_colorbar:
        #divider = make_axes_locatable(ax1)
        #cax = divider.append_axes('right', size='5%', pad = 0.05)
        #fig.colorbar(im1, cax=cax, orientation='vertical')
    ax1.set_xlim(0, x.shape[1])
    ax1.set_ylim(0, x.shape[0])
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    if title:
        fig.suptitle(title, fontsize=20)
    if plot_colorbar:
        fig.savefig(out_path, bbox_inches='tight')
    else:
        fig.savefig(out_path, bbox_inches='tight', transparent=True, pad_inches=0)
    plt.close(fig)

# Enum
#--------------------------------
# mirror class of GridBase::GridType -> keep in sync
class MantaGridType(Enum):
    TypeNone = 0
    TypeReal = 1
    TypeInt = 2
    TypeVec3 = 4
    TypeMAC = 8
    TypeLevelset = 16
    TypeFlags = 32
    # Manually added permutations
    TypeMACVec3 = 12
    TypeLevelsetReal = 17
    TypeFlagsInt = 34

def get_manta_type(enum_type):
    if enum_type is MantaGridType.TypeReal:
        return RealGrid
    if enum_type is MantaGridType.TypeLevelset or enum_type is MantaGridType.TypeLevelsetReal:
        return LevelsetGrid
    if enum_type is MantaGridType.TypeVec3:
        return MACGrid
    if enum_type is MantaGridType.TypeMAC or enum_type is MantaGridType.TypeMACVec3:
        return MACGrid
    assert False, "Not supported"


# Meshes
#--------------------------------
mesh_types=["duck.obj", "bunny.obj", "simpletorus.obj", "anvil.obj"]

def vec_to_str(vec):
    return "{}, {}, {}".format(vec.x, vec.y, vec.z)

def spawn_mesh(mesh_dir, scene):
    # use random mesh from meshes directory with random position, rotation and scaling
    mesh = mesh_dir + "/" + mesh_types[randint( 0, len(mesh_types)-1 )]
    center = v.random_vec3(low=[0.05, 0.5, 0.05], high=[0.95, 0.95, 0.95])
    scale = np.random.uniform(low=0.4, high=0.9)
    rotation = v.random_vec3(low=[-pi, -pi, -pi], high=[pi, pi, pi])
    print("Spawning '{}': Pos {} Rot {} Scale {}".format(mesh, vec_to_str(center), vec_to_str(rotation), scale ))
    mesh = v.MeshVolume(
        mesh,
        center=center,
        scale=scale,
        rotation=rotation
    )
    scene.add_fluid(mesh)
    scene.set_velocity(mesh, v.random_velocity(vel_min=[-3.5, -3.5, -3.5], vel_max=[3.5, 0, 3.5]))


#----------------------------------------------------------------------------------
# Drop (Box)
def spawn_drop(scene):
    box = v.random_box(center_min=[0, 0.7, 0], center_max=[1, 0.9, 1], size_min=[0.005, 0.005, 0.005], size_max=[0.25, 0.25, 0.25])
    scene.add_fluid(box)
    drop_velo = v.random_velocity(vel_min=[-2.5, -2.5, -2.5], vel_max=[2.5, 0, 2.5])
    scene.set_velocity(box, drop_velo)
    # print("Drop {}:\n\t{},{},{}\n\t{},{},{}\n\t{},{},{}".format(i, box._center.x, box._center.y, box._center.z,
    #                                                         box._size.x, box._size.y, box._size.z,
    #                                                         drop_velo.x, drop_velo.y, drop_velo.z))

#----------------------------------------------------------------------------------
# test
def initialize_simple_scene(scene):
    source_count = randint(1, 2)
    for i in range(source_count):
        box = v.random_box(center_min=[0.2, 0.1, 0.2], center_max=[0.8, 0.6, 0.8], size_min=[0.005, 0.005, 0.005], size_max=[0.2, 0.2, 0.2])
        # TODO scene.add_source( box )

#----------------------------------------------------------------------------------
def initialize_smoke_scene(scene):
    source_count = randint(4, 10)
    for i in range(source_count):
        box = v.random_box(center_min=[0.2, 0.1, 0.2], center_max=[0.8, 0.6, 0.8], size_min=[0.005, 0.005, 0.005], size_max=[0.2, 0.2, 0.2])
        scene.add_source( box )
        #scene.set_velocity(box, random_vec3(vel_min=[-2.5, -2.5, -2.5], vel_max=[2.5, 2.5, 2.5]))
        #scene.add_source( box.shape(scene.solver) )

#----------------------------------------------------------------------------------
def initialize_liquid_scene(scene, simple=False, obstacles=False, meshes=False, mesh_path="meshes"):
    if simple:
        # low basin is the large volume of liquid at the bottom of the scene. in general adds most of the liquid in the scene
        low_basin = v.random_box(center_min=[0.5, 0.0, 0.5], center_max=[0.5, 0.0, 0.5], size_min=[1, 0.3, 1], size_max=[1, 0.5, 1])
        scene.add_fluid(low_basin)
        scene.set_velocity(low_basin, v.random_velocity(vel_min=[-0.5, 0, -0.5], vel_max=[0.5, 0, 0.5]))
        # the high basin is a tall but narrow block of liquid, adding motion to the scene
        high_basin = v.random_box(center_min=[0.1, 0.3, 0.1], center_max=[0.9, 0.4, 0.9], size_min=[0.1, 0.2, 0.1], size_max=[0.3, 0.4, 0.3])
        scene.add_fluid(high_basin)
        scene.set_velocity(high_basin, v.random_velocity(vel_min=[-2.0, -2.0, -2.0], vel_max=[2.0, 0, 2.0]))
    else:
        # low basin is the large volume of liquid at the bottom of the scene. in general adds most of the liquid in the scene
        low_basin = v.random_box(center_min=[0.5, 0.0, 0.5], center_max=[0.5, 0.0, 0.5], size_min=[1, 0.3, 1], size_max=[1, 0.6, 1])
        scene.add_fluid(low_basin)
        low_basin_velo = v.random_velocity(vel_min=[-1.5, 0, -1.5], vel_max=[1.5, 0, 1.5])
        scene.set_velocity(low_basin, low_basin_velo)
        # print("Low:\n\t{},{},{}\n\t{},{},{}\n\t{},{},{}".format(low_basin._center.x, low_basin._center.y, low_basin._center.z,
        #                                                         low_basin._size.x, low_basin._size.y, low_basin._size.z,
        #                                                         low_basin_velo.x, low_basin_velo.y, low_basin_velo.z))

        # the high basin is a tall but narrow block of liquid, adding motion to the scene
        high_basin = v.random_box(center_min=[0.0, 0.3, 0.0], center_max=[1.0, .3, 1.0], size_min=[0.1, 0.2, 0.1], size_max=[0.4, 0.5, 0.4])
        scene.add_fluid(high_basin)
        high_basin_velo =  v.random_velocity(vel_min=[-2.5, -2.5, -2.5], vel_max=[2.5, 0, 2.5])
        scene.set_velocity(high_basin, high_basin_velo)
        # print("High:\n\t{},{},{}\n\t{},{},{}\n\t{},{},{}".format(high_basin._center.x, high_basin._center.y, high_basin._center.z,
        #                                                         high_basin._size.x, high_basin._size.y, high_basin._size.z,
        #                                                         high_basin_velo.x, high_basin_velo.y, high_basin_velo.z))

        # drops: spawn at most 3 drops or meshes
        drop_count = randint(0, 3)
        for i in range(drop_count):
            if meshes and bool(getrandbits(1)):
                spawn_mesh(mesh_path, scene)
            else:
                spawn_drop(scene)

    # optional: add solid boxes to scene. more realistic, harder to learn
    if obstacles:
        obstacle = v.random_box(center_min=[0.1,0.1,0.1], center_max=[1.0,0.5,1.0], size_min=[0.2,0.2,0.2], size_max=[0.3,1.0,0.3])
        print("Added obstacle to the scene")
        scene.add_obstacle(obstacle)

#----------------------------------------------------------------------------------
def scene_selection(scene_type, scene, obstacles=False, meshes=False):
	print("Scene Type: {}".format(scene_type))
	if scene_type == "liquid":
		initialize_liquid_scene(scene, simple=False, obstacles=obstacles)
	elif scene_type == "smoke":
		initialize_smoke_scene(scene)
	elif scene_type == "simple":
		initialize_simple_scene(scene)
	else:
		assert False, "Unknown scene type " + scene_type

