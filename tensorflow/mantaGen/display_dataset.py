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
from random import randint, getrandbits, seed
import datetime
from enum import Enum
import subprocess

from scenes.display_scene import DisplayScene
from scene_setup import MantaGridType
from util import uniio
from util.path import find_dir, make_dir, get_unique_path
from util import git
from util import arguments

# Functions
#--------------------------------
def load_description(directory):
    desc = None
    if os.path.isfile(directory + "description.json"):
        with open(directory + "description.json", 'r') as f:
            desc = json.load(f)
    assert desc is not None, ("Description '" + directory + "description.json" + "' not found")
    return desc

# Arguments
#--------------------------------
args = arguments.display_dataset()

# Data Set
#--------------------------------
if not args.datasets_path:
    dataset_path = find_dir("datasets", 2)
else:
    dataset_path = args.datasets_path
assert os.path.exists(dataset_path), ("Datasets directory {} does not exist".format(dataset_path)) 
# set output parent directory name
dataset_path += "/" + args.dataset
dataset_path = dataset_path + "/"

if args.video:
    dataset_video_path = dataset_path+"dataset_display/"
    make_dir(dataset_video_path)

print("Reading data set from: {}".format(dataset_path))

# Description
#--------------------------------
dataset_desc = load_description(dataset_path)

# Grids
#--------------------------------
dataset_grids = []
for grid in args.grids:
    for grid2 in dataset_desc["grids"]:
        if grid == grid2:
            dataset_grids.append(grid)
            break
print("Found grids: {}".format(dataset_grids))

# Scenes
#--------------------------------
dataset_start_scene = args.start_scene
assert 0 <= dataset_start_scene < dataset_desc["num_scenes"], "start_scene {} is beyond data set range [0,{}]".format(dataset_start_scene, dataset_desc.num_scenes)
dataset_end_scene = dataset_desc["num_scenes"] if args.num_scenes < 0 else args.num_scenes
dataset_end_scene = min(dataset_start_scene + dataset_end_scene, dataset_desc["num_scenes"])
assert dataset_start_scene < dataset_end_scene, "start_scene {} must be smaller than end_scene {}".format(dataset_start_scene, dataset_end_scene)
print("Scenes to display: {} -> {}".format(dataset_start_scene, dataset_end_scene))

dataset_simulation_steps = dataset_desc["simulation_steps"] if args.simulation_steps < 0 else args.simulation_steps
dataset_simulation_steps = min(dataset_desc["simulation_steps"], dataset_simulation_steps)

# Scene initialization
#--------------------------------
scene = DisplayScene(resolution=dataset_desc["resolution"], dimension=dataset_desc.get("dimension", 3), time_to_wait=args.time_to_wait, skip_steps=args.skip_steps)

# check scene for errors
grids = []
for grid_name in dataset_grids:
    try:
        getattr(scene, grid_name)
        grids.append(grid_name)
    except AttributeError as err:
        print("Scene has no grid named {}".format(grid_name))
dataset_grids = grids
print("Remaining grids: {}".format(dataset_grids))

# Scene simulate
#--------------------------------
def on_grid_copy_step(scene, t):
    # iterate grids and transfer data
    for name, content in scene.display_grids.items():
        scene_grid = getattr(scene, name)
        grid_type = MantaGridType(scene_grid.getGridType())
        if grid_type == MantaGridType.TypeReal:
            copyArrayToGridReal(content[t], scene_grid)
        elif grid_type == MantaGridType.TypeInt:
            assert False, "Not supported"
        elif grid_type == MantaGridType.TypeVec3:
            copyArrayToGridVec3(content[t], scene_grid)
        elif grid_type == MantaGridType.TypeMAC or grid_type == MantaGridType.TypeMACVec3:
            copyArrayToGridMAC(content[t], scene_grid)
        elif grid_type == MantaGridType.TypeLevelset or grid_type == MantaGridType.TypeLevelsetReal:
            copyArrayToGridReal(content[t], scene_grid)

#--------------------------------
def on_solver_step(scene, t):
    # screenshot of the mantaflow gui
    if scene.show_gui and args.video:
        scene._gui.screenshot(dataset_path + "/dataset_display/screen_{:06d}.jpg".format(scene.file_num))

#--------------------------------
# Display Loop
for scene_num in range(dataset_start_scene, dataset_end_scene):
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    print("Scene {} ({})\n".format(scene_num, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    # load npz data
    for grid_name in dataset_grids:
        cur_npz = np.load(dataset_path + grid_name + "/{:06d}.npz".format(scene_num))
        scene.display_grids[grid_name] = cur_npz["data"]
        cur_npz.close()

    # restart the simulation with the new scene and write out the grids as .uni
    scene.simulate(num_steps=dataset_simulation_steps, on_grid_copy_step=on_grid_copy_step, on_solver_step=on_solver_step)

    # cleanup loaded grids
    scene.display_grids.clear()


if args.video:
    subprocess.call(['ffmpeg',
    '-r', '30',
    '-f', 'image2',
    '-start_number', '0',
    '-i', dataset_video_path + 'screen_%06d.jpg',
    '-vcodec', 'libx264',
    '-crf', '18',
    '-pix_fmt', 'yuv420p',
    '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
    dataset_video_path+"dataset.mp4"])

    files_in_dir = os.listdir(dataset_video_path)
    for item in files_in_dir:
       if item.endswith(".jpg"):
           os.remove(os.path.join(dataset_video_path, item))