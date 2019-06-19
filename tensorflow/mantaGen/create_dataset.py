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
import json
import random
import datetime
from math import pi
import numpy as np

import scenes.volumes as v

from util.path import find_dir, make_dir, get_unique_path
from util import arguments
from util.io import GridIO
from util.logger import *

# Example Call
# ./build/manta create_dataset.py --name=TESTSIM -t smoke_buoyant -n 2 -s 10 -w 15 -d 2 -g --datasets_path="./datasets/"

# Arguments
#--------------------------------
args, unknown_args = arguments.create_dataset()

# Scene initialization
#--------------------------------
np.random.seed(seed=args.seed)
random.seed(args.seed)

step_file_num = 0

simulation_steps = args.simulation_steps
simulation_steps *= args.skip_steps
simulation_steps += args.warmup

if simulation_steps < 2 or args.skip_steps != 1:
    warning("It will not be possible to train the LSTM with this dataset")

# load scene class from string via module of the same name, e.g., "SmokeScene" loads "./scenes/SmokeScene.py"
#--------------------------------
import importlib
moduleName = "scenes."+args.type+"_scene"
info("Loading module "+moduleName)
sceneModule = importlib.import_module(moduleName)

scene = sceneModule.instantiate_scene(resolution=Vec3(args.resolution_x, args.resolution_y, args.resolution_z), boundary=1, dimension=args.dimension, timestep=0.5, gravity=Vec3(0,-0.003,0), merge_ghost_fluid=True, show_gui=args.gui, pause_on_start=args.pause, **unknown_args)

# create output folder
#--------------------------------
if not args.no_output:
    # get the datasets directory, in all datasets should reside
    if not args.datasets_path:
        dataset_path = find_dir("datasets", 1)
    else:
        dataset_path = args.datasets_path
    assert os.path.exists(dataset_path), ("Datasets directory {} does not exist".format(dataset_path)) 
    # set output parent directory name
    dataset_path += "/" + args.name
    # always append index
    dataset_path = get_unique_path(dataset_path)
    #warning("Dataset with name '{}' already exists. [{}]".format(args.name, dataset_path))
    

# check for grids to be saved as outputs
#--------------------------------
grids = []
for grid_name in args.grids:
    try:
        getattr(scene, grid_name)
    except AttributeError as err:
        warning("Scene has no grid named {}".format(grid_name))
    else:
        grids.append(grid_name)
args.grids = grids
info("Grids selected for output: {}".format(grids))

grid_io = GridIO()

# Simulate scene
#--------------------------------
def on_sim_step(scene, t):
    global step_file_num
    
    if args.no_output:
        return
    # screenshot of the mantaflow gui
    if t == args.warmup and scene.show_gui:
        scene._gui.screenshot(output_path + "/screenshots/screen_{:06d}.jpg".format(step_file_num))
    # warmup / skip steps should not be written to file
    if t < args.warmup or t % args.skip_steps != 0:
        return

    # write grids to a file
    step_file_num = t - args.warmup
    output_name = "{:06d}".format(step_file_num)

    for grid_name in grids:
        # it was already checked if the attribute is present in the scene
        grid = getattr(scene, grid_name)
        # save the grid to a npz file
        grid_io.write_grid(grid, output_path + "/" + grid_name + "_" + output_name)

# reset seed -> same starting condition for scene creation -> reproducibility
np.random.seed(seed=args.seed)
random.seed(args.seed)

for scene_num in range(args.num_scenes):
    # Output dirs, one per simulation
    #--------------------------------
    if not args.no_output:
        output_path = dataset_path + "/sim"
        output_path = get_unique_path(output_path)
        make_dir(output_path)
        if args.gui:
            make_dir(output_path + "/screenshots/")
        # keep track of stored scenes
        stored_scenes_num = 0

    # Dataset description
    #--------------------------------
    description = {}
    description["version"] = "v0.01"
    description["creation_date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    description["grids"] = grids
    description["resolution_x"] = args.resolution_x
    description["resolution_y"] = args.resolution_y
    description["resolution_z"] = args.resolution_z
    description.update(vars(args)) # insert args

    if not args.no_output:
        with open(output_path + "/description.json", 'w') as f:
            json.dump(description, f, indent=4)

    info("<><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    info("Scene {} ({})".format(scene_num, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    info("Output path: {}".format(output_path))
    # increase scene count
    stored_scenes_num = scene_num

    # [the following "scene_selection" was moved to _create_scene all of derived class]
    # scene_setup.scene_selection(args.type, scene, obstacles=args.obstacles, meshes=args.meshes)

    # restart the simulation with the new scene and write out the grids as .uni
    scene.simulate(num_steps=simulation_steps, on_simulation_step=on_sim_step)
