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

import argparse

def create_dataset():
    parser = argparse.ArgumentParser()
    output = parser.add_argument_group("Output")
    output.add_argument("--no_output", action="store_true", help="do not generate output")
    output.add_argument("--name", type=str, required=True, help="A distinguishable name for the dataset")
    output.add_argument("--datasets_path", type=str, required=False, help="Optional path for the datasets directory. Especially useful, if datasets should be written to second disk")
    output.add_argument("--grids", nargs="*", default=["pressure", "pressure_static", "pressure_dynamic", "vel", "phi_fluid", "density"], choices=["pressure", "pressure_static", "pressure_dynamic", "vel", "phi_fluid", "density"], help="Specify which grids should be written to uni files.")
    output.add_argument("--quantization", type=int, default=0, help="Decimal places for quantization of grids. 0 means off")

    scene = parser.add_argument_group("Scene")
    scene.add_argument("-n", "--num_scenes", type=int, default=1, help="number of scenes used in data generation")
    parser.add_argument('-t', "--type", default="flip", \
        choices=["flip", "smoke_buoyant", "smoke_simple"], help='simulation type (default: %(default)s)')
    scene.add_argument("-s", "--simulation_steps", type=int, default=100, help="number of simulation steps for a scene. must be constant for all scenes")
    scene.add_argument("-w", "--warmup", type=int, default=75, help="number of steps to discard in the beginning of the scene")
    scene.add_argument("--seed", type=int, default=1337, help="Seed for the random number generator")
    scene.add_argument("--skip_steps", type=int, default=1, help="Write out interval. Value of 1 writes every step")
    scene.add_argument("--obstacles", action="store_true", help="add obstacles to the scene")
    scene.add_argument("--meshes", action="store_true", help="add meshes as water surface scene")
    scene.add_argument("-rx", "--resolution_x", type=int, default=64, help="simulation resolution x")
    scene.add_argument("-ry", "--resolution_y", type=int, default=64, help="simulation resolution y")
    scene.add_argument("-rz", "--resolution_z", type=int, default=64, help="simulation resolution z")
    scene.add_argument("-d", "--dimension", type=int, default=3, help="simulation dimension (usually 2D or 3D)")

    mantaflow = parser.add_argument_group("Mantaflow")
    mantaflow.add_argument("-g", "--gui", action="store_true")

    debug = parser.add_argument_group("Debug")
    debug.add_argument("-p", "--pause", action="store_true", help="Pause on start")

    #args = parser.parse_args()
    args, unknown_args = parser.parse_known_args()

    def process_arg_name(arg_name):
        while arg_name[0] is "-":
            arg_name = arg_name[1:]
        return arg_name
    def unpack_list(arg_list):
        if type(arg_list) is list and len(arg_list) == 1:
            return arg_list[0]
        return arg_list

    # parse unknown arguments
    unknown_args_dict = {}
    last_split = 0
    for i, b in enumerate(unknown_args):
        if type(b) is str and (b.startswith("-") or b.startswith("--")):
            list_part = unknown_args[last_split:i]
            if list_part:
                # remove the leading "-" or "--"
                arg_name = process_arg_name(list_part[0])
                # if the list_part has only one value, it is seen as a bool var otherwise as list
                unknown_args_dict[arg_name] = unpack_list(list_part[1:] if len(list_part) > 1 else True)
            last_split = i

    # add last arg
    if unknown_args:
        arg_name = process_arg_name(unknown_args[last_split])
        unknown_args_dict[arg_name] = unpack_list(unknown_args[last_split+1:] if len(unknown_args[last_split:]) > 1 else True)

    #args=parser.parse_args()
    del output
    del scene
    del mantaflow
    del debug

    return args, unknown_args_dict

# ---------------------------------------------------------------------------------------------------------
def display_dataset():
    parser = argparse.ArgumentParser()
    dataset = parser.add_argument_group("Data Set")
    dataset.add_argument("--dataset", type=str, required=True, help="Path to already existing data set")
    dataset.add_argument("--datasets_path", type=str, required=False, help="Optional path for the datasets directory. Especially useful, if datasets should be read from second disk")
    dataset.add_argument("--grids", nargs="*", default=["pressure", "pressure_static", "pressure_dynamic", "vel", "phi_fluid", "density"], choices=["pressure", "pressure_static", "pressure_dynamic", "vel", "phi_fluid", "density"], help="specify the displayed grids.")

    scene = parser.add_argument_group("Scene")
    scene.add_argument("-s", "--start_scene", type=int, default=0, help="scene number to start")
    scene.add_argument("-n", "--num_scenes", type=int, default=-1, help="number of scenes to display beginning from start_scene. Default of -1 displays all scenes")
    scene.add_argument("-t", "--time_to_wait", type=float, default=0.0, help="time to wait in each iteration (additionally to IO)")
    scene.add_argument("--simulation_steps", type=int, default=-1, help="number of simulation steps to display for a scene")
    scene.add_argument("--skip_steps", type=int, default=1, help="display interval. Value of 1 displays all stored time steps")

    general = parser.add_argument_group("General")
    general.add_argument("--video", action="store_true", help="create a video of the dataset content")

    args = parser.parse_args()
    del dataset
    del scene
    del general

    return args
