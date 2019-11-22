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

import os

def find_dir(dirname, parent_levels=0):
    """ 
    find_dir searches for a directory in the proximity of the working directory with the specified name
    It will always search in the subdirectories of the working directory, but can be allowed to search in
    parent directories as well.
    * __dirname__: the name of the directory to search for
    * __parent_levels__: a value of 0 to disable search in parent directories, > 0 max levels of steps upwards in the directory tree
    """
    dirs_to_search = ["./"]
    for i in range(1, parent_levels + 1):
        dirs_to_search.append("../" * i)
    for start_dir in map(os.path.abspath, dirs_to_search):
        for dirpath, _, _ in os.walk(start_dir): 
            if dirpath.split(os.path.sep)[-1] == dirname:
                return dirpath
    raise RuntimeError("Could not find directory '{}'. Check if programm is executed in the right place.".format(dirname))

def make_dir(directory):
    """ 
    check if directory exists, otherwise makedir
    
    (workaround for python2 incompatibility)
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_unique_path(path):
    """
    create directory with unique name

    if directory already exists, count up until unique name
    """
    output_num = 0
    unique_path = path + "_{:06d}".format(output_num)
    while os.path.exists(unique_path):
        output_num += 1
        unique_path = path + "_{:06d}".format(output_num)
    return unique_path
