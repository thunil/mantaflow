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

from manta import *
import numpy as np
import os.path
from util.path import find_dir, make_dir, get_unique_path

# IO
#--------------------------------
class GridIO(object):
    """ Store manta grids as npz """
    #----------------------------------------------------------------------------------
    def __init__(self, **kwargs):
        self._buffer_list = []

    #----------------------------------------------------------------------------------
    def _print_grid_info(self, grid):
        # e.g. grid_name  phi_fluid                                 density
        # grid            <LevelsetGrid object at 0x7f3c32e6a0d0>   <Grid_Real object at 0x7f3c32e6a130>
        # grid.__class__  <class 'LevelsetGrid'>                    <class 'Grid_Real'>
        # grid._class     LevelsetGrid                              Grid
        # grid._cname     LevelsetGrid                              Grid<Real>
        # grid._T         --EMPTY--                                 Real
        # grid.getSize()  [+64,000000,+64,000000,+1,000000]         [+64,000000,+64,000000,+1,000000]
        # grid.is3D()     False                                     False
        print(grid)
        print(grid.__class__)
        print(grid._class)
        print(grid._cname)
        print(grid._T)
        print(grid.getSize())
        print(grid.is3D())
        #print(inspect.getmembers(grid))

    #----------------------------------------------------------------------------------
    def _get_grid_description(self, grid):
        grid_desc = {}
        grid_desc["__class__"]  = str(grid.__class__)
        grid_desc["__doc__"]    = str(grid.__doc__)
        grid_desc["_class"]     = str(grid._class)
        grid_desc["_cname"]     = str(grid._cname)
        grid_desc["_T"]         = str(grid._T)
        grid_desc["size"]       = str(grid.getSize())
        grid_desc["is3D"]       = grid.is3D()
        grid_desc["is4D"]       = grid.is4D()
        return grid_desc

    #----------------------------------------------------------------------------------
    def _get_buffer_handle(self, grid):
        #name = grid._cname
        size = grid.getSize()

        dimension = 0
        data_type = np.float32
        if grid._class == "Grid":
            if grid._T == "Real":
                dimension = 1
                data_type = np.float32
            elif grid._T == "int":
                dimension = 1
                data_type = np.int32
            elif grid._T == "Vec3":
                dimension = 3
                data_type = np.float32
            elif grid._T == "Vec4":
                dimension = 4
                data_type = np.float32
        elif grid._class == "FlagGrid":
            dimension = 1
            data_type = np.int32
        elif grid._class == "LevelsetGrid":
            dimension = 1
            data_type = np.float32
        elif grid._class == "MACGrid":
            dimension = 3
            data_type = np.float32
        else:
            assert False, "grid._class {} is not supported".format(grid._class)
            return

        # search for buffer by comparing shapes (Numpy Format: z,y,x,d | Manta Format: x,y,z,d)
        search_shape = (int(size.z), int(size.y), int(size.x), dimension) # resulting tuple for vel, e.g.: (64,64,128,3)

        for i in range(len(self._buffer_list)):
            if search_shape == self._buffer_list[i].shape and data_type == self._buffer_list[i].dtype:
                # return if exists
                return i

        # create new, add to list and return
        new_buffer = np.zeros(shape=search_shape, order='C', dtype=data_type)
        self._buffer_list.append(new_buffer)
        return len(self._buffer_list) - 1

    #----------------------------------------------------------------------------------
    # Store grids as npz
    def write_grid(self, grid, path):
        # check for support
        if grid._class != "Grid" and grid._class != "FlagGrid" and grid._class != "LevelsetGrid" and grid._class != "MACGrid":
            assert False, "grid._class {} is not supported".format(grid._class)
            return
        # check if path exists
        make_dir(os.path.dirname(path))
        # find buffer handle
        handle = self._get_buffer_handle(grid)
        vec3content = False
        # copy grid to buffer
        if grid._class == "Grid":
            if grid._T == "Real":
                copyGridToArrayReal(grid, self._buffer_list[handle])
            elif grid._T == "Vec3":
                copyGridToArrayVec3(grid, self._buffer_list[handle])
                vec3content = True
        elif grid._class == "FlagGrid":
            copyGridToArrayFlag(grid, self._buffer_list[handle])
        elif grid._class == "LevelsetGrid":
            copyGridToArrayLevelset(grid, self._buffer_list[handle])
        elif grid._class == "MACGrid":
            copyGridToArrayMAC(grid, self._buffer_list[handle])
            vec3content = True
        else:
            print("Grid is not supported")
            self._print_grid_info(grid)
            assert False

        # store buffer in npz
        grid_desc = self._get_grid_description(grid)
        np_grid = self._buffer_list[handle]
        is2d = (np_grid.shape[0] == 1)
        if vec3content and is2d:
            np_grid = np_grid[...,0:2] # remove z component of vectors
        if is2d:
            np_grid = np.squeeze(np_grid, axis=0) # remove z axis

        np.savez_compressed(path, data=np_grid, header=grid_desc)
