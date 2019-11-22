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

#!/usr/bin/python3

from subprocess import call, check_output
import argparse
import shutil
import glob
import numpy as np
import scipy as sp
import re

import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir+"/src/")



def make_dir(directory):
    """ 
    check if directory exists, otherwise makedir
    
    (workaround for python2 incompatibility)
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_path", type=str, default="./uni/", help="Input file path")
parser.add_argument("-if", "--input_flip_file_name", type=str, default="flipParts_{:04d}.uni", help="Input flip file name")
parser.add_argument("-il", "--input_levelset_file_name", type=str, default="levelset_{:04d}.uni", help="Input levelset file name")
parser.add_argument("-o", "--output_path", type=str, default="./obj/", help="Output file path ")
parser.add_argument("-of", "--output_file_name", type=str, default="fluidsurface_final_{:04d}.bobj.gz", help="Output file name")
parser.add_argument("-mb", "--motion_blur", type=int, default=0, help="Motion blur for particles to reduce surface generation flickering")

args = parser.parse_args()

#----------------------------------------------------------------------------------
def main():
    # input file 
    partfile        = args.input_path + args.input_flip_file_name
    levelsetfile    = args.input_path + args.input_levelset_file_name
    interval        = 1

    # how much larger?
    upres = 2.0

    # create output dir
    make_dir(args.output_path)
    print("Writing output to '{}'".format(args.output_path))
    # output file name so that blender can directly read it...
    meshfile = args.output_path + args.output_file_name

    # resolution for level set / output mesh
    refName = args.input_path + "ref_" + args.input_flip_file_name.format(0)
    gs = getUniFileSize(refName)
    if gs.x<=0: 
        mantaMsg("Warning! File '%s' not found, cannot determine size...\n"%refName, 0)
        exit(1)

    # low res solver
    s_LR = Solver(name='lowres', gridSize = gs, dim=3)

    # high res solver
    gs.x = int(gs.x*upres)
    gs.y = int(gs.y*upres)
    gs.z = int(gs.z*upres)
    s = Solver(name='main', gridSize = gs , dim=3)

    # kernel radius for surface creation
    radiusFactor = 2.0

    # triangle scale relative to cell size
    #scale = 0.5

    # counters
    outCnt = 0
    frame = 0

    # prepare grids and particles
    flags    = s.create(FlagGrid)
    phi      = s.create(LevelsetGrid)
    phiParts = s.create(LevelsetGrid)
    pp       = s.create(BasicParticleSystem) 
    ppMb     = s.create(BasicParticleSystem) # optional 
    mesh     = s.create(Mesh)
    # low res
    phi_LR   = s_LR.create(LevelsetGrid)

    # acceleration data for particle nbs
    pindex = s.create(ParticleIndexSystem) 
    gpi    = s.create(IntGrid)

    # scene setup
    flags.initDomain(boundaryWidth=0)

    # main loop
    endFrame = len(glob.glob( os.path.dirname(args.input_path)+"/*.uni"))

    while frame < endFrame:
        meshfileCurr = meshfile.format(outCnt)
        #phi.setBound(value=0.5, boundaryWidth=1)

        # generate mesh; first read input sim particles
        if os.path.isfile( partfile.format(frame) ) and os.path.isfile( levelsetfile.format(frame) ):
            pp.load( partfile.format(frame) )
            if args.motion_blur>0:
                for mb in range(args.motion_blur):
                    frameMb = frame-mb if frame-mb>=0 else 0
                    ppMb.load( partfile.format(frameMb) )
                    pp.appendParticles(ppMb)
            # load low res phi from file and upscale to main solver
            phi_LR.load( levelsetfile.format(frame) )
            interpolateGrid(phi, phi_LR)
            flags.updateFromLevelset(phi)
            phi.reinitMarching(flags, 4.0)

            # create surface
            gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
            #unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor ) # faster, but not as smooth
            #averagedParticleLevelset( pp, pindex, flags, gpi, phiParts , radiusFactor , 2, 2 ) 
            improvedParticleLevelset( pp, pindex, flags, gpi, phiParts , radiusFactor , 1, 1 )  
            #phi.setBound(value=0.5, boundaryWidth=1)

            # Merge levelset and particles
            phi.addConst(1.); # shrink slightly
            phi.join( phiParts )
            # extrapolateLsSimple(phi=phi, distance=narrowBand+2, inside=True ) 
            # extrapolateLsSimple(phi=phi, distance=3 )
            # phi.setBoundNeumann(1) # make sure no particles are placed at outer boundary
            phi.setBound(0.0,3)

            phi.createMesh(mesh)
            # beautify mesh, too slow right now!
            #subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=3*scale, cutTubes=False) 
            # perform smoothing
            #for iters in range(10):
                #smoothMesh(mesh=mesh, strength=1e-3, steps=10) 
                #subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=3*scale, cutTubes=True)
            # write output file:
            mesh.save( meshfileCurr )

        outCnt += 1
        frame  += interval
        s.step()

#----------------------------------------------------------------------------------
# execute the prediction
if __name__ == "__main__":
    main()
