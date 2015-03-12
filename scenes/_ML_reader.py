#
# Read in generated data, display & compare it
# Needs manual file paths and grid dimension settings for now...
# 
import sys
import os 
import glob
from array import *
from manta import *

# solver params
res = 64
gs = vec3(res,res,1)
dim = 2

# solver setup
s = Solver(name='main', gridSize = gs, dim=dim)
flags = s.create(FlagGrid)
real = s.create(RealGrid)
mac  = s.create(MACGrid)
#flags.initDomain(boundaryWidth=0)

if GUI:
    gui = Gui()
    gui.show()
    gui.pause()
    
# try to load uni file if it exists
def tryToLoad(grid, name, suffix):
    rfile = name + suffix 
    if(os.path.isfile(rfile)):
        grid.load(rfile)
        print("Did load " + rfile)
        return 1
        
path = 'Z:\\MantaflowMaster\\mantaflowgit\\build\\'
names = ['flagsApplyMatrix2D','flagsClosedBounds','flagsIterate1','flagsOpenBounds','flagsSolvePressure']
# iterate over time of screenshots
os.chdir(path)
for file in glob.glob('flags*.uni'):
    if (os.path.isfile(path+file)):
        if tryToLoad(flags, path + file[:-4], ".uni"):
            s.step()
            gui.screenshot(path + file[:-4] + '.png')
        flags.clear()
            