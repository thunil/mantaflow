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
mac = s.create(MACGrid)
pp = s.create(BasicParticleSystem) 
flags.initDomain(boundaryWidth=1)

if 0 and GUI:
    gui = Gui()
    gui.show()
    #gui.pause()
    
# try to load uni file if it exists
def tryToLoad(grid, name):
    rfile = name 
    if(os.path.isfile(rfile)):
        grid.load(rfile)
        #print("Did load " + rfile)
        return 1
        
def numToString(num):
    if(num<10):
        string = '000'+str(num)
    elif(num<100):
        string = '00' + str(num)
    elif(num<1000):
        string = '0' + str(num)
    else:
        string = str(num)
    return string
    
def diffGrids(grid, parent, name2, i):
    if ( type(grid).__name__ == "MACGrid" ):
        gridTmpMac = parent.create(VecGrid)
        copyMacToVec3(grid , gridTmpMac )
        return diffGrids(gridTmpMac, parent, name2, i)
    if ( type(grid).__name__ == "LevelsetGrid" ):
        gridTmpLs = parent.create(RealGrid)
        copyLevelsetToReal(grid , gridTmpLs )
        return diffGrids(gridTmpMac, parent, name2, i)
        
    if ( grid._class == "Grid" and grid._T == "Real" ):
        compareTmpGrid = parent.create(RealGrid)
        tryToLoad(compareTmpGrid, name2)
        errVal = gridMaxDiff    ( grid, compareTmpGrid )
    elif ( grid._class == "Grid" and grid._T == "Vec3" ):
        compareTmpGrid = parent.create(VecGrid)
        tryToLoad(compareTmpGrid, name2)
        errVal = gridMaxDiffVec3( grid, compareTmpGrid )
    elif ( grid._class == "Grid" and grid._T == "int" ):
        compareTmpGrid = parent.create(IntGrid)
        tryToLoad(compareTmpGrid, name2)
        errVal = gridMaxDiffInt ( grid, compareTmpGrid )
    elif ( grid._class == "ParticleDataImpl"):
        if(grid._T == "Real" ):
            compareTmpGrid = parent.create(PdataReal)
        elif (grid._T == "Vec3" ):
            compareTmpGrid = parent.create(PdataVec3)
        elif (grid._T == "int" ):
            compareTmpGrid = parent.create(PdataInt)
        tryToLoad(compareTmpGrid, name2)
        errVal = pdataMaxDiff ( grid, compareTmpGrid )
    else:
        print( "Error doTestGrid - unknown grid type " + type(grid).__name__+ " class:"+grid._class+ " T:"+grid._T  )
        return 1
    
    
    maxVal = grid.getMaxAbsValue() + 1e-15
    
    if 1:
        minVal1 = grid.getMinValue()
        maxVal1 = grid.getMaxValue()
        minVal2 = compareTmpGrid.getMinValue()
        maxVal2 = compareTmpGrid.getMaxValue()
        if (minVal1==0 and maxVal1==0) or (minVal2==0 and maxVal2==0):
            print("timestep "+str(i)+" min/max curr "+str(minVal1)+" to "+str(maxVal1)+" min/max ref "+str(minVal2)+" to "+str(maxVal2) );

    errValRel = errVal/maxVal
    if (errValRel!=0):
        print("timestep "+str(i)+", errValRel: "+str(errValRel))
    
    #gridSub(grid,compareTmpGrid)
    #s.step()
    #gui.screenshot(path + file[:-4] + '.png')
    #compareTmpGrid.clear()
        
path = "Z:\\MantaflowMaster\\mantaflowgit\\build\\Release\\"
for i in range(100):
    names = ['vel1_'+numToString(i)+'.uni','vel2_'+numToString(i)+'.uni']
    if tryToLoad(mac, path + names[0]):
        diffGrids(mac,s,path+names[1],i)