/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Main class for the fluid solver
 *
 ******************************************************************************/

#ifndef _FLUIDSOLVER_H
#define _FLUIDSOLVER_H

#include "pclass.h"
#include "vectorbase.h"
#include <vector>
#include <map>

namespace Manta {
    
//! Encodes grid size, timstep etc.
PYTHON(name=Solver) 
class FluidSolver : public PbClass {
public:
    PYTHON FluidSolver(Vec3i gridSize, int dim=3);
    virtual ~FluidSolver();
    
    // accessors
    PYTHON Vec3i getGridSize() { return mGridSize; }
    inline Real getDt() { return mDt; }
    inline Real getTime() { return mTimeTotal; }
    inline Real getDx() { return 1.0 / mGridSize.max(); }
    inline Real getScale() { return mScale; }
    //! Check dimensionality
    inline bool is2D() const { return mDim==2; }
    //! Check dimensionality
    inline bool is3D() const { return mDim==3; }
    
    // Python callable methods    
    //! output performace statistics
    PYTHON void printTimings();
    PYTHON void saveMeanTimings(std::string filename);
    
    //! Advance the solver one timestep, update GUI if present
    PYTHON void step();
    
    //! create a object with the solver as its parent
    PYTHON PbClass* create(PbType type, const std::string& name = "");
    
    // temp grid and plugin stuff: you shouldn't call this manually
    template<class T> T* getGridPointer();
    template<class T> void freeGridPointer(T* ptr);    
    void pluginStart(const std::string& name);
    void pluginStop(const std::string& name);        
protected:
    //! subclass for managing grid memory
    //! stored as a stack to allow fast allocation
    template<class T> struct GridStorage {
        GridStorage() : used(0) {}
        T* get(Vec3i size);
        void free();
        void release(T* ptr);
        
        std::vector<T*> grids;
        size_t used;        
    };
    
    Vec3i mGridSize;
    const int mDim;
    PYTHON(name=timestep) Real mDt;
    Real mTimeTotal, mScale;
    int mFrame;
        
    GridStorage<int> mGridsInt;
    GridStorage<Real> mGridsReal;
    GridStorage<Vec3> mGridsVec;

    // for timing plugins
    MuTime mPluginTimer;
    std::string mLastPlugin;
    std::vector<std::pair<std::string, MuTime> > mTimings;
    std::map<std::string, std::pair<int,MuTime> > mTimingsTotal;
};

}

#endif