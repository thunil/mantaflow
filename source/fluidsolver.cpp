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

#include "fluidsolver.h"
#include "grid.h"
#include <sstream>
#include <fstream>

using namespace std;
namespace Manta {

#ifdef GUI
    // defined in qtmain.cpp
    extern void updateQtGui(bool full, int frame, const std::string& curPlugin);
#else
    inline void updateQtGui(bool full, int frame, const std::string& curPlugin) {}
#endif

//******************************************************************************
// Gridstorage-related members

template<class T>
void FluidSolver::GridStorage<T>::free() {
    if (used != 0)
        errMsg("can't clean grid cache, some grids are still in use");
    for(size_t i = 0; i<grids.size(); i++)
        delete[] grids[i];
    grids.clear();
}
template<class T>
T* FluidSolver::GridStorage<T>::get(Vec3i size) {
    if ((int)grids.size() <= used) {
        grids.push_back(new T[size.x * size.y * size.z]);
    }
    if (used > 200)
        errMsg("too many temp grids used -- are they released properly ?");
    return grids[used++];
}
template<class T>
void FluidSolver::GridStorage<T>::release(T* ptr) {
    // rewrite pointer, as it may have changed due to swap operations
    used--;
    if (used < 0)
        errMsg("temp grid inconsistency");
    grids[used] = ptr;
}

template<> int* FluidSolver::getGridPointer<int>() {
    return mGridsInt.get(mGridSize);    
}
template<> Real* FluidSolver::getGridPointer<Real>() {
    return mGridsReal.get(mGridSize);    
}
template<> Vec3* FluidSolver::getGridPointer<Vec3>() {
    return mGridsVec.get(mGridSize);    
}
template<> void FluidSolver::freeGridPointer<int>(int *ptr) {
    mGridsInt.release(ptr);
}
template<> void FluidSolver::freeGridPointer<Real>(Real* ptr) {
    mGridsReal.release(ptr);
}
template<> void FluidSolver::freeGridPointer<Vec3>(Vec3* ptr) {
    mGridsVec.release(ptr);
}

void FluidSolver::pluginStart(const string& name) {
    mLastPlugin = name;
    mPluginTimer.get();
}

void FluidSolver::pluginStop(const string& name) {
    if (mLastPlugin == name && name != "FluidSolver::step") {
        MuTime diff = mPluginTimer.update();
        mTimings.push_back(pair<string,MuTime>(name, diff));
    }    
}

//******************************************************************************
// FluidSolver members

FluidSolver::FluidSolver(Vec3i gridsize, int dim)
    : PbClass(this), mDt(1.0), mScale(1.0), mFrame(0), mGridSize(gridsize), mDim(dim), mTimeTotal(0.)
{    
    assertMsg(dim==2 || dim==3, "Can only create 2D and 3D solvers");
    assertMsg(dim!=2 || gridsize.z == 1, "Trying to create 2D solver with size.z != 1");
}

FluidSolver::~FluidSolver() {
    mGridsInt.free();
    mGridsReal.free();
    mGridsVec.free();
}

PbClass* FluidSolver::create(PbType t, const string& name) {        
    _args.add("nocheck",true);
    if (t.str == "")
        errMsg("Need to specify object type. Use e.g. Solver.create(FlagGrid, ...) or Solver.create(type=FlagGrid, ...)");
    
    return PbClass::createPyObject(t.str, name, _args, this);
}

void FluidSolver::step() {
    mTimeTotal += mDt;
    mFrame++;
    updateQtGui(true, mFrame, "FluidSolver::step");
    
    // update timings
    for(size_t i=0;i<mTimings.size(); i++) {
        const string name = mTimings[i].first;
        if (mTimingsTotal.find(name) == mTimingsTotal.end())
            mTimingsTotal[name].second.clear();
        mTimingsTotal[name].first++;
        mTimingsTotal[name].second+=mTimings[i].second;
    }
    mTimings.clear();
}
 
void FluidSolver::printTimings() {
    MuTime total;
    total.clear();
    for(size_t i=0; i<mTimings.size(); i++)
        total += mTimings[i].second;
    
    printf("\n-- STEP %d -----------------------------\n", mFrame);
    for(size_t i=0; i<mTimings.size(); i++)
        printf("[%4.1f%%] %s (%s)\n", 100.0*((Real)mTimings[i].second.time / (Real)total.time), mTimings[i].first.c_str(), mTimings[i].second.toString().c_str());
    printf("----------------------------------------\n");
    printf("Total : %s\n\n", total.toString().c_str());
}

void FluidSolver::saveMeanTimings(string filename) {
    ofstream ofs(filename.c_str());
    if (!ofs.good())
        errMsg("can't open " + filename + " as timing log");
    ofs << "Mean timings of " << mFrame << " steps :" <<endl <<endl;
    MuTime total;
    total.clear();
    for(map<string, pair<int,MuTime> >::iterator it=mTimingsTotal.begin(); it!=mTimingsTotal.end(); it++) {
        total += it->second.second;
    }    
    for(map<string, pair<int,MuTime> >::iterator it=mTimingsTotal.begin(); it!=mTimingsTotal.end(); it++) {
        ofs << it->first << ": " << it->second.second / it->second.first << endl;        
    }
    ofs << endl << "Total : " << total << " (mean " << total/mFrame << ")" << endl;
    ofs.close();
}
 
}

