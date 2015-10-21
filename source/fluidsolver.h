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

#include "manta.h"
#include "vectorbase.h"
#include <vector>
#include <map>

namespace Manta { 
	
//! Encodes grid size, timstep etc.
PYTHON(name=Solver) 
class FluidSolver : public PbClass {
public:
	PYTHON() FluidSolver(Vec3i gridSize, int dim=3);
	virtual ~FluidSolver();
	
	// accessors
	PYTHON() Vec3i getGridSize() { return mGridSize; }
	inline Real  getDt()       { return mDt; }
	inline Real  getDx()       { return 1.0 / mGridSize.max(); }
	inline Real  getTime()     { return mTimeTotal; }

	//! Check dimensionality
	inline bool is2D() const { return mDim==2; }
	//! Check dimensionality
	inline bool is3D() const { return mDim==3; }
	
	PYTHON() void printMemInfo();
	
	//! Advance the solver one timestep, update GUI if present
	PYTHON() void step();
	
	//! Update the timestep size based on given maximal velocity magnitude 
	PYTHON() void adaptTimestep(Real maxVel);
	
	//! create a object with the solver as its parent
	PYTHON() PbClass* create(PbType type, PbTypeVec T=PbTypeVec(),const std::string& name = "");
	
	// temp grid and plugin stuff: you shouldn't call this manually
	template<class T> T* getGridPointer();
	template<class T> void freeGridPointer(T* ptr);    

	//! expose animation time to python
	PYTHON(name=timestep)  Real mDt;  
	PYTHON(name=timeTotal) Real mTimeTotal;
	PYTHON(name=frame)     int  mFrame;
	//! parameters for adaptive time stepping
	PYTHON(name=cfl)          Real mCflCond;  
	PYTHON(name=timestepMin)  Real mDtMin;  
	PYTHON(name=timestepMax)  Real mDtMax;  
	PYTHON(name=frameLength)  Real mFrameLength;

protected:
	Vec3i     mGridSize;
	const int mDim;
	Real      mTimePerFrame;
	bool      mLockDt;
	bool      mAdaptDt;
		
	//! subclass for managing grid memory
	//! stored as a stack to allow fast allocation
	template<class T> struct GridStorage {
		GridStorage() : used(0) {}
		T* get(Vec3i size);
		void free();
		void release(T* ptr);
		
		std::vector<T*> grids;
		int used;
	};
	
	GridStorage<int> mGridsInt;
	GridStorage<Real> mGridsReal;
	GridStorage<Vec3> mGridsVec;
};

}

#endif
