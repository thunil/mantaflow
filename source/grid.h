/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Grid representation
 *
 ******************************************************************************/

#ifndef _GRID_H
#define _GRID_H

#include "manta.h"
#include "vectorbase.h"
#include "interpol.h"
#include "interpolHigh.h"
#include "kernel.h"

namespace Manta {
class LevelsetGrid;
	
//! Base class for all grids
PYTHON class GridBase : public PbClass {
public:
	enum GridType { TypeNone = 0, TypeReal = 1, TypeInt = 2, TypeVec3 = 4, TypeMAC = 8, TypeLevelset = 16, TypeFlags = 32 };
		
	PYTHON GridBase(FluidSolver* parent);
	
	//! Get the grids X dimension
	inline int getSizeX() const { return mSize.x; }
	//! Get the grids Y dimension
	inline int getSizeY() const { return mSize.y; }
	//! Get the grids Z dimension
	inline int getSizeZ() const { return mSize.z; }
	//! Get the grids dimensions
	inline Vec3i getSize() const { return mSize; }
	
	//! Get Stride in X dimension
	inline int getStrideX() const { return 1; }
	//! Get Stride in Y dimension
	inline int getStrideY() const { return mSize.x; }
	//! Get Stride in Z dimension
	inline int getStrideZ() const { return mStrideZ; }
	
	inline Real getDx() { return mDx; }
	
	//! Check if indices are within bounds, otherwise error (should only be called when debugging)
	inline void checkIndex(int i, int j, int k) const;
	//! Check if indices are within bounds, otherwise error (should only be called when debugging)
	inline void checkIndex(int idx) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec3i& p, int bnd) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec3i& p) const;
	//! Check if index is within given boundaries
	inline bool isInBounds(const Vec3& p, int bnd = 0) const { return isInBounds(toVec3i(p), bnd); }
	//! Check if linear index is in the range of the array
	inline bool isInBounds(int idx) const;
	
	//! Get the type of grid
	inline GridType getType() const { return mType; }
	//! Check dimensionality
	inline bool is2D() const { return !m3D; }
	//! Check dimensionality
	inline bool is3D() const { return m3D; }
	
	//! Get index into the data
	inline int index(int i, int j, int k) const { DEBUG_ONLY(checkIndex(i,j,k)); return i + mSize.x * j + mStrideZ * k; }
	//! Get index into the data
	inline int index(const Vec3i& pos) const    { DEBUG_ONLY(checkIndex(pos.x,pos.y,pos.z)); return pos.x + mSize.x * pos.y + mStrideZ * pos.z; }
protected:
	
	GridType mType;
	Vec3i mSize;
	Real mDx;
	bool m3D;
	// precomputed Z shift: to ensure 2D compatibility, always use this instead of sx*sy !
	int mStrideZ; 
};

//! Grid class
PYTHON template<class T>
class Grid : public GridBase {
public:
	//! init new grid, values are set to zero
	PYTHON Grid(FluidSolver* parent, bool show = true);
	//! create new & copy content from another grid
	Grid(const Grid<T>& a);
	//! return memory to solver
	virtual ~Grid();
	
	typedef T BASETYPE;
	
	PYTHON void save(std::string name);
	PYTHON void load(std::string name);
	
	//! set all cells to zero
	PYTHON void clear();
	
	//! all kinds of access functions, use grid(), grid[] or grid.get()
	//! access data
	inline T get(int i,int j, int k) const         { return mData[index(i,j,k)]; }
	//! access data
	inline T& get(int i,int j, int k)              { return mData[index(i,j,k)]; }
	//! access data
	inline T get(int idx) const                    { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T get(const Vec3i& pos) const           { return mData[index(pos)]; }
	//! access data
	inline T& operator()(int i, int j, int k)      { return mData[index(i, j, k)]; }
	//! access data
	inline T operator()(int i, int j, int k) const { return mData[index(i, j, k)]; }
	//! access data
	inline T& operator()(int idx)                  { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T operator()(int idx) const             { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline T& operator()(const Vec3i& pos)         { return mData[index(pos)]; }
	//! access data
	inline T operator()(const Vec3i& pos) const    { return mData[index(pos)]; }
	//! access data
	inline T& operator[](int idx)                  { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	//! access data
	inline const T operator[](int idx) const       { DEBUG_ONLY(checkIndex(idx)); return mData[idx]; }
	
	// interpolated access
	inline T    getInterpolated(const Vec3& pos) const { return interpol<T>(mData, mSize, mStrideZ, pos); }
	inline void setInterpolated(const Vec3& pos, const T& val, Grid<Real>& sumBuffer) const { setInterpol<T>(mData, mSize, mStrideZ, pos, val, &sumBuffer[0]); }
	// higher order interpolation (1=linear, 2=cubic)
	inline T getInterpolatedHi(const Vec3& pos, int order) const { 
		switch(order) {
		case 1:  return interpol     <T>(mData, mSize, mStrideZ, pos); 
		case 2:  return interpolCubic<T>(mData, mSize, mStrideZ, pos); 
		default: assertMsg(false, "Unknown interpolation order "<<order); }
	}
	
	// assignment / copy

	//! warning - do not use "=" for grids in python, this copies the reference! not the grid content...
	//Grid<T>& operator=(const Grid<T>& a);
	//! copy content from other grid (use this one instead of operator= !)
	PYTHON Grid<T>& copyFrom(const Grid<T>& a); // { *this = a; }

	// helper functions to work with grids in scene files 

	//! add/subtract other grid
	PYTHON void add(const Grid<T>& a);
	PYTHON void sub(const Grid<T>& a);
	//! set all cells to constant value
	PYTHON void setConst(T s);
	//! add constant to all grid cells
	PYTHON void addConst(T s);
	//! add scaled other grid to current one (note, only "Real" factor, "T" type not supported here!)
	PYTHON void addScaled(const Grid<T>& a, const T& factor); 
	//! multiply contents of grid
	PYTHON void mult( const Grid<T>& a);
	//! multiply each cell by a constant scalar value
	PYTHON void multConst(T s);
	//! clamp content to range (for vec3, clamps each component separately)
	PYTHON void clamp(Real min, Real max);
	
	// common compound operators
	//! get absolute max value in grid 
	PYTHON Real getMaxAbs();
	//! get max value in grid 
	PYTHON Real getMax();
	//! get min value in grid 
	PYTHON Real getMin();    
	//! set all boundary cells to constant value (Dirichlet)
	PYTHON void setBound(T value, int boundaryWidth=1);
	//! set all boundary cells to last inner value (Neumann)
	PYTHON void setBoundNeumann(int boundaryWidth=1);

	//! for compatibility, old names:
	PYTHON Real getMaxAbsValue() { return getMaxAbs(); }
	PYTHON Real getMaxValue()    { return getMax(); }
	PYTHON Real getMinValue()    { return getMin(); }

	//! debugging helper, print grid from python
	PYTHON void printGrid(int zSlice=-1,  bool printIndex=false); 

	// c++ only operators
	template<class S> Grid<T>& operator+=(const Grid<S>& a);
	template<class S> Grid<T>& operator+=(const S& a);
	template<class S> Grid<T>& operator-=(const Grid<S>& a);
	template<class S> Grid<T>& operator-=(const S& a);
	template<class S> Grid<T>& operator*=(const Grid<S>& a);
	template<class S> Grid<T>& operator*=(const S& a);
	template<class S> Grid<T>& operator/=(const Grid<S>& a);
	template<class S> Grid<T>& operator/=(const S& a);
	Grid<T>& safeDivide(const Grid<T>& a);    
	
	//! Swap data with another grid (no actual data is moved)
	void swap(Grid<T>& other);

protected:
	T* mData;
};

// Python doesn't know about templates: explicit aliases needed
PYTHON alias Grid<int>  IntGrid;
PYTHON alias Grid<Real> RealGrid;
PYTHON alias Grid<Vec3> VecGrid;

//! Special function for staggered grids
PYTHON class MACGrid : public Grid<Vec3> {
public:
	PYTHON MACGrid(FluidSolver* parent, bool show=true) : Grid<Vec3>(parent, show) { 
		mType = (GridType)(TypeMAC | TypeVec3); }
	
	// specialized functions for interpolating MAC information
	inline Vec3 getCentered(int i, int j, int k) const;
	inline Vec3 getCentered(const Vec3i& pos) const { return getCentered(pos.x, pos.y, pos.z); }
	inline Vec3 getAtMACX(int i, int j, int k) const;
	inline Vec3 getAtMACY(int i, int j, int k) const;
	inline Vec3 getAtMACZ(int i, int j, int k) const;
	// interpolation
	inline Vec3 getInterpolated(const Vec3& pos) const { return interpolMAC(mData, mSize, mStrideZ, pos); }
	inline void setInterpolated(const Vec3& pos, const Vec3& val, Vec3* tmp) { return setInterpolMAC(mData, mSize, mStrideZ, pos, val, tmp); }
	inline Vec3 getInterpolatedHi(const Vec3& pos, int order) const { 
		switch(order) {
		case 1:  return interpolMAC     (mData, mSize, mStrideZ, pos); 
		case 2:  return interpolCubicMAC(mData, mSize, mStrideZ, pos); 
		default: assertMsg(false, "Unknown interpolation order "<<order); }
	}
	// specials for mac grid:
	template<int comp> inline Real getInterpolatedComponent(Vec3 pos) const { return interpolComponent<comp>(mData, mSize, mStrideZ, pos); }
	template<int comp> inline Real getInterpolatedComponentHi(const Vec3& pos, int order) const { 
		switch(order) {
		case 1:  return interpolComponent<comp>(mData, mSize, mStrideZ, pos); 
		case 2:  return interpolCubicMAC(mData, mSize, mStrideZ, pos)[comp];  // warning - not yet optimized
		default: assertMsg(false, "Unknown interpolation order "<<order); }
	}
	
protected:
};

//! Special functions for FlagGrid
PYTHON class FlagGrid : public Grid<int> {
public:
	PYTHON FlagGrid(FluidSolver* parent, int dim=3, bool show=true) : Grid<int>(parent, show), mBoundaryWidth(0) { 
		mType = (GridType)(TypeFlags | TypeInt); }
	
	//! types of cells, in/outflow can be combined, e.g., TypeFluid|TypeInflow
	enum CellType { 
		TypeNone = 0,
		TypeFluid = 1,
		TypeObstacle = 2,
		TypeEmpty = 4,
		TypeInflow = 8,
		TypeOutflow = 16,
		TypeStick = 128,
		TypeReserved = 256,
		// 2^10 - 2^14 reserved for moving obstacles
		TypeZeroPressure = (1<<15) 
	};
		
	//! access for particles
	inline int getAt(const Vec3& pos) const { return mData[index((int)pos.x, (int)pos.y, (int)pos.z)]; }
			
	//! check for different flag types
	inline bool isObstacle(int idx) const { return get(idx) & TypeObstacle; }
	inline bool isObstacle(int i, int j, int k) const { return get(i,j,k) & TypeObstacle; }
	inline bool isObstacle(const Vec3i& pos) const { return get(pos) & TypeObstacle; }
	inline bool isObstacle(const Vec3& pos) const { return getAt(pos) & TypeObstacle; }
	inline bool isFluid(int idx) const { return get(idx) & TypeFluid; }
	inline bool isFluid(int i, int j, int k) const { return get(i,j,k) & TypeFluid; }
	inline bool isFluid(const Vec3i& pos) const { return get(pos) & TypeFluid; }
	inline bool isFluid(const Vec3& pos) const { return getAt(pos) & TypeFluid; }
	inline bool isInflow(int idx) const { return get(idx) & TypeInflow; }
	inline bool isInflow(int i, int j, int k) const { return get(i,j,k) & TypeInflow; }
	inline bool isInflow(const Vec3i& pos) const { return get(pos) & TypeInflow; }
	inline bool isInflow(const Vec3& pos) const { return getAt(pos) & TypeInflow; }
	inline bool isEmpty(int idx) const { return get(idx) & TypeEmpty; }
	inline bool isEmpty(int i, int j, int k) const { return get(i,j,k) & TypeEmpty; }
	inline bool isEmpty(const Vec3i& pos) const { return get(pos) & TypeEmpty; }
	inline bool isEmpty(const Vec3& pos) const { return getAt(pos) & TypeEmpty; }
	inline bool isStick(int idx) const { return get(idx) & TypeStick; }
	inline bool isStick(int i, int j, int k) const { return get(i,j,k) & TypeStick; }
	inline bool isStick(const Vec3i& pos) const { return get(pos) & TypeStick; }
	inline bool isStick(const Vec3& pos) const { return getAt(pos) & TypeStick; }
	
	inline int getBoundaryWidth() const {return mBoundaryWidth;}

	// Python callables
	PYTHON void initDomain(int boundaryWidth=0);
	PYTHON void initBoundaries(int boundaryWidth=0);
	PYTHON void updateFromLevelset(LevelsetGrid& levelset);    
	PYTHON void fillGrid(int type=TypeFluid);

protected:
	int mBoundaryWidth;

};

//! helper to compute grid conversion factor between local coordinates of two grids
inline Vec3 calcGridSizeFactor(Vec3i s1, Vec3i s2) {
	return Vec3( Real(s1[0])/s2[0], Real(s1[1])/s2[1], Real(s1[2])/s2[2] );
}


//******************************************************************************
// enable compilation of a more complicated test data type
// for grids... note - this also enables code parts in fileio.cpp!
// the code below is meant only as an example for a grid with a more complex data type
// and illustrates which functions need to be implemented; it's not needed
// to run any simulations in mantaflow!

#define ENABLE_GRID_TEST_DATATYPE 0

#if ENABLE_GRID_TEST_DATATYPE==1

typedef std::vector<int> nbVectorBaseType;
class nbVector : public nbVectorBaseType {
	public: 
		inline nbVector() : nbVectorBaseType() {};
		inline ~nbVector() {};

		// grid operators require certain functions
		inline nbVector(Real v) : nbVectorBaseType() { this->push_back( (int)v ); };

		inline const nbVector& operator+= ( const nbVector &v1 ) {
			assertMsg(false,"Never call!"); return *this; 
		}
		inline const nbVector& operator-= ( const nbVector &v1 ) {
			assertMsg(false,"Never call!"); return *this; 
		}
		inline const nbVector& operator*= ( const nbVector &v1 ) {
			assertMsg(false,"Never call!"); return *this; 
		}
};

template<> inline nbVector* FluidSolver::getGridPointer<nbVector>() {
	return new nbVector[mGridSize.x * mGridSize.y * mGridSize.z];
}
template<> inline void FluidSolver::freeGridPointer<nbVector>(nbVector* ptr) {
	return delete[] ptr;
}

inline nbVector operator+ ( const nbVector &v1, const nbVector &v2 ) {
	assertMsg(false,"Never call!"); return nbVector(); 
}
inline nbVector operator* ( const nbVector &v1, const nbVector &v2 ) {
	assertMsg(false,"Never call!"); return nbVector(); 
}
template<class S>
inline nbVector operator* ( const nbVector& v, S s ) {
	assertMsg(false,"Never call!"); return nbVector(); 
} 
template<class S> 
inline nbVector operator* ( S s, const nbVector& v ) {
	assertMsg(false,"Never call!"); return nbVector(); 
}

template<> inline nbVector safeDivide<nbVector>(const nbVector &a, const nbVector& b) { 
	assertMsg(false,"Never call!"); return nbVector(); 
}

std::ostream& operator<< ( std::ostream& os, const nbVectorBaseType& i ) {
	os << " nbVectorBaseType NYI ";
	return os;
}

// make data type known to python
// (python keyword changed here, because the preprocessor does not yet parse #ifdefs correctly)
PY THON alias Grid<nbVector> TestDataGrid;
// ? PY THON alias nbVector TestDatatype;

#endif // end ENABLE_GRID_TEST_DATATYPE



//******************************************************************************
// Implementation of inline functions

inline void GridBase::checkIndex(int i, int j, int k) const {
	//if (i<0 || j<0  || i>=mSize.x || j>=mSize.y || (is3D() && (k<0|| k>= mSize.z))) {
	if (i<0 || j<0  || i>=mSize.x || j>=mSize.y || k<0|| k>= mSize.z ) {
		std::ostringstream s;
		s << "Grid " << mName << " dim " << mSize << " : index " << i << "," << j << "," << k << " out of bound ";
		errMsg(s.str());
	}
}

inline void GridBase::checkIndex(int idx) const {
	if (idx<0 || idx >= mSize.x * mSize.y * mSize.z) {
		std::ostringstream s;
		s << "Grid " << mName << " dim " << mSize << " : index " << idx << " out of bound ";
		errMsg(s.str());
	}
}

bool GridBase::isInBounds(const Vec3i& p) const { 
	return (p.x >= 0 && p.y >= 0 && p.z >= 0 && p.x < mSize.x && p.y < mSize.y && p.z < mSize.z); 
}

bool GridBase::isInBounds(const Vec3i& p, int bnd) const { 
	bool ret = (p.x >= bnd && p.y >= bnd && p.x < mSize.x-bnd && p.y < mSize.y-bnd);
	if(this->is3D()) {
		ret &= (p.z >= bnd && p.z < mSize.z-bnd); 
	} else {
		ret &= (p.z == 0);
	}
	return ret;
}
//! Check if linear index is in the range of the array
bool GridBase::isInBounds(int idx) const {
	if (idx<0 || idx >= mSize.x * mSize.y * mSize.z) {
		return false;
	}
	return true;
}

inline Vec3 MACGrid::getCentered(int i, int j, int k) const {
	DEBUG_ONLY(checkIndex(i+1,j+1,k));
	const int idx = index(i,j,k);
	Vec3 v = Vec3(0.5* (mData[idx].x + mData[idx+1].x),
				  0.5* (mData[idx].y + mData[idx+mSize.x].y),
				  0.);
	if( this->is3D() ) {
		DEBUG_ONLY(checkIndex(idx+mStrideZ));
		v[2] =    0.5* (mData[idx].z + mData[idx+mStrideZ].z);
	}
	return v;
}

inline Vec3 MACGrid::getAtMACX(int i, int j, int k) const {
	DEBUG_ONLY(checkIndex(i-1,j+1,k));
	const int idx = index(i,j,k);
	Vec3 v =  Vec3(   (mData[idx].x),
				0.25* (mData[idx].y + mData[idx-1].y + mData[idx+mSize.x].y + mData[idx+mSize.x-1].y),
				0.);
	if( this->is3D() ) {
		DEBUG_ONLY(checkIndex(idx+mStrideZ-1));
		v[2] = 0.25* (mData[idx].z + mData[idx-1].z + mData[idx+mStrideZ].z + mData[idx+mStrideZ-1].z);
	}
	return v;
}

inline Vec3 MACGrid::getAtMACY(int i, int j, int k) const {
	DEBUG_ONLY(checkIndex(i+1,j-1,k));
	const int idx = index(i,j,k);
	Vec3 v =  Vec3(0.25* (mData[idx].x + mData[idx-mSize.x].x + mData[idx+1].x + mData[idx+1-mSize.x].x),
						 (mData[idx].y),   0. );
	if( this->is3D() ) {
		DEBUG_ONLY(checkIndex(idx+mStrideZ-mSize.x));
		v[2] = 0.25* (mData[idx].z + mData[idx-mSize.x].z + mData[idx+mStrideZ].z + mData[idx+mStrideZ-mSize.x].z);
	}
	return v;
}

inline Vec3 MACGrid::getAtMACZ(int i, int j, int k) const {
	const int idx = index(i,j,k);
	DEBUG_ONLY(checkIndex(idx-mStrideZ));
	DEBUG_ONLY(checkIndex(idx+mSize.x-mStrideZ));
	Vec3 v =  Vec3(0.25* (mData[idx].x + mData[idx-mStrideZ].x + mData[idx+1].x + mData[idx+1-mStrideZ].x),
				   0.25* (mData[idx].y + mData[idx-mStrideZ].y + mData[idx+mSize.x].y + mData[idx+mSize.x-mStrideZ].y),
						 (mData[idx].z) );
	return v;
}

KERNEL(idx) template<class T, class S> void gridAdd  (Grid<T>& me, const Grid<S>& other) { me[idx] += other[idx]; }
KERNEL(idx) template<class T, class S> void gridSub  (Grid<T>& me, const Grid<S>& other) { me[idx] -= other[idx]; }
KERNEL(idx) template<class T, class S> void gridMult (Grid<T>& me, const Grid<S>& other) { me[idx] *= other[idx]; }
KERNEL(idx) template<class T, class S> void gridDiv  (Grid<T>& me, const Grid<S>& other) { me[idx] /= other[idx]; }
KERNEL(idx) template<class T, class S> void gridAddScalar (Grid<T>& me, const S& other)  { me[idx] += other; }
KERNEL(idx) template<class T, class S> void gridMultScalar(Grid<T>& me, const S& other)  { me[idx] *= other; }
KERNEL(idx) template<class T, class S> void gridScaledAdd (Grid<T>& me, const Grid<T>& other, const S& factor) { me[idx] += factor * other[idx]; }

KERNEL(idx) template<class T> void gridSafeDiv (Grid<T>& me, const Grid<T>& other) { me[idx] = safeDivide(me[idx], other[idx]); }
KERNEL(idx) template<class T> void gridSetConst(Grid<T>& grid, T value) { grid[idx] = value; }

template<class T> template<class S> Grid<T>& Grid<T>::operator+= (const Grid<S>& a) {
	gridAdd<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator+= (const S& a) {
	gridAddScalar<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator-= (const Grid<S>& a) {
	gridSub<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator-= (const S& a) {
	gridAddScalar<T,S> (*this, -a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator*= (const Grid<S>& a) {
	gridMult<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator*= (const S& a) {
	gridMultScalar<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator/= (const Grid<S>& a) {
	gridDiv<T,S> (*this, a);
	return *this;
}
template<class T> template<class S> Grid<T>& Grid<T>::operator/= (const S& a) {
	S rez((S)1.0 / a);
	gridMultScalar<T,S> (*this, rez);
	return *this;
}


//******************************************************************************
// Other helper functions

// compute gradient of a scalar grid
inline Vec3 getGradient(const Grid<Real>& data, int i, int j, int k) {
	Vec3 v;

	if (i > data.getSizeX()-2) i= data.getSizeX()-2;
	if (j > data.getSizeY()-2) j= data.getSizeY()-2;
	if (i < 1) i = 1;
	if (j < 1) j = 1;
	v = Vec3( data(i+1,j  ,k  ) - data(i-1,j  ,k  ) ,
			  data(i  ,j+1,k  ) - data(i  ,j-1,k  ) , 0. );

	if(data.is3D()) {
		if (k > data.getSizeZ()-2) k= data.getSizeZ()-2;
		if (k < 1) k = 1;
		v[2]= data(i  ,j  ,k+1) - data(i  ,j  ,k-1);
	} 

	return v;
}

// interpolate grid from one size to another size
KERNEL template<class S>
void knInterpolateGridTempl(Grid<S>& target, Grid<S>& source, const Vec3& sourceFactor , Vec3 offset, int orderSpace=1 ) {
	Vec3 pos = Vec3(i,j,k) * sourceFactor + offset;
	if(!source.is3D()) pos[2] = 0; // allow 2d -> 3d
	target(i,j,k) = source.getInterpolated(pos);
} 
// template glue code - choose interpolation based on template arguments
template<class GRID>
void interpolGridTempl( GRID& target, GRID& source ) {
		errMsg("interpolGridTempl - Only valid for specific instantiations");
}


} //namespace
#endif
