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

#include "pclass.h"
#include "vectorbase.h"
#include "interpol.h"
#include "kernel.h"

namespace Manta {
template<int DIM> class LevelsetGrid;
    
//! Base class for all grids
PYTHON class GridBase : public PbClass {
public:
    enum GridType { TypeNone = 0, TypeReal = 1, TypeInt = 2, TypeVec3 = 4, TypeMAC = 8, TypeLevelset = 16, TypeFlags = 32, Type2D = 64};
        
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
    
    //! Check if indices are within bounds, otherwise error
    void checkIndex(int i, int j, int k) const;
    //! Check if indices are within bounds, otherwise error
    void checkIndex(int idx) const;
    //! Check if index is within given boundaries
    bool isInBounds(const Vec3i& pos, int bnd);
    //! Check if index is within given boundaries
    bool isInBounds(const Vec3i& pos);
    //! Check if index is within given boundaries
    bool isInBounds(const Vec3& pos, int bnd = 0);
    
    //! Get the type of grid
    inline GridType getType() { return mType; }
    
    //! Get index into the data
    inline int index(int i, int j, int k) const;
    //! Get index into the data
    inline int index(const Vec3i& pos) const;
protected:
    
    GridType mType;
    Vec3i mSize;
    Real mDx;
    int mStrideZ; // precomputed for speed
};

//! Grid class
PYTHON template<int DIM, class T>
class Grid : public GridBase {
public:
    PYTHON Grid(FluidSolver* parent, bool show = true);
    virtual ~Grid();
    
    typedef T BASETYPE;
    
    PYTHON void save(std::string name);
    
    //! set all cells to zero
    void clear();
    
    //! access data
    inline T& get(int i,int j, int k) { return mData[index(i,j,k)]; }
    //! access data
    inline T get(int idx) const;
    //! access data
    inline T get(const Vec3i& pos) const { return mData[index(pos)]; }
    //! access data
    inline T& operator()(int i, int j, int k) { return mData[index(i, j, k)]; }
    //! access data
    inline T operator()(int i, int j, int k) const { return mData[index(i, j, k)]; }
    //! access data
    inline T& operator()(const Vec3i& pos) { return mData[index(pos)]; }
    //! access data
    inline T operator()(const Vec3i& pos) const { return mData[index(pos)]; }
    //! access data
    inline T& operator[](int idx);
    //! access data
    inline const T operator[](int idx) const;
    
    // interpolated access
    inline T getInterpolated(const Vec3& pos) const { return interpol<DIM,T>(mData, mSize, pos); }
    inline void setInterpolated(const Vec3& pos, const T& val, Grid<DIM,Real>& sumBuffer) const { setInterpol<DIM,T>(mData, mSize, pos, val, &sumBuffer[0]); }
    
    // operators
    template<class S> Grid<DIM,T>& operator+=(const Grid<DIM,S>& a);
    Grid<DIM,T>& operator+=(const T& a);
    template<class S> Grid<DIM,T>& operator-=(const Grid<DIM,S>& a);
    template<class S> Grid<DIM,T>& operator-=(const S& a);
    template<class S> Grid<DIM,T>& operator*=(const Grid<DIM,S>& a);
    Grid<DIM,T>& operator*=(const T& a);
    template<class S> Grid<DIM,T>& operator/=(const Grid<DIM,S>& a);
    template<class S> Grid<DIM,T>& operator/=(const S& a);
    Grid<DIM,T>& operator=(const T& a);
    Grid<DIM,T>& operator=(const Grid<DIM,T>& a);
    Grid<DIM,T>& safeDivide(const Grid<DIM,T>& a);    
    PYTHON void add(const Grid<DIM,T>& a, const Grid<DIM,T>& b);
    PYTHON void scale(T a) { (*this) *= a; }
    PYTHON void copyFrom(const Grid<DIM,T>& a) { *this = a; }
    
    // common compound operators
    //! Grid += a*factor
    void scaledAdd(const Grid<DIM,T>& a, const T& factor);
    //! get absolute max value in grid (only Real grids)
    Real getMaxAbsValue();
    //! get max value in grid (only Real grids)
    Real getMaxValue();
    //! get min value in grid (only Real grids)
    Real getMinValue();    
    //! Swap data with another grid (no actual data is moved)
    void swap(Grid<DIM,T>& other);
    
protected:
    T* mData;
};

// Python doesnt have templates: need to create aliases
PYTHON alias Grid<2,int> IntGrid2;
PYTHON alias Grid<2,Real> RealGrid2;
PYTHON alias Grid<2,Vec3> VecGrid2;
PYTHON alias Grid<3,int> IntGrid3;
PYTHON alias Grid<3,Real> RealGrid3;
PYTHON alias Grid<3,Vec3> VecGrid3;

//! Special function for staggered grids
PYTHON template<int DIM>
class MACGrid : public Grid<DIM, Vec3> {
public:
    PYTHON MACGrid(FluidSolver* parent, bool show=true) : Grid<DIM,Vec3>(parent, show) { this->mType = (GridBase::GridType)(this->TypeMAC | this->TypeVec3 | ((DIM==2) ? this->Type2D : 0)); }
    
    // specialized functions for interpolating MAC information
    inline Vec3 getCentered(int i, int j, int k);
    inline Vec3 getCentered(const Vec3i& pos) { return getCentered(pos.x, pos.y, pos.z); }
    inline Vec3 getAtMACX(int i, int j, int k);
    inline Vec3 getAtMACY(int i, int j, int k);
    inline Vec3 getAtMACZ(int i, int j, int k);
    template<int comp> inline Real getInterpolatedComponent(Vec3 pos) const { return interpolComponent<DIM,comp>(Grid<DIM,Vec3>::mData, GridBase::mSize, pos); }
    inline Vec3 getInterpolated(const Vec3& pos) const { return interpolMAC<DIM>(Grid<DIM,Vec3>::mData, GridBase::mSize, pos); }
    inline void setInterpolated(const Vec3& pos, const Vec3& val, Vec3* tmp) { return setInterpolMAC<DIM>(Grid<DIM,Vec3>::mData, GridBase::mSize, pos, val, tmp); }
    
protected:
};

//! Special functions for FlagGrid
PYTHON template<int DIM>
class FlagGrid : public Grid<DIM, int> {
public:
    PYTHON FlagGrid(FluidSolver* parent, bool show=true) : Grid<DIM,int>(parent, show) {  this->mType = (GridBase::GridType)(this->TypeFlags | this->TypeInt | ((DIM==2) ? this->Type2D : 0)); }
    
	//! types of cells, in/outflow can be combined, e.g., TypeFluid|TypeInflow
    enum CellType { 
        TypeNone = 0,
        TypeFluid = 1,
        TypeObstacle = 2,
        TypeEmpty = 4,
        TypeInflow = 8,
        TypeOutflow = 16,
		TypeStick = 128
	};
        
    //! access for particles
    inline int getAt(const Vec3& pos) const { return Grid<DIM,int>::mData[Grid<DIM,int>::index((int)pos.x, (int)pos.y, (int)pos.z)]; }
            
	//! check for different flag types
	using Grid<DIM,int>::get;
    inline bool isObstacle(int idx) { return get(idx) & TypeObstacle; }
    inline bool isObstacle(int i, int j, int k) { return get(i,j,k) & TypeObstacle; }
    inline bool isObstacle(const Vec3i& pos) { return get(pos) & TypeObstacle; }
    inline bool isObstacle(const Vec3& pos) { return getAt(pos) & TypeObstacle; }
    inline bool isFluid(int idx) { return get(idx) & TypeFluid; }
    inline bool isFluid(int i, int j, int k) { return get(i,j,k) & TypeFluid; }
    inline bool isFluid(const Vec3i& pos) { return get(pos) & TypeFluid; }
    inline bool isFluid(const Vec3& pos) { return getAt(pos) & TypeFluid; }
    inline bool isEmpty(int idx) { return get(idx) & TypeEmpty; }
    inline bool isEmpty(int i, int j, int k) { return get(i,j,k) & TypeEmpty; }
    inline bool isEmpty(const Vec3i& pos) { return get(pos) & TypeEmpty; }
    inline bool isEmpty(const Vec3& pos) { return getAt(pos) & TypeEmpty; }
    inline bool isStick(int idx) { return get(idx) & TypeStick; }
    inline bool isStick(int i, int j, int k) { return get(i,j,k) & TypeStick; }
    inline bool isStick(const Vec3i& pos) { return get(pos) & TypeStick; }
    inline bool isStick(const Vec3& pos) { return getAt(pos) & TypeStick; }
    
    // Python callables
    PYTHON void initDomain(int boundaryWidth=1);
    PYTHON void initBoundaries(int boundaryWidth=1);
    PYTHON void updateFromLevelset(LevelsetGrid<DIM>& levelset);    
    PYTHON void fillGrid();
};

// Python doesnt have templates: need to create aliases
PYTHON alias MACGrid<2> MACGrid2;
PYTHON alias MACGrid<3> MACGrid3;
PYTHON alias FlagGrid<2> FlagGrid2;
PYTHON alias FlagGrid<3> FlagGrid3;

// g++/MSVC don't support this C++11 thing yet... therefore ugly inheritance for now
// -> template <typename T> using Grid3 = Grid<3,T>;
PYTHON template<class T>
class Grid3 : public Grid<3,T> {
public:
    PYTHON Grid3(FluidSolver* parent, bool show = true) : Grid<3,T>(parent,show) {}
};
PYTHON template<class T>
class Grid2 : public Grid<2,T> {
public:
    PYTHON Grid2(FluidSolver* parent, bool show = true) : Grid<2,T>(parent,show) {}
};
PYTHON alias Grid2<int> IntGrid2;
PYTHON alias Grid2<Real> RealGrid2;
PYTHON alias Grid2<Vec3> VecGrid2;
PYTHON alias Grid3<int> IntGrid3;
PYTHON alias Grid3<Real> RealGrid3;
PYTHON alias Grid3<Vec3> VecGrid3;


//******************************************************************************
// Implementation of inline functions

inline int GridBase::index(int i, int j, int k) const { 
#ifdef DEBUG
    checkIndex(i, j, k);
#endif
    return i + mSize.x * (j + mSize.y * k);
}

inline int GridBase::index(const Vec3i& pos) const { 
#ifdef DEBUG
    checkIndex(pos.x, pos.y, pos.z);
#endif
    return pos.x + mSize.x * (pos.y + mSize.y * pos.z);
}

inline bool GridBase::isInBounds(const Vec3i& p, int bnd) {
    return (p.x >= bnd && p.y >= bnd && p.z >= bnd &&
            p.x < mSize.x-bnd && p.y < mSize.y-bnd && p.z < mSize.z-bnd);
}

inline bool GridBase::isInBounds(const Vec3i& p) {
    return (p.x >= 0 && p.y >= 0 && p.z >= 0 &&
            p.x < mSize.x && p.y < mSize.y && p.z < mSize.z);
}

inline bool GridBase::isInBounds(const Vec3& p, int bnd) {
    return isInBounds(toVec3i(p), bnd);
}

template<int DIM,class T> inline T Grid<DIM,T>::get(int idx) const { 
#ifdef DEBUG
    checkIndex(idx); 
#endif
    return mData[idx];         
}   

template<int DIM,class T> inline T& Grid<DIM,T>::operator[](int idx) { 
#ifdef DEBUG
    checkIndex(idx); 
#endif
    return mData[idx];         
}

template<int DIM, class T> inline const T Grid<DIM,T>::operator[](int idx) const { 
#ifdef DEBUG
    checkIndex(idx); 
#endif
    return mData[idx];         
}

template<> inline Vec3 MACGrid<3>::getCentered(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i+1,j+1,k+1);
#endif
    const int idx = index(i,j,k);
    return Vec3(0.5* (mData[idx].x + mData[idx+1].x),
                0.5* (mData[idx].y + mData[idx+mSize.x].y),
                0.5* (mData[idx].z + mData[idx+mStrideZ].z) );
}

template<> inline Vec3 MACGrid<2>::getCentered(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i+1,j+1,k);
#endif
    const int idx = index(i,j,k);
    return Vec3(0.5* (mData[idx].x + mData[idx+1].x),
                0.5* (mData[idx].y + mData[idx+mSize.x].y), 0);
}

template<> inline Vec3 MACGrid<3>::getAtMACX(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i-1,j+1,k+1);
#endif
    const int idx = index(i,j,k);
    return Vec3(      (mData[idx].x),
                0.25* (mData[idx].y + mData[idx-1].y + mData[idx+mSize.x].y + mData[idx+mSize.x-1].y),
                0.25* (mData[idx].z + mData[idx-1].z + mData[idx+mStrideZ].z + mData[idx+mStrideZ-1].z) );
}

template<> inline Vec3 MACGrid<2>::getAtMACX(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i-1,j+1,k+1);
#endif
    const int idx = index(i,j,k);
    return Vec3(      (mData[idx].x),
                0.25* (mData[idx].y + mData[idx-1].y + mData[idx+mSize.x].y + mData[idx+mSize.x-1].y), 0);
}

template<> inline Vec3 MACGrid<3>::getAtMACY(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i+1,j-1,k+1);
#endif
    const int idx = index(i,j,k);
    return Vec3(0.25* (mData[idx].x + mData[idx-mSize.x].x + mData[idx+1].x + mData[idx+1-mSize.x].x),
                      (mData[idx].y),
                0.25* (mData[idx].z + mData[idx-mSize.x].z + mData[idx+mStrideZ].z + mData[idx+mStrideZ-mSize.x].z) );
}

template<> inline Vec3 MACGrid<2>::getAtMACY(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i+1,j-1,k);
#endif
    const int idx = index(i,j,k);
    return Vec3(0.25* (mData[idx].x + mData[idx-mSize.x].x + mData[idx+1].x + mData[idx+1-mSize.x].x),
                      (mData[idx].y), 0);
}

template<> inline Vec3 MACGrid<3>::getAtMACZ(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i+1,j+1,k-1);
#endif
    const int idx = index(i,j,k);
    return Vec3(0.25* (mData[idx].x + mData[idx-mStrideZ].x + mData[idx+1].x + mData[idx+1-mStrideZ].x),
                0.25* (mData[idx].y + mData[idx-mStrideZ].y + mData[idx+mSize.x].y + mData[idx+mSize.x-mStrideZ].y),
                      (mData[idx].z) );
}

template<> inline Vec3 MACGrid<2>::getAtMACZ(int i, int j, int k) {
#ifdef DEBUG
    checkIndex(i+1,j+1,k);
#endif
    return Vec3(0,0,0);
}


KERNEL(idx) template<int DIM, class T> gridAdd2 (Grid<DIM,T>& me, const Grid<DIM,T>& a, const Grid<DIM,T>& b) { me[idx] = a[idx] + b[idx]; }
KERNEL(idx) template<int DIM, class T, class S> gridAdd (Grid<DIM,T>& me, const Grid<DIM,S>& other) { me[idx] += other[idx]; }
KERNEL(idx) template<int DIM, class T, class S> gridSub (Grid<DIM,T>& me, const Grid<DIM,S>& other) { me[idx] -= other[idx]; }
KERNEL(idx) template<int DIM, class T, class S> gridMult (Grid<DIM,T>& me, const Grid<DIM,S>& other) { me[idx] *= other[idx]; }
KERNEL(idx) template<int DIM, class T, class S> gridDiv (Grid<DIM,T>& me, const Grid<DIM,S>& other) { me[idx] /= other[idx]; }
KERNEL(idx) template<int DIM, class T> gridSafeDiv (Grid<DIM,T>& me, const Grid<DIM,T>& other) { me[idx] = safeDivide(me[idx], other[idx]); }
KERNEL(idx) template<int DIM, class T, class S> gridAddScalar (Grid<DIM,T>& me, const S& other) { me[idx] += other; }
KERNEL(idx) template<int DIM, class T, class S> gridMultScalar (Grid<DIM,T>& me, const S& other) { me[idx] *= other; }
KERNEL(idx) template<int DIM, class T> gridScaleAdd (Grid<DIM,T>& me, const Grid<DIM,T>& other, const T& factor) { me[idx] += factor * other[idx]; }

template<int DIM, class T> template<class S> Grid<DIM,T>& Grid<DIM,T>::operator+= (const Grid<DIM,S>& a) {
    gridAdd<DIM,T,S> (*this, a);
    return *this;
}
template<int DIM, class T> Grid<DIM,T>& Grid<DIM,T>::operator+= (const T& a) {
    gridAddScalar<DIM,T,T> (*this, a);
    return *this;
}
template<int DIM, class T> template<class S> Grid<DIM,T>& Grid<DIM,T>::operator-= (const Grid<DIM,S>& a) {
    gridSub<DIM,T,S> (*this, a);
    return *this;
}
template<int DIM, class T> template<class S> Grid<DIM,T>& Grid<DIM,T>::operator-= (const S& a) {
    gridAddScalar<DIM,T,S> (*this, -a);
    return *this;
}
template<int DIM, class T> template<class S> Grid<DIM,T>& Grid<DIM,T>::operator*= (const Grid<DIM,S>& a) {
    gridMult<DIM,T,S> (*this, a);
    return *this;
}
template<int DIM, class T> Grid<DIM,T>& Grid<DIM,T>::operator*= (const T& a) {
    gridMultScalar<DIM,T,T> (*this, a);
    return *this;
}
template<int DIM, class T> template<class S> Grid<DIM,T>& Grid<DIM,T>::operator/= (const Grid<DIM,S>& a) {
    gridDiv<DIM,T,S> (*this, a);
    return *this;
}
template<int DIM, class T> template<class S> Grid<DIM,T>& Grid<DIM,T>::operator/= (const S& a) {
    S rez((S)1.0 / a);
    gridMultScalar<DIM,T,S> (*this, rez);
    return *this;
}

} //namespace
#endif