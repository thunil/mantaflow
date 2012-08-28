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

#include "grid.h"
#include "levelset.h"
#include "kernel.h"
#include <limits>
#include <sstream>
#include <cstring>
#include "fileio.h"

using namespace std;
namespace Manta {

//******************************************************************************
// GridBase members

GridBase::GridBase (FluidSolver* parent) 
    : PbClass(parent), mType(TypeNone)
{
    checkParent();
    mDim = getParent()->is2D() ? 2:3;
}

void GridBase::checkIndex(int i, int j, int k) const {
    if (i<0 || j<0 || k<0 || i>=mSize.x || j>=mSize.y || k>= mSize.z) {
        std::stringstream s;
        s << "Grid " << mName << " dim " << mSize << " : index " << i << "," << j << "," << k << " out of bound ";
        errMsg(s.str());
    }
}

void GridBase::checkIndex(int idx) const {
    if (idx<0 || idx > mSize.x * mSize.y * mSize.z) {
        std::stringstream s;
        s << "Grid " << mName << " dim " << mSize << " : index " << idx << " out of bound ";
        errMsg(s.str());
    }
}


//******************************************************************************
// Grid<T> members

// helpers to set type
template<class T> inline GridBase::GridType typeList() { return GridBase::TypeNone; }
template<> inline GridBase::GridType typeList<Real>() { return GridBase::TypeReal; }
template<> inline GridBase::GridType typeList<int>() { return GridBase::TypeInt; }
template<> inline GridBase::GridType typeList<Vec3>() { return GridBase::TypeVec3; }

template<class T>
Grid<T>::Grid(FluidSolver* parent, bool show)
    : GridBase(parent)
{
    mType = typeList<T>();
    mSize = parent->getGridSize();
    mData = parent->getGridPointer<T>();
    
    mStrideZ = parent->is2D() ? 0 : (mSize.x * mSize.y);
    mDx = 1.0 / mSize.max();
    clear();
    setHidden(!show);
}

template<class T>
Grid<T>::~Grid() {
    mParent->freeGridPointer<T>(mData);    
}

template<class T>
void Grid<T>::clear() {
    memset(mData, 0, sizeof(T) * mSize.x * mSize.y * mSize.z);    
}

template<class T>
void Grid<T>::swap(Grid<T>& other) {
    if (other.getSizeX() != getSizeX() || other.getSizeY() != getSizeY() || other.getSizeZ() != getSizeZ())
        errMsg("Grid::swap(): Grid dimensions mismatch.");
    
    T* dswap = other.mData;
    other.mData = mData;
    mData = dswap;
}

template<class T>
void Grid<T>::save(string name) {
    if (name.find_last_of('.') == string::npos)
        errMsg("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if (ext == ".raw")
        writeGridRaw(name, this);
    else if (ext == ".uni")
        writeGridUni(name, this);
    else
        errMsg("file '" + name +"' filetype not supported");
}

//******************************************************************************
// Grid<T> operators

//! Kernel: Compute minmax value of Real grid
KERNEL(idx, reduce) struct CompMinmaxReal (Grid<Real>& val) {
    Real minVal, maxVal;
    
    void operator()(int idx) {
        if (val[idx] < minVal)
            minVal = val[idx];
        if (val[idx] > maxVal)
            maxVal = val[idx];
    }
    void setup() {
        minVal = std::numeric_limits<Real>::max();
        maxVal = -std::numeric_limits<Real>::max();
    }
    void join(const CompMinmaxReal& a) {
        minVal = std::min(a.minVal, minVal);
        maxVal = std::max(a.maxVal, maxVal);
    }
};

//! Kernel: Compute minmax value of int grid
KERNEL(idx, reduce) struct CompMinmaxInt (Grid<int>& val) {
    int minVal, maxVal;
    
    void operator()(int idx) {
        if (val[idx] < minVal)
            minVal = val[idx];
        if (val[idx] > maxVal)
            maxVal = val[idx];
    }
    void setup() {
        minVal = std::numeric_limits<int>::max();
        maxVal = std::numeric_limits<int>::min();
    }
    void join(const CompMinmaxInt& a) {
        minVal = std::min(a.minVal, minVal);
        maxVal = std::max(a.maxVal, maxVal);
    }
};

//! Kernel: Compute minmax squared norm of Vec3 grid
KERNEL(idx, reduce) struct CompMinmaxVec3 (Grid<Vec3>& val) {
    Real minVal2, maxVal2;
    
    void operator()(int idx) {
        const Real s = normSquare(val[idx]);
        if (s < minVal2)
            minVal2 = s;
        if (s > maxVal2)
            maxVal2 = s;
    }
    void setup() {
        minVal2 = std::numeric_limits<Real>::max();
        maxVal2 = 0;
    }
    void join(const CompMinmaxVec3& a) {
        minVal2 = std::min(a.minVal2, minVal2);
        maxVal2 = std::max(a.maxVal2, maxVal2);
    }
};

template<class T> Grid<T>& Grid<T>::safeDivide (const Grid<T>& a) {
    gridSafeDiv<T> (*this, a);
    return *this;
}
template<class T> Grid<T>& Grid<T>::operator= (const Grid<T>& a) {
    memcpy(mData, a.mData, sizeof(T) * mSize.x * mSize.y * mSize.z);
    mType = a.mType; // copy type marker
    return *this;
}
template<class T> Grid<T>& Grid<T>::operator= (const T& a) {
    FOR_IDX(*this) { mData[idx] = a; }
    return *this;
}
template<class T> void Grid<T>::scaledAdd(const Grid<T>& a, const T& factor) {
    gridScaleAdd<T> (*this, a, factor);
}
template<> Real Grid<Real>::getMaxValue() {
    return CompMinmaxReal (*this).maxVal;
}
template<> Real Grid<Real>::getMinValue() {
    return CompMinmaxReal (*this).minVal;
}
template<> Real Grid<Real>::getMaxAbsValue() {
    CompMinmaxReal op (*this);
    return max( fabs(op.minVal), fabs(op.maxVal));
}
template<> Real Grid<Vec3>::getMaxValue() {
    return sqrt(CompMinmaxVec3 (*this).maxVal2);
}
template<> Real Grid<Vec3>::getMinValue() { 
    return sqrt(CompMinmaxVec3 (*this).minVal2);
}
template<> Real Grid<Vec3>::getMaxAbsValue() {
    return sqrt(CompMinmaxVec3 (*this).maxVal2);
}
template<> Real Grid<int>::getMaxValue() {
    return (Real) CompMinmaxInt (*this).maxVal;
}
template<> Real Grid<int>::getMinValue() {
    return (Real) CompMinmaxInt (*this).minVal;
}
template<> Real Grid<int>::getMaxAbsValue() {
    CompMinmaxInt op (*this);
    return max( fabs((Real)op.minVal), fabs((Real)op.maxVal));
}
template<class T> void Grid<T>::add(const Grid<T>& a, const Grid<T>& b) {
    gridAdd2<T>(*this, a, b);
}

//******************************************************************************
// Specialization classes

void FlagGrid::initDomain(int boundaryWidth) {
    memset(mData, TypeEmpty, sizeof(int) * mSize.x * mSize.y * mSize.z);    
    initBoundaries(boundaryWidth);
}

void FlagGrid::initBoundaries(int boundaryWidth) {
    const int w = boundaryWidth;
    FOR_IJK(*this) {
        bool bnd = (i<=w || i>=mSize.x-1-w || j<=w || j>=mSize.y-1-w || (is3D() && (k<=w || k>=mSize.z-1-w)));
        if (bnd) 
            mData[index(i,j,k)] = TypeObstacle;
    }
}

void FlagGrid::updateFromLevelset(LevelsetGrid& levelset) {
    FOR_IDX(*this) {
        if (!isObstacle(idx)) {
            const Real phi = levelset[idx];
            if (phi <= levelset.invalidTimeValue()) continue;
            
            mData[idx] &= ~(TypeEmpty | TypeFluid); // clear empty/fluid flags
            mData[idx] |= (phi <= 0) ? TypeFluid : TypeEmpty; // set resepctive flag
        }
    }
}   

void FlagGrid::fillGrid() {
    FOR_IDX(*this) {
        if ((mData[idx] & TypeObstacle)==0)
            mData[idx] = (mData[idx] & ~TypeEmpty) | TypeFluid;
    }
}

// explicit instantiation
template class Grid<int>;
template class Grid<Real>;
template class Grid<Vec3>;

} //namespace