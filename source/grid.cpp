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
}

void GridBase::checkIndex(int i, int j, int k) const {
    if (i<0 || j<0 || k<0 || i>=mSize.x || j>=mSize.y || k>= mSize.z) {
        std::stringstream s;
        s << "Grid " << mName << " dim " << mSize << " : index " << i << "," << j << "," << k << " out of bound ";
        throw Error(s.str());
    }
}

void GridBase::checkIndex(int idx) const {
    if (idx<0 || idx > mSize.x * mSize.y * mSize.z) {
        std::stringstream s;
        s << "Grid " << mName << " dim " << mSize << " : index " << idx << " out of bound ";
        throw Error(s.str());
    }
}


//******************************************************************************
// Grid<DIM,T> members

// helpers to set type
template<class T> inline GridBase::GridType typeList() { return GridBase::TypeNone; }
template<> inline GridBase::GridType typeList<Real>() { return GridBase::TypeReal; }
template<> inline GridBase::GridType typeList<int>() { return GridBase::TypeInt; }
template<> inline GridBase::GridType typeList<Vec3>() { return GridBase::TypeVec3; }

template<int DIM,class T>
Grid<DIM,T>::Grid(FluidSolver* parent, bool show)
    : GridBase(parent)
{     
    mType = typeList<T>();
    mSize = parent->getGridSize();
    mData = parent->getGridPointer<T>();
    
    mStrideZ = mSize.x * mSize.y;
    mDx = 1.0 / mSize.max();
    clear();
    setHidden(!show);
}

template<int DIM,class T>
Grid<DIM,T>::~Grid() {
    mParent->freeGridPointer<T>(mData);    
}

template<int DIM,class T>
void Grid<DIM,T>::clear() {
    memset(mData, 0, sizeof(T) * mSize.x * mSize.y * mSize.z);    
}

template<int DIM,class T>
void Grid<DIM,T>::swap(Grid<DIM,T>& other) {
    if (other.getSizeX() != getSizeX() || other.getSizeY() != getSizeY() || other.getSizeZ() != getSizeZ())
        throw Error("Grid::swap(): Grid dimensions mismatch.");
    
    T* dswap = other.mData;
    other.mData = mData;
    mData = dswap;
}

template<int DIM,class T>
void Grid<DIM,T>::save(string name) {
    if (name.find_last_of('.') == string::npos)
        throw Error("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if (ext == ".raw")
        writeGridRaw(name, this);
    else if (ext == ".uni")
        writeGridUni(name, this);
    else
        throw Error("file '" + name +"' filetype not supported");
}

//******************************************************************************
// Grid<DIM,T> operators

//! Kernel: Compute minmax value of Real grid
KERNEL(idx, reduce) template<int DIM> 
struct CompMinmaxReal (Grid<DIM, Real>& val) {
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
KERNEL(idx, reduce) template<int DIM> 
struct CompMinmaxInt (Grid<DIM,int>& val) {
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
template<int DIM> KERNEL(idx, reduce)
struct CompMinmaxVec3 (Grid<DIM,Vec3>& val) {
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

template<int DIM,class T> Grid<DIM,T>& Grid<DIM,T>::safeDivide (const Grid<DIM,T>& a) {
    gridSafeDiv<DIM,T> (*this, a);
    return *this;
}
template<int DIM,class T> Grid<DIM,T>& Grid<DIM,T>::operator= (const Grid<DIM,T>& a) {
    memcpy(mData, a.mData, sizeof(T) * mSize.x * mSize.y * mSize.z);
    mType = a.mType; // copy type marker
    return *this;
}
template<int DIM,class T> Grid<DIM,T>& Grid<DIM,T>::operator= (const T& a) {
    FOR_IDX(*this) { mData[idx] = a; }
    return *this;
}
template<int DIM,class T> void Grid<DIM,T>::scaledAdd(const Grid<DIM,T>& a, const T& factor) {
    gridScaleAdd<DIM,T> (*this, a, factor);
}
template<int DIM, class T> void Grid<DIM,T>::add(const Grid<DIM,T>& a, const Grid<DIM,T>& b) {
    gridAdd2<DIM,T>(*this, a, b);
}

template<> Real Grid<3,Real>::getMaxValue() { return CompMinmaxReal<3> (*this).maxVal; }
template<> Real Grid<2,Real>::getMaxValue() { return CompMinmaxReal<2> (*this).maxVal; }
template<> Real Grid<3,Real>::getMinValue() { return CompMinmaxReal<3> (*this).minVal; }
template<> Real Grid<2,Real>::getMinValue() { return CompMinmaxReal<2> (*this).minVal; }
template<> Real Grid<3,int>::getMaxValue() { return (Real) CompMinmaxInt<3> (*this).maxVal; }
template<> Real Grid<2,int>::getMaxValue() { return (Real) CompMinmaxInt<2> (*this).maxVal; }
template<> Real Grid<3,int>::getMinValue() { return (Real) CompMinmaxInt<3> (*this).minVal; }
template<> Real Grid<2,int>::getMinValue() { return (Real) CompMinmaxInt<2> (*this).minVal; }
template<> Real Grid<3,Vec3>::getMaxValue() { return sqrt(CompMinmaxVec3<3> (*this).maxVal2); }
template<> Real Grid<2,Vec3>::getMaxValue() { return sqrt(CompMinmaxVec3<2> (*this).maxVal2); }
template<> Real Grid<3,Vec3>::getMinValue() { return sqrt(CompMinmaxVec3<3> (*this).minVal2); }
template<> Real Grid<2,Vec3>::getMinValue() { return sqrt(CompMinmaxVec3<2> (*this).minVal2); }
template<> Real Grid<3,Vec3>::getMaxAbsValue() { return sqrt(CompMinmaxVec3<3> (*this).maxVal2); }
template<> Real Grid<2,Vec3>::getMaxAbsValue() { return sqrt(CompMinmaxVec3<2> (*this).maxVal2); }
template<> Real Grid<3,Real>::getMaxAbsValue() { CompMinmaxReal<3> op (*this); return max( fabs(op.minVal), fabs(op.maxVal)); }
template<> Real Grid<2,Real>::getMaxAbsValue() { CompMinmaxReal<2> op (*this); return max( fabs(op.minVal), fabs(op.maxVal)); }
template<> Real Grid<3,int>::getMaxAbsValue() { CompMinmaxInt<3> op (*this); return max( fabs((Real)op.minVal), fabs((Real)op.maxVal)); }
template<> Real Grid<2,int>::getMaxAbsValue() { CompMinmaxInt<2> op (*this); return max( fabs((Real)op.minVal), fabs((Real)op.maxVal)); }

//******************************************************************************
// Specialization classes

template<int DIM>
void FlagGrid<DIM>::initDomain(int boundaryWidth) {
    memset(this->mData, TypeEmpty, sizeof(int) * this->mSize.x * this->mSize.y * this->mSize.z);    
    initBoundaries(boundaryWidth);
}

template<int DIM>
void FlagGrid<DIM>::initBoundaries(int boundaryWidth) {
    FOR_IJK(*this) {
        if (!this->isInBounds(Vec3i(i,j,k),boundaryWidth))
            (*this)(i,j,k) = TypeObstacle;
    }
}

KERNEL(idx) template<int DIM>
KnUpdateFromLevelset(FlagGrid<DIM>& flags, LevelsetGrid<DIM>& levelset) {
    if (!flags.isObstacle(idx)) {
        const Real phi = levelset[idx];
        if (phi <= levelset.invalidTimeValue()) continue;
        
        flags[idx] &= ~(flags.TypeEmpty | flags.TypeFluid); // clear empty/fluid flags
        flags[idx] |= (phi <= 0) ? flags.TypeFluid : flags.TypeEmpty; // set resepctive flag
    }
}

template<int DIM>
void FlagGrid<DIM>::updateFromLevelset(LevelsetGrid<DIM>& levelset) {
    KnUpdateFromLevelset<DIM>(*this,levelset);
}   

template<int DIM>
void FlagGrid<DIM>::fillGrid() {
    FOR_IDX(*this) {
        if (((*this)[idx] & TypeObstacle)==0)
            (*this)[idx] = ((*this)[idx] & ~TypeEmpty) | TypeFluid;
    }
}

// explicit instantiation
template class Grid<3,int>;
template class Grid<3,Real>;
template class Grid<3,Vec3>;
template class Grid<2,int>;
template class Grid<2,Real>;
template class Grid<2,Vec3>;
template class FlagGrid<2>;
template class FlagGrid<3>;

} //namespace