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
    m3D = getParent()->is3D();
}

//******************************************************************************
// Grid<T> members

// helpers to set type
template<class T> inline GridBase::GridType typeList() { return GridBase::TypeNone; }
template<> inline GridBase::GridType typeList<Real>()  { return GridBase::TypeReal; }
template<> inline GridBase::GridType typeList<int>()   { return GridBase::TypeInt;  }
template<> inline GridBase::GridType typeList<Vec3>()  { return GridBase::TypeVec3; }

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
Grid<T>::Grid(const Grid<T>& a) : GridBase(a.getParent()) {
    mSize = a.mSize;
    mType = a.mType;
    mStrideZ = a.mStrideZ;
    mDx = a.mDx;
    FluidSolver *gp = a.getParent();
    mData = gp->getGridPointer<T>();
    memcpy(mData, a.mData, sizeof(T) * a.mSize.x * a.mSize.y * a.mSize.z);
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
void Grid<T>::load(string name) {
    if (name.find_last_of('.') == string::npos)
        errMsg("file '" + name + "' does not have an extension");
    string ext = name.substr(name.find_last_of('.'));
    if (ext == ".raw")
        readGridRaw(name, this);
    else if (ext == ".uni")
        readGridUni(name, this);
    else
        errMsg("file '" + name +"' filetype not supported");
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
	else if (ext == ".vol")
		writeGridVol(name, this);
	else if (ext == ".txt")
		writeGridTxt(name, this);
    else
        errMsg("file '" + name +"' filetype not supported");
}

//******************************************************************************
// Grid<T> operators

//! Kernel: Compute min value of Real grid
KERNEL(idx, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real CompMinReal(Grid<Real>& val) {
    if (val[idx] < minVal)
        minVal = val[idx];
}

//! Kernel: Compute max value of Real grid
KERNEL(idx, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
Real CompMaxReal(Grid<Real>& val) {
    if (val[idx] > maxVal)
        maxVal = val[idx];
}

//! Kernel: Compute min value of int grid
KERNEL(idx, reduce=min) returns(int minVal=std::numeric_limits<int>::max())
int CompMinInt(Grid<int>& val) {
    if (val[idx] < minVal)
        minVal = val[idx];
}

//! Kernel: Compute max value of int grid
KERNEL(idx, reduce=max) returns(int maxVal=-std::numeric_limits<int>::min())
int CompMaxInt(Grid<int>& val) {
    if (val[idx] > maxVal)
        maxVal = val[idx];
}

//! Kernel: Compute min norm of vec grid
KERNEL(idx, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real CompMinVec(Grid<Vec3>& val) {
    const Real s = normSquare(val[idx]);
    if (s < minVal)
        minVal = s;
}

//! Kernel: Compute max norm of vec grid
KERNEL(idx, reduce=max) returns(Real maxVal=0)
Real CompMaxVec(Grid<Vec3>& val) {
    const Real s = normSquare(val[idx]);
    if (s > maxVal)
        maxVal = s;
}


template<class T> Grid<T>& Grid<T>::safeDivide (const Grid<T>& a) {
    gridSafeDiv<T> (*this, a);
    return *this;
}
template<class T> Grid<T>& Grid<T>::operator= (const Grid<T>& a) {
    assertMsg (a.mSize.x == mSize.x && a.mSize.y == mSize.y && a.mSize.z == mSize.z, "different grid resolutions");
    memcpy(mData, a.mData, sizeof(T) * mSize.x * mSize.y * mSize.z);
    mType = a.mType; // copy type marker
    return *this;
}
/*template<class T> Grid<T>& Grid<T>::operator= (const T& a) {
    FOR_IDX(*this) { mData[idx] = a; }
    return *this;
}*/

//! Grid a += b*factor (note, shouldnt be part of the grid class! can cause problems with python instantiation)
// also the python integration doesnt support templated functions for now (only classes)
PYTHON void scaledAddReal(Grid<Real>& a, const Grid<Real>& b, const Real& factor) {
    gridScaledAdd<Real,Real> (a, b, factor);
}
PYTHON void scaledAddVec3(Grid<Vec3>& a, const Grid<Vec3>& b, const Vec3& factor) {
    gridScaledAdd<Vec3,Vec3> (a, b, factor);
}
template<> Real Grid<Real>::getMaxValue() {
    return CompMaxReal (*this);
}
template<> Real Grid<Real>::getMinValue() {
    return CompMinReal (*this);
}
template<> Real Grid<Real>::getMaxAbsValue() {
    Real amin = CompMinReal (*this);
    Real amax = CompMaxReal (*this);
    return max( fabs(amin), fabs(amax));
}
template<> Real Grid<Vec3>::getMaxValue() {
    return sqrt(CompMaxVec (*this));
}
template<> Real Grid<Vec3>::getMinValue() { 
    return sqrt(CompMinVec (*this));
}
template<> Real Grid<Vec3>::getMaxAbsValue() {
    return sqrt(CompMaxVec (*this));
}
template<> Real Grid<int>::getMaxValue() {
    return (Real) CompMaxInt (*this);
}
template<> Real Grid<int>::getMinValue() {
    return (Real) CompMinInt (*this);
}
template<> Real Grid<int>::getMaxAbsValue() {
    int amin = CompMinInt (*this);
    int amax = CompMaxInt (*this);
    return max( fabs((Real)amin), fabs((Real)amax));
}
template<class T> void Grid<T>::add(const Grid<T>& a, const Grid<T>& b) {
    gridAdd2<T>(*this, a, b);
}

PYTHON void setConstant    (Grid<Real>& grid, Real value=0.) { gridSetConst<Real>(grid,value); }
PYTHON void setConstantVec3(Grid<Vec3>& grid, Vec3 value=0.) { gridSetConst<Vec3>(grid,value); }
PYTHON void setConstantInt (Grid<int >& grid, int  value=0.) { gridSetConst<int>(grid,value); }

// compute maximal diference of two cells in the grid
// used for testing
PYTHON Real gridMaxDiff(Grid<Real>& g1, Grid<Real>& g2 )
{
	double maxVal = 0.;
    FOR_IJK(g1) {
		maxVal = std::max(maxVal, (double)fabs( g1(i,j,k)-g2(i,j,k) ));
	}
	return maxVal; 
}
PYTHON Real gridMaxDiffVec3(Grid<Vec3>& g1, Grid<Vec3>& g2 )
{
	double maxVal = 0.;
    FOR_IJK(g1) {
		maxVal = std::max(maxVal, (double)fabs( norm(g1(i,j,k)-g2(i,j,k)) ));
	}
	return maxVal; 
}

// simple helper functions to convert mac to vec3 , and levelset to real grids
// (are assumed to be the same for running the test cases - in general they're not!)
PYTHON void convertMacToVec3 (MACGrid &source, Grid<Vec3>& target)
{
    FOR_IJK(target) {
		target(i,j,k) = source(i,j,k);
	}
}

PYTHON void convertLevelsetToReal (LevelsetGrid &source , Grid<Real> &target)
{
    FOR_IJK(target) {
		target(i,j,k) = source(i,j,k);
	}
}


//******************************************************************************
// Specialization classes

void FlagGrid::initDomain(int boundaryWidth) {
    FOR_IDX(*this)
        mData[idx] = TypeEmpty;
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

void FlagGrid::fillGrid(int type) {
    FOR_IDX(*this) {
        if ((mData[idx] & TypeObstacle)==0)
            mData[idx] = (mData[idx] & ~(TypeEmpty | TypeFluid)) | type;
    }
}

// explicit instantiation
template class Grid<int>;
template class Grid<Real>;
template class Grid<Vec3>;

//template void scaledAdd<Real,Real>(const Grid<Real>& a, const Grid<Real>& b, const Real& factor);

#if ENABLE_GRID_TEST_DATATYPE==1
// instantiate test datatype , not really required for simulations, mostly here for demonstration purposes
template class Grid<nbVector>;
#endif // ENABLE_GRID_TEST_DATATYPE


} //namespace
