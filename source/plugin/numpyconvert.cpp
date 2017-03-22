/******************************************************************************
*
* MantaFlow fluid solver framework
* Copyright 2017 Steffen Wiewel, Moritz Baecher 
*
* This program is free software, distributed under the terms of the
* GNU General Public License (GPL) 
* http://www.gnu.org/licenses
*
* Plugins to convert mantaflow grids to/from numpy arrays, also support pdata fields  
# (only compiled if NUMPY is enabled)
*
******************************************************************************/

#include "manta.h"
#include "kernel.h"
#include "grid.h"
#include "particle.h"
#include "levelset.h"

using namespace std;
namespace Manta
{

//====================================================================================================
// Grid numpy conversion
//----------------------------------------------------------------------------------------------------

template<typename T>
void copyArrayToGridScalar(const PyArrayContainer _Source, T& _Target)
{
	_Target.setConst(0.0f);
	unsigned int uGridSize = _Target.getSizeX() * _Target.getSizeY() * _Target.getSizeZ();
	assertMsg(_Source.TotalSize == uGridSize, "The size of the numpy array doesn't match the size of the Grid!");
	
	NumpyTypes eDataType  = _Source.DataType; 

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(_Target) { _Target(idx) = (reinterpret_cast<float*>(_Source.pData))[idx]; }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(_Target) { _Target(idx) = (reinterpret_cast<double*>(_Source.pData))[idx]; } 
			break;
		default:
			errMsg("unknown/unsupported type of Numpy array");
			return;
	}
}

template<typename T>
void copyGridToArrayScalar(const T& _Source, PyArrayContainer _Target)
{
	unsigned int uGridsize = _Source.getSizeX() * _Source.getSizeY() * _Source.getSizeZ();
	assertMsg(_Target.TotalSize == uGridsize, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = _Target.DataType;

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IDX(_Source) { reinterpret_cast<float*>(_Target.pData)[idx] = _Source(idx); }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IDX(_Source) { reinterpret_cast<double*>(_Target.pData)[idx] = _Source(idx); }
			break;
		default:
			break;
	}
}

template<typename T>
void copyArrayToGridVector(const PyArrayContainer _Source, T& _Target)
{
	unsigned int uSizeX = _Target.getSizeX();
	unsigned int uSizeY = _Target.getSizeY();
	unsigned int uSizeZ = _Target.getSizeZ();
	unsigned int uSizeW = 3u;
	
	assertMsg(_Source.TotalSize == uSizeX * uSizeY * uSizeZ * uSizeW, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = _Source.DataType;

	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IJK(_Target) { for(int w = 0; w < 3; ++w) {
				_Target(i,j,k).value[w] = reinterpret_cast<float*>(_Source.pData)[w + uSizeW * (k +  uSizeZ * (i + uSizeX * j))]; } }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IJK(_Target) { for(int w = 0; w < 3; ++w) {
				_Target(i,j,k).value[w] = reinterpret_cast<double*>(_Source.pData)[w + uSizeW * (k +  uSizeZ * (i + uSizeX * j))]; } }
			break;
		default:
			break;
	}
}

template<typename T>
void copyGridToArrayVector(const T& _Source, PyArrayContainer _Target)
{
	unsigned int uSizeX = _Source.getSizeX();
	unsigned int uSizeY = _Source.getSizeY();
	unsigned int uSizeZ = _Source.getSizeZ();
	unsigned int uSizeW = 3u;

	assertMsg(_Target.TotalSize == uSizeX * uSizeY * uSizeZ * uSizeW, "The size of the numpy array doesn't match the size of the grid!");
	
	NumpyTypes eDataType = _Target.DataType;
	
	switch (eDataType)
	{
		case NumpyTypes::N_FLOAT:
			FOR_IJK(_Source) { for(int w = 0; w < 3; ++w) {
				reinterpret_cast<float*>(_Target.pData)[w + uSizeW * (k +  uSizeZ * (i + uSizeX * j))] = _Source(i,j,k).value[w]; } }
			break;
		case NumpyTypes::N_DOUBLE:
			FOR_IJK(_Source) { for(int w = 0; w < 3; ++w) {
				reinterpret_cast<double*>(_Target.pData)[w + uSizeW * (k +  uSizeZ * (i + uSizeX * j))] = _Source(i,j,k).value[w]; } }
			break;
		default:
			break;
	}
}

//====================================================================================================
// Python interface
//----------------------------------------------------------------------------------------------------

PYTHON() void copyArrayToGridReal(const PyArrayContainer source, Grid<Real>& target) {
	copyArrayToGridScalar<Grid<Real>>(source, target);
}

PYTHON() void copyGridToArrayReal(const Grid<Real>& source, PyArrayContainer target) {
	copyGridToArrayScalar<Grid<Real>>(source, target);
}

PYTHON() void copyArrayToGridLevelset(const PyArrayContainer source, LevelsetGrid& target) {
	copyArrayToGridScalar<LevelsetGrid>(source, target);
}

PYTHON() void copyGridToArrayLevelset(const LevelsetGrid& source, PyArrayContainer target) {
	copyGridToArrayScalar<LevelsetGrid>(source, target);
}

PYTHON() void copyArrayToGridVec3(const PyArrayContainer source, Grid<Vec3>& target) {
	copyArrayToGridVector<Grid<Vec3>>(source, target);
}

PYTHON() void copyGridToArrayVec3(const Grid<Vec3>& source, PyArrayContainer target) {
	copyGridToArrayVector<Grid<Vec3>>(source, target);
}

PYTHON() void copyArrayToGridMAC(const PyArrayContainer source, MACGrid& target) {
	copyArrayToGridVector<MACGrid>(source, target);
}

PYTHON() void copyGridToArrayMAC(const MACGrid& _Source, PyArrayContainer target) {
	copyGridToArrayVector<MACGrid>(_Source, target);
}

//====================================================================================================
// pdata conversion functions
//----------------------------------------------------------------------------------------------------

template<typename T>
void numpyToParticleDataImpl(ParticleDataImpl<T> &p, const PyArrayContainer n) {
	assertMsg(n.TotalSize == p.size(), "Sizes are different!");
	std::copy(reinterpret_cast<const T*>(n.pData), reinterpret_cast<const T*>(n.pData)+n.TotalSize,  &(p[0]));
}
template<typename T>
void particleDataImplToNumpy(PyArrayContainer n, const ParticleDataImpl<T> &p) {
	assertMsg(n.TotalSize == p.size(), "Sizes are different!");
	std::copy(&(p[0]), &(p[0])+n.TotalSize, reinterpret_cast<T*>(n.pData));
}

// python interface

PYTHON() void copyArrayToPdataInt(ParticleDataImpl<int> &p, const PyArrayContainer n) { 
	numpyToParticleDataImpl<int>(p, n); 
}
PYTHON() void copyPdataToArrayInt(PyArrayContainer n, const ParticleDataImpl<int> &p) { 
	particleDataImplToNumpy<int>(n, p); 
}

PYTHON() void copyArrayToPdataReal(ParticleDataImpl<Real> &p, const PyArrayContainer n) { 
	numpyToParticleDataImpl<Real>(p, n); 
}
PYTHON() void copyPdataToArrayReal(PyArrayContainer n, const ParticleDataImpl<Real> &p) { 
	particleDataImplToNumpy<Real>(n, p); 
}

PYTHON() void copyArrayToPdataVec3(ParticleDataImpl<Vec3> &p, const PyArrayContainer n) {
	assertMsg(n.TotalSize == p.size()*3, "Sizes are different!");
	std::copy(reinterpret_cast<const Real*>(n.pData), reinterpret_cast<const Real*>(n.pData)+n.TotalSize,  &(p[0][0]));
}
PYTHON() void copyPdataToArrayVec3(PyArrayContainer n, const ParticleDataImpl<Vec3> &p) {
	assertMsg(n.TotalSize == p.size()*3, "Sizes are different!");
	std::copy(&(p[0][0]), &(p[0][0])+n.TotalSize, reinterpret_cast<Real*>(n.pData));
}

} // manta

