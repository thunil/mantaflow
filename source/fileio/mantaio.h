/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Loading and writing grids and meshes to disk
 *
 ******************************************************************************/

#ifndef _FILEIO_H
#define _FILEIO_H

#include <string>

#include "manta.h"

namespace Manta {

// Forward declations
class Mesh;
class FlagGrid;
class GridBase;
template<class T> class Grid;
template<class T> class Grid4d;
class BasicParticleSystem;
template<class T> class ParticleDataImpl;
template<class T> class MeshDataImpl;

// Obj format
void writeObjFile(const std::string& name, Mesh* mesh);
void writeBobjFile(const std::string& name, Mesh* mesh);
void readObjFile(const std::string& name, Mesh* mesh, bool append);
void readBobjFile(const std::string& name, Mesh* mesh, bool append);

// Other formats (Raw, Uni, Vol)
template<class T> void readGridUni (const std::string& name, Grid<T>* grid);
template<class T> void readGridRaw (const std::string& name, Grid<T>* grid);
template<class T> void readGridVol (const std::string& name, Grid<T>* grid);
int readGridsRaw(const std::string& name, std::vector<PbClass*>* grids);
int readGridsUni(const std::string& name, std::vector<PbClass*>* grids);
int readGridsVol(const std::string& name, std::vector<PbClass*>* grids);
int readGridsTxt(const std::string& name, std::vector<PbClass*>* grids);

template<class T> void writeGridRaw(const std::string& name, Grid<T>* grid);
template<class T> void writeGridUni(const std::string& name, Grid<T>* grid);
template<class T> void writeGridVol(const std::string& name, Grid<T>* grid);
template<class T> void writeGridTxt(const std::string& name, Grid<T>* grid);
int writeGridsRaw(const std::string& name, std::vector<PbClass*>* grids);
int writeGridsUni(const std::string& name, std::vector<PbClass*>* grids);
int writeGridsVol(const std::string& name, std::vector<PbClass*>* grids);
int writeGridsTxt(const std::string& name, std::vector<PbClass*>* grids);

// OpenVDB
#if OPENVDB==1
int writeGridsVDB(const std::string& filename, std::vector<PbClass*>* grids, float scale=1.0);
int readGridsVDB(const std::string& filename, std::vector<PbClass*>* grids, float scale=1.0);
#endif

// Numpy
template<class T> void writeGridNumpy(const std::string& name, Grid<T>* grid);
template<class T> void readGridNumpy (const std::string& name, Grid<T>* grid);

int writeGridsNumpy(const std::string& name, std::vector<PbClass*>* grids);
int readGridsNumpy(const std::string& name, std::vector<PbClass*>* grids);

// 4D Grids
template<class T> void writeGrid4dUni(const std::string& name, Grid4d<T>* grid);
template<class T> void readGrid4dUni (const std::string& name, Grid4d<T>* grid, int readTslice=-1, Grid4d<T>* slice=NULL, void** fileHandle=NULL);
void readGrid4dUniCleanup(void** fileHandle);
template<class T> void writeGrid4dRaw(const std::string& name, Grid4d<T>* grid);
template<class T> void readGrid4dRaw (const std::string& name, Grid4d<T>* grid);

// Particles + particle data
void writeParticlesUni(const std::string& name, const BasicParticleSystem* parts );
void readParticlesUni (const std::string& name, BasicParticleSystem* parts );

template <class T> void writePdataUni(const std::string& name, ParticleDataImpl<T>* pdata );
template <class T> void readPdataUni (const std::string& name, ParticleDataImpl<T>* pdata );

// Mesh data
template <class T> void writeMdataUni(const std::string& name, MeshDataImpl<T>* mdata );
template <class T> void readMdataUni (const std::string& name, MeshDataImpl<T>* mdata );

// Helpers
void getUniFileSize(const std::string& name, int& x, int& y, int& z, int* t = NULL, std::string* info = NULL);
void *safeGzopen(const char *filename, const char *mode);
#if OPENVDB==1
template<class S, class T> void convertFrom(S& in, T* out);
template<class S, class T> void convertTo(S* out, T& in);
#endif

} // namespace

#endif
