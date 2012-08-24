/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Loading and writing grids and meshes to disk
 *
 ******************************************************************************/

#ifndef _FILEIO_H
#define _FILEIO_H

#include <string>

namespace Manta {
// forward decl.
template<int DIM, class T> class Grid;
class Mesh;

void writeObjFile(const std::string& name, Mesh* mesh);
void writeBobjFile(const std::string& name, Mesh* mesh);
void readObjFile(const std::string& name, Mesh* mesh, bool append);
template<int DIM, class T> void writeGridRaw(const std::string& name, Grid<DIM,T>* grid);
template<int DIM, class T> void writeGridUni(const std::string& name, Grid<DIM,T>* grid);

} // namespace

#endif