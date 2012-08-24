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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <zlib.h>
#include "fileio.h"
#include "grid.h"
#include "mesh.h"
#include "vortexsheet.h"

using namespace std;

namespace Manta {
/*
void trim(string& s) {
    size_t pl = s.find_first_not_of("\t\r ");
    if (pl != string::npos) s=s.substr(pl);
    size_t pr = s.find_last_not_of("\t\r ");
    if (pr != string::npos) s=s.substr(0,pr);
}*/

void writeObjFile(const string& name, Mesh* mesh) {
    throw Error("obj exporter not yet implemented");
}

void writeBobjFile(const string& name, Mesh* mesh) {
    cout << "writing mesh file " << name << endl;
    const Real dx = mesh->getParent()->getDx();
    const Vec3i gs = mesh->getParent()->getGridSize();
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf)
        throw Error("writeBobj: unable to open file");
    
    // write vertices
    int numVerts = mesh->numNodes();
    gzwrite(gzf, &numVerts, sizeof(int));
    for (int i=0; i<numVerts; i++) {
        Vector3D<float> pos = toVec3f(mesh->nodes(i).pos);
        // normalize
        pos -= toVec3f(gs)*0.5;
        pos *= dx;
        gzwrite(gzf, &pos.value[0], sizeof(float)*3);
    }
    
    // normals
    mesh->computeVertexNormals();
    gzwrite(gzf, &numVerts, sizeof(int));
    for (int i=0; i<numVerts; i++) {
        Vector3D<float> pos = toVec3f(mesh->nodes(i).normal);
        gzwrite(gzf, &pos.value[0], sizeof(float)*3);
    }
    
    // write tris
    int numTris = mesh->numTris();
    gzwrite(gzf, &numTris, sizeof(int));
    for(int t=0; t<numTris; t++) {
        for(int j=0; j<3; j++) { 
            int trip = mesh->tris(t).c[j];
            gzwrite(gzf, &trip, sizeof(int)); 
        }
    }
    
    // per vertex smoke densities
    if (mesh->getType() == Mesh::TypeVortexSheet) {
        VortexSheetMesh* vmesh = (VortexSheetMesh*) mesh;
        int densId[4] = {0, 'v','d','e'};
        gzwrite(gzf, &densId[0], sizeof(int) * 4); 

        // compute densities
        vector<float> triDensity(numTris);
        for (int tri=0; tri < numTris; tri++) {
            Real area = vmesh->getFaceArea(tri);
            if (area>0)
                triDensity[tri] = vmesh->sheet(tri).smokeAmount;
        }
        
        // project triangle data to vertex
        vector<int> triPerVertex(numVerts);
        vector<float> density(numVerts);
        for (int tri=0; tri < numTris; tri++) {
            for (int c=0; c<3; c++) {
                int vertex = mesh->tris(tri).c[c];
                density[vertex] += triDensity[tri];
                triPerVertex[vertex]++;
            }
        }
        
        // averaged smoke densities
        for(int point=0; point<numVerts; point++) {
            float dens = 0;
            if (triPerVertex[point]>0)
                dens = density[point] / triPerVertex[point];
            gzwrite(gzf, &dens, sizeof(float));             
        }
    }
    
    // vertex flags
    if (mesh->getType() == Mesh::TypeVortexSheet) {
        int Id[4] = {0, 'v','x','f'};
        gzwrite(gzf, &Id[0], sizeof(int) * 4); 

        // averaged smoke densities
        for(int point=0; point<numVerts; point++) {
            float alpha = (mesh->nodes(point).flags & Mesh::NfMarked) ? 1: 0;
            gzwrite(gzf, &alpha, sizeof(float));             
        }
    }

    gzclose( gzf );    
}

void readObjFile(const std::string& name, Mesh* mesh, bool append) {
    ifstream ifs (name.c_str());
    
    if (!ifs.good())
        throw Error("can't open file '" + name + "'");
    
    if (!append)
        mesh->clear();
    int nodebase = mesh->numNodes();
    while(ifs.good() && !ifs.eof()) {
        string id;
        ifs >> id;
        
        if (id[0] == '#') {
            // comment
            getline(ifs, id);
            continue;
        }
        if (id == "vt") {
            // tex coord, ignore            
        } else if (id == "vn") {
            // normals, ignore            
        } else if (id == "v") {
            // vertex
            Node n;
            ifs >> n.pos.x >> n.pos.y >> n.pos.z;
            mesh->addNode(n);
        } else if (id == "g") {
            // group
            string group;
            ifs >> group;
        } else if (id == "f") {
            // face
            string face;
            Triangle t;
            for (int i=0; i<3; i++) {
                ifs >> face;
                if (face.find('/') != string::npos)
                    face = face.substr(0, face.find('/')); // ignore other indices
                int idx = atoi(face.c_str()) - 1;
                if (idx < 0)
                    throw Error("invalid face encountered");
                idx += nodebase;
                t.c[i] = idx;
            }
            mesh->addTri(t);
        } else {
            // whatever, ignore
        }
        // kill rest of line
        getline(ifs, id);   
    }
    ifs.close();    
}

template<int DIM, class T>
void writeGridRaw(const string& name, Grid<DIM,T>* grid) {
    cout << "writing grid " << grid->getName() << " to raw file " << name << endl;
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf) throw Error("can't open file");
    gzwrite(gzf, &((*grid)[0]), sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ());
    gzclose(gzf);
}

typedef struct {
    char id[4];
    int dimX, dimY, dimZ;
    int gridType, elementType, bytesPerElement;
} UniHeader;

template <int DIM, class T>
void writeGridUni(const string& name, Grid<DIM,T>* grid) {
    cout << "writing grid " << grid->getName() << " to uni file " << name << endl;
    
    UniHeader head;
	head.id[0]='M';
    head.id[1]='N';
    head.id[2]='T';
    head.id[3]='1';
    head.dimX = grid->getSizeX();
    head.dimY = grid->getSizeY();
    head.dimZ = grid->getSizeZ();
    head.gridType = grid->getType();
    head.bytesPerElement = sizeof(T);
    
    if (grid->getType() & GridBase::TypeInt)
        head.elementType = 0;
    else if (grid->getType() & GridBase::TypeReal)
        head.elementType = 1;
    else if (grid->getType() & GridBase::TypeVec3)
        head.elementType = 2;
    else 
        throw Error("unknown element type");
    
    gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
    if (!gzf) throw Error("can't open file");
    
    gzwrite(gzf, &head, sizeof(UniHeader));
    gzwrite(gzf, &((*grid)[0]), sizeof(T)*head.dimX*head.dimY*head.dimZ);
    gzclose(gzf);
};

template void writeGridRaw<2,int>(const string& name, Grid<2,int>* grid);
template void writeGridRaw<2,Real>(const string& name, Grid<2,Real>* grid);
template void writeGridRaw<2,Vec3>(const string& name, Grid<2,Vec3>* grid);
template void writeGridUni<2,int>(const string& name, Grid<2,int>* grid);
template void writeGridUni<2,Real>(const string& name, Grid<2,Real>* grid);
template void writeGridUni<2,Vec3>(const string& name, Grid<2,Vec3>* grid);
template void writeGridRaw<3,int>(const string& name, Grid<3,int>* grid);
template void writeGridRaw<3,Real>(const string& name, Grid<3,Real>* grid);
template void writeGridRaw<3,Vec3>(const string& name, Grid<3,Vec3>* grid);
template void writeGridUni<3,int>(const string& name, Grid<3,int>* grid);
template void writeGridUni<3,Real>(const string& name, Grid<3,Real>* grid);
template void writeGridUni<3,Vec3>(const string& name, Grid<3,Vec3>* grid);


} //namespace
