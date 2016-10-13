/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2016 Tobias Pfaff, Nils Thuerey  
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
#if NO_ZLIB!=1
extern "C" { 
#include <zlib.h>
}
#endif

#include "fileio.h"
#include "grid.h"
#include "mesh.h"
#include "vortexsheet.h"
#include "particle.h"
#include "vector4d.h"
#include "grid4d.h"
#include  <cstring>

using namespace std;

namespace Manta {

static const int STR_LEN_GRID  = 252;
static const int STR_LEN_PDATA = 256;

//! uni file header, v4
typedef struct {
	int dimX, dimY, dimZ; // grid size
	int gridType, elementType, bytesPerElement; // data type info
	char info[252]; // mantaflow build information
	int dimT;       // optionally store forth dimension for 4d grids
	unsigned long long timestamp; // creation time
} UniHeader;

//! pdata uni header, v3  (similar to grid header)
typedef struct {
	int dim; // number of partilces
	int dimX, dimY, dimZ; // underlying solver resolution (all data in local coordinates!)
	int elementType, bytesPerElement; // type id and byte size
	char info[STR_LEN_PDATA]; // mantaflow build information
	unsigned long long timestamp; // creation time
} UniPartHeader;

//*****************************************************************************
// mesh data
//*****************************************************************************

void readBobjFile(const string& name, Mesh* mesh, bool append) {
	debMsg( "reading mesh file " << name ,1);
	if (!append)
		mesh->clear();
	else
		errMsg("readBobj: append not yet implemented!");

#	if NO_ZLIB!=1
	const Real dx = mesh->getParent()->getDx();
	const Vec3 gs = toVec3( mesh->getParent()->getGridSize() );

	gzFile gzf = gzopen(name.c_str(), "rb1"); // do some compression
	if (!gzf)
		errMsg("readBobj: unable to open file");
	
	// read vertices
	int num = 0;
	gzread(gzf, &num, sizeof(int));
	mesh->resizeNodes(num);
	debMsg( "read mesh , verts "<<num,1);
	for (int i=0; i<num; i++) {
		Vector3D<float> pos;
		gzread(gzf, &pos.value[0], sizeof(float)*3);
	   	mesh->nodes(i).pos = toVec3(pos);

		// convert to grid space
		mesh->nodes(i).pos /= dx;
		mesh->nodes(i).pos += gs*0.5;
	}
	
	// normals
	num = 0;
	gzread(gzf, &num, sizeof(int));
	for (int i=0; i<num; i++) {
		Vector3D<float> pos;
		gzread(gzf, &pos.value[0], sizeof(float)*3);
	   	mesh->nodes(i).normal = toVec3(pos);
	}
	
	// read tris
	num = 0;
	gzread(gzf, &num, sizeof(int));
	mesh->resizeTris( num );
	for(int t=0; t<num; t++) {
		for(int j=0; j<3; j++) { 
			int trip = 0;
			gzread(gzf, &trip, sizeof(int)); 
			mesh->tris(t).c[j] = trip;
		}
	} 
	// note - vortex sheet info ignored for now... (see writeBobj)
	gzclose( gzf );    
	debMsg( "read mesh , triangles "<<mesh->numTris()<<", vertices "<<mesh->numNodes()<<" ",1 );
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

void writeBobjFile(const string& name, Mesh* mesh) {
	debMsg( "writing mesh file " << name ,1);
#	if NO_ZLIB!=1
	const Real  dx = mesh->getParent()->getDx();
	const Vec3i gs = mesh->getParent()->getGridSize();
	
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf)
		errMsg("writeBobj: unable to open file");
	
	// write vertices
	int numVerts = mesh->numNodes();
	gzwrite(gzf, &numVerts, sizeof(int));
	for (int i=0; i<numVerts; i++) {
		Vector3D<float> pos = toVec3f(mesh->nodes(i).pos);
		// normalize to unit cube around 0
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
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

void readObjFile(const std::string& name, Mesh* mesh, bool append) {
	ifstream ifs (name.c_str());
	
	if (!ifs.good())
		errMsg("can't open file '" + name + "'");
	
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
					errMsg("invalid face encountered");
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

// write regular .obj file, in line with bobj.gz output (but only verts & tris for now)
void writeObjFile(const string& name, Mesh* mesh) {
	const Real  dx = mesh->getParent()->getDx();
	const Vec3i gs = mesh->getParent()->getGridSize();

	ofstream ofs(name.c_str());
	if (!ofs.good())
		errMsg("writeObjFile: can't open file " << name);

	ofs << "o MantaMesh\n";
	
	// write vertices
	int numVerts = mesh->numNodes();
	for (int i=0; i<numVerts; i++) {
		Vector3D<float> pos = toVec3f(mesh->nodes(i).pos);
		// normalize to unit cube around 0
		pos -= toVec3f(gs)*0.5;
		pos *= dx;
		ofs << "v "<< pos.value[0] <<" "<< pos.value[1] <<" "<< pos.value[2] <<" "<< "\n";
	}
	
	// no normals for now
	
	// write tris
	int numTris = mesh->numTris();
	for(int t=0; t<numTris; t++) {
		ofs << "f "<< (mesh->tris(t).c[0]+1) <<" "<< (mesh->tris(t).c[1]+1) <<" "<< (mesh->tris(t).c[2]+1) <<" "<< "\n";
	}

	ofs.close();
}

//*****************************************************************************
// conversion functions for double precision
// (note - uni files always store single prec. values)
//*****************************************************************************

#if NO_ZLIB!=1
template <class GRIDT> 
void gridConvertWrite(gzFile& gzf, GRIDT& grid, void* ptr, UniHeader& head) {
	errMsg("gridConvertWrite: unknown type, not yet supported");
}

template <>
void gridConvertWrite(gzFile& gzf, Grid<int>& grid, void* ptr, UniHeader& head) {
	gzwrite(gzf, &head,    sizeof(UniHeader));
	gzwrite(gzf, &grid[0], sizeof(int)*head.dimX*head.dimY*head.dimZ);
} 
template <>
void gridConvertWrite(gzFile& gzf, Grid<double>& grid, void* ptr, UniHeader& head) {
	head.bytesPerElement = sizeof(float);
	gzwrite(gzf, &head, sizeof(UniHeader));
	float* ptrf = (float*)ptr;
	for(int i=0; i<grid.getSizeX()*grid.getSizeY()*grid.getSizeZ(); ++i,++ptrf) {
		*ptrf = (float)grid[i];
	} 
	gzwrite(gzf, ptr, sizeof(float)* head.dimX*head.dimY*head.dimZ);
} 
template <>
void gridConvertWrite(gzFile& gzf, Grid<Vector3D<double> >& grid, void* ptr, UniHeader& head) {
	head.bytesPerElement = sizeof(Vector3D<float>);
	gzwrite(gzf, &head, sizeof(UniHeader));
	float* ptrf = (float*)ptr;
	for(int i=0; i<grid.getSizeX()*grid.getSizeY()*grid.getSizeZ(); ++i) {
		for(int c=0; c<3; ++c) { *ptrf = (float)grid[i][c]; ptrf++; }
	} 
	gzwrite(gzf, ptr, sizeof(Vector3D<float>) *head.dimX*head.dimY*head.dimZ);
}

template <>
void gridConvertWrite(gzFile& gzf, Grid4d<int>& grid, void* ptr, UniHeader& head) {
	gzwrite(gzf, &head,    sizeof(UniHeader));
	gzwrite(gzf, &grid[0], sizeof(int)*head.dimX*head.dimY*head.dimZ*head.dimT);
}
template <>
void gridConvertWrite(gzFile& gzf, Grid4d<double>& grid, void* ptr, UniHeader& head) {
	head.bytesPerElement = sizeof(float);
	gzwrite(gzf, &head, sizeof(UniHeader));
	float* ptrf = (float*)ptr;
	IndexInt s = grid.getStrideT()*grid.getSizeT();
	for(IndexInt i=0; i<s; ++i,++ptrf) {
		*ptrf = (float)grid[i];
	} 
	gzwrite(gzf, ptr, sizeof(float)* s );
} 
template <>
void gridConvertWrite(gzFile& gzf, Grid4d<Vector3D<double> >& grid, void* ptr, UniHeader& head) {
	head.bytesPerElement = sizeof(Vector3D<float>);
	gzwrite(gzf, &head, sizeof(UniHeader));
	float* ptrf = (float*)ptr;
	IndexInt s = grid.getStrideT()*grid.getSizeT();
	for(IndexInt i=0; i<s; ++i) {
		for(int c=0; c<3; ++c) { *ptrf = (float)grid[i][c]; ptrf++; }
	} 
	gzwrite(gzf, ptr, sizeof(Vector3D<float>) *s);
}
template <>
void gridConvertWrite(gzFile& gzf, Grid4d<Vector4D<double> >& grid, void* ptr, UniHeader& head) {
	head.bytesPerElement = sizeof(Vector4D<float>);
	gzwrite(gzf, &head, sizeof(UniHeader));
	float* ptrf = (float*)ptr;
	IndexInt s = grid.getStrideT()*grid.getSizeT();
	for(IndexInt i=0; i<s; ++i) {
		for(int c=0; c<4; ++c) { *ptrf = (float)grid[i][c]; ptrf++; }
	} 
	gzwrite(gzf, ptr, sizeof(Vector4D<float>) *s);
}

template <class T>
void pdataConvertWrite( gzFile& gzf, ParticleDataImpl<T>& pdata, void* ptr, UniPartHeader& head) {
	errMsg("pdataConvertWrite: unknown type, not yet supported");
}

template <>
void pdataConvertWrite( gzFile& gzf, ParticleDataImpl<int>& pdata, void* ptr, UniPartHeader& head) {
	gzwrite(gzf, &head,     sizeof(UniPartHeader));
	gzwrite(gzf, &pdata[0], sizeof(int)*head.dim);
} 
template <>
void pdataConvertWrite( gzFile& gzf, ParticleDataImpl<double>& pdata, void* ptr, UniPartHeader& head) {
	head.bytesPerElement = sizeof(float);
	gzwrite(gzf, &head, sizeof(UniPartHeader));
	float* ptrf = (float*)ptr;
	for(int i=0; i<pdata.size(); ++i,++ptrf) {
		*ptrf = (float)pdata[i];
	} 
	gzwrite(gzf, ptr, sizeof(float)* head.dim);
} 
template <>
void pdataConvertWrite( gzFile& gzf, ParticleDataImpl<Vec3>& pdata, void* ptr, UniPartHeader& head) {
	head.bytesPerElement = sizeof(Vector3D<float>);
	gzwrite(gzf, &head, sizeof(UniPartHeader));
	float* ptrf = (float*)ptr;
	for(int i=0; i<pdata.size(); ++i) {
		for(int c=0; c<3; ++c) { *ptrf = (float)pdata[i][c]; ptrf++; }
	} 
	gzwrite(gzf, ptr, sizeof(Vector3D<float>) *head.dim);
}

template <class T>
void gridReadConvert(gzFile& gzf, Grid<T>& grid, void* ptr, int bytesPerElement) {
	errMsg("gridReadConvert: unknown type, not yet supported");
}

template <>
void gridReadConvert<int>(gzFile& gzf, Grid<int>& grid, void* ptr, int bytesPerElement) {
	gzread(gzf, ptr, sizeof(int)*grid.getSizeX()*grid.getSizeY()*grid.getSizeZ());
	assertMsg (bytesPerElement == sizeof(int), "grid element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(int) );
	// easy, nothing to do for ints
	memcpy(&(grid[0]), ptr, sizeof(int) * grid.getSizeX()*grid.getSizeY()*grid.getSizeZ() );
}

template <>
void gridReadConvert<double>(gzFile& gzf, Grid<double>& grid, void* ptr, int bytesPerElement) {
	gzread(gzf, ptr, sizeof(float)*grid.getSizeX()*grid.getSizeY()*grid.getSizeZ());
	assertMsg (bytesPerElement == sizeof(float), "grid element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(float) );
	float* ptrf = (float*)ptr;
	for(int i=0; i<grid.getSizeX()*grid.getSizeY()*grid.getSizeZ(); ++i,++ptrf) {
		grid[i] = (double)(*ptrf);
	} 
}

template <>
void gridReadConvert<Vec3>(gzFile& gzf, Grid<Vec3>& grid, void* ptr, int bytesPerElement) {
	gzread(gzf, ptr, sizeof(Vector3D<float>)*grid.getSizeX()*grid.getSizeY()*grid.getSizeZ());
	assertMsg (bytesPerElement == sizeof(Vector3D<float>), "grid element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(Vector3D<float>) );
	float* ptrf = (float*)ptr;
	for(int i=0; i<grid.getSizeX()*grid.getSizeY()*grid.getSizeZ(); ++i) {
		Vec3 v;
		for(int c=0; c<3; ++c) { v[c] = double(*ptrf); ptrf++; }
		grid[i] = v;
	} 
}

template <class T>
void pdataReadConvert(gzFile& gzf, ParticleDataImpl<T>& grid, void* ptr, int bytesPerElement) {
	errMsg("pdataReadConvert: unknown pdata type, not yet supported");
}

template <>
void pdataReadConvert<int>(gzFile& gzf, ParticleDataImpl<int>& pdata, void* ptr, int bytesPerElement) {
	gzread(gzf, ptr, sizeof(int)*pdata.size());
	assertMsg (bytesPerElement == sizeof(int), "pdata element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(int) );
	// int dont change in double precision mode - copy over
	memcpy(&(pdata[0]), ptr, sizeof(int) * pdata.size() );
}

template <>
void pdataReadConvert<double>(gzFile& gzf, ParticleDataImpl<double>& pdata, void* ptr, int bytesPerElement) {
	gzread(gzf, ptr, sizeof(float)*pdata.size());
	assertMsg (bytesPerElement == sizeof(float), "pdata element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(float) );
	float* ptrf = (float*)ptr;
	for(int i=0; i<pdata.size(); ++i,++ptrf) {
		pdata[i] = double(*ptrf); 
	} 
}

template <>
void pdataReadConvert<Vec3>(gzFile& gzf, ParticleDataImpl<Vec3>& pdata, void* ptr, int bytesPerElement) {
	gzread(gzf, ptr, sizeof(Vector3D<float>)*pdata.size());
	assertMsg (bytesPerElement == sizeof(Vector3D<float>), "pdata element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(Vector3D<float>) );
	float* ptrf = (float*)ptr;
	for(int i=0; i<pdata.size(); ++i) {
		Vec3 v;
		for(int c=0; c<3; ++c) { v[c] = double(*ptrf); ptrf++; }
		pdata[i] = v;
	} 
}

template <class T>
void gridReadConvert4d(gzFile& gzf, Grid4d<T>& grid, void* ptr, int bytesPerElement, int t) {
	errMsg("gridReadConvert4d: unknown type, not yet supported");
}

template <>
void gridReadConvert4d<int>(gzFile& gzf, Grid4d<int>& grid, void* ptr, int bytesPerElement, int t) { 
	gzread(gzf, ptr, sizeof(int)*grid.getSizeX()*grid.getSizeY()*grid.getSizeZ() );
	assertMsg (bytesPerElement == sizeof(int), "grid element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(int) );
	// nothing to do for ints
	memcpy(&(grid[grid.getSizeX()*grid.getSizeY()*grid.getSizeZ() *t]), ptr, sizeof(int) * grid.getSizeX()*grid.getSizeY()*grid.getSizeZ() );
}

template <>
void gridReadConvert4d<double>(gzFile& gzf, Grid4d<double>& grid, void* ptr, int bytesPerElement, int t) {
	assertMsg (bytesPerElement == sizeof(float), "grid element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(float) );

	float* ptrf = (float*)ptr;
	gzread(gzf, ptr, sizeof(float)*grid.getSizeX()*grid.getSizeY()*grid.getSizeZ());
	for(IndexInt i=0; i<grid.getSizeX()*grid.getSizeY()*grid.getSizeZ(); ++i , ++ptrf ) {
		grid[ grid.getSizeX()*grid.getSizeY()*grid.getSizeZ()*t + i] = (double)(*ptrf);
	} 
}

template <>
void gridReadConvert4d<Vec3>(gzFile& gzf, Grid4d<Vec3>& grid, void* ptr, int bytesPerElement, int t) {
	assertMsg (bytesPerElement == sizeof(Vector3D<float>), "grid element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(float) );

	gzread(gzf, ptr, sizeof(Vector3D<float>)*grid.getSizeX()*grid.getSizeY()*grid.getSizeZ());
	float* ptrf = (float*)ptr;
	for(IndexInt i=0; i<grid.getSizeX()*grid.getSizeY()*grid.getSizeZ(); ++i) {
		Vec3 v;
		for(int c=0; c<3; ++c) { v[c] = double(*ptrf); ptrf++; }
		grid[grid.getSizeX()*grid.getSizeY()*grid.getSizeZ()*t + i] = v;
	} 
}

template <>
void gridReadConvert4d<Vec4>(gzFile& gzf, Grid4d<Vec4>& grid, void* ptr, int bytesPerElement, int t) {
	assertMsg (bytesPerElement == sizeof(Vector4D<float>), "grid element size doesn't match "<< bytesPerElement <<" vs "<< sizeof(float) );

	gzread(gzf, ptr, sizeof(Vector4D<float>)*grid.getSizeX()*grid.getSizeY()*grid.getSizeZ());
	float* ptrf = (float*)ptr;
	for(IndexInt i=0; i<grid.getSizeX()*grid.getSizeY()*grid.getSizeZ(); ++i) {
		Vec4 v;
		for(int c=0; c<4; ++c) { v[c] = double(*ptrf); ptrf++; }
		grid[grid.getSizeX()*grid.getSizeY()*grid.getSizeZ()*t + i] = v;
	} 
}



// make sure compatible grid types dont lead to errors...
static int unifyGridType(int type) {
	// real <> levelset
	if(type & GridBase::TypeReal)     type |= GridBase::TypeLevelset;
	if(type & GridBase::TypeLevelset) type |= GridBase::TypeReal;
	// vec3 <> mac
	if(type & GridBase::TypeVec3)     type |= GridBase::TypeMAC;
	if(type & GridBase::TypeMAC)      type |= GridBase::TypeVec3;
	return type;
}

#endif // NO_ZLIB!=1

//*****************************************************************************
// grid data
//*****************************************************************************

template<class T>
void writeGridTxt(const string& name, Grid<T>* grid) {
	debMsg( "writing grid " << grid->getName() << " to text file " << name ,1);

	ofstream ofs(name.c_str());
	if (!ofs.good())
		errMsg("can't open file " << name);
	FOR_IJK(*grid) {
		ofs << Vec3i(i,j,k) <<" = "<< (*grid)(i,j,k) <<"\n";
	}
	ofs.close();
}

template<class T>
void writeGridRaw(const string& name, Grid<T>* grid) {
	debMsg( "writing grid " << grid->getName() << " to raw file " << name ,1);
	
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf) errMsg("can't open file " << name);
	gzwrite(gzf, &((*grid)[0]), sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ());
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

template<class T>
void readGridRaw(const string& name, Grid<T>* grid) {
	debMsg( "reading grid " << grid->getName() << " from raw file " << name ,1);
	
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (!gzf) errMsg("can't open file " << name);
	
	IndexInt bytes = sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ();
	IndexInt readBytes = gzread(gzf, &((*grid)[0]), bytes);
	assertMsg(bytes==readBytes, "can't read raw file, stream length does not match, "<<bytes<<" vs "<<readBytes);
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

//! legacy headers for reading old files
typedef struct {
	int dimX, dimY, dimZ;
	int frames, elements, elementType, bytesPerElement, bytesPerFrame;
} UniLegacyHeader;

typedef struct {
	int dimX, dimY, dimZ;
	int gridType, elementType, bytesPerElement;
} UniLegacyHeader2;

typedef struct {
	int dimX, dimY, dimZ;
	int gridType, elementType, bytesPerElement;
	char info[256];
	unsigned long long timestamp;
} UniLegacyHeader3;

//! for auto-init & check of results of test runs , optionally returns info string of header
void getUniFileSize(const string& name, int& x, int& y, int& z, int* t, std::string* info) {
	x = y = z = 0;
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (gzf) { 
		char ID[5]={0,0,0,0,0};
		gzread(gzf, ID, 4); 

		// v3
		if ( (!strcmp(ID, "MNT2")) || (!strcmp(ID, "M4T2")) ) {
			UniLegacyHeader3 head;
			assertMsg (gzread(gzf, &head, sizeof(UniLegacyHeader3)) == sizeof(UniLegacyHeader3), "can't read file, no header present"); 
			x = head.dimX;
			y = head.dimY;
			z = head.dimZ;

			// optionally , read fourth dim
			if ((!strcmp(ID, "M4T2")) && t) {
				int dimT = 0;
				gzread(gzf, &dimT, sizeof(int) );
				(*t) = dimT;
			}
		}

		// v4
		if ( (!strcmp(ID, "MNT3")) || (!strcmp(ID, "M4T3")) ) {
			UniHeader head;
			assertMsg (gzread(gzf, &head, sizeof(UniHeader)) == sizeof(UniHeader), "can't read file, no header present"); 
			x = head.dimX;
			y = head.dimY;
			z = head.dimZ;
			if(t) (*t) = head.dimT;
		}

		gzclose(gzf);
	}
#	endif
}
PYTHON() Vec3 getUniFileSize(const string& name) {
	int x,y,z;
	getUniFileSize(name, x,y,z);
	return Vec3( Real(x), Real(y), Real(z) );
}

//! for test run debugging
PYTHON() void printUniFileInfoString(const string& name) {
	std::string info("<file not found>");
	int x=-1,y=-1,z=-1,t=-1;
	// use getUniFileSize to parse the different headers
	getUniFileSize(name, x,y,z,&t, &info);
	debMsg("File '"<<name<<"' info: "<< info ,1);
}


// actual read/write functions

template <class T>
void writeGridUni(const string& name, Grid<T>* grid) {
	debMsg( "Writing grid " << grid->getName() << " to uni file " << name ,1);
	
#	if NO_ZLIB!=1
	char ID[5] = "MNT3";
	UniHeader head;
	head.dimX = grid->getSizeX();
	head.dimY = grid->getSizeY();
	head.dimZ = grid->getSizeZ();
	head.dimT = 0;
	head.gridType = grid->getType();
	head.bytesPerElement = sizeof(T);
	snprintf( head.info, STR_LEN_GRID, "%s", buildInfoString().c_str() );	
	MuTime stamp;
	head.timestamp = stamp.time;
	
	if (grid->getType() & GridBase::TypeInt)
		head.elementType = 0;
	else if (grid->getType() & GridBase::TypeReal)
		head.elementType = 1;
	else if (grid->getType() & GridBase::TypeVec3)
		head.elementType = 2;
	else 
		errMsg("unknown element type");
	
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf) errMsg("can't open file " << name);
	
	gzwrite(gzf, ID, 4);
#	if FLOATINGPOINT_PRECISION!=1
	// always write float values, even if compiled with double precision...
	Grid<T> temp(grid->getParent());
	// "misuse" temp grid as storage for floating point values (we have double, so it will always fit)
	gridConvertWrite( gzf, *grid, &(temp[0]), head);
#	else
	void* ptr = &((*grid)[0]);
	gzwrite(gzf, &head, sizeof(UniHeader));
	gzwrite(gzf, ptr, sizeof(T)*head.dimX*head.dimY*head.dimZ);
#	endif
	gzclose(gzf);

#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};


template <class T>
void readGridUni(const string& name, Grid<T>* grid) {
	debMsg( "Reading grid " << grid->getName() << " from uni file " << name ,1);

#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (!gzf) errMsg("can't open file " << name);

	char ID[5]={0,0,0,0,0};
	gzread(gzf, ID, 4);
	
	if (!strcmp(ID, "DDF2")) {
		// legacy file format
		UniLegacyHeader head;
		assertMsg (gzread(gzf, &head, sizeof(UniLegacyHeader)) == sizeof(UniLegacyHeader), "can't read file, no header present");
		assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match");
		assertMsg (head.bytesPerElement * head.elements == sizeof(T), "grid type doesn't match");
		// skip flags
		int numEl = head.dimX*head.dimY*head.dimZ;
		gzseek(gzf, numEl, SEEK_CUR);
		// actual grid read
		gzread(gzf, &((*grid)[0]), sizeof(T)*numEl);
	} 
	else if (!strcmp(ID, "MNT1")) {
		// legacy file format 2
		UniLegacyHeader2 head;
		assertMsg (gzread(gzf, &head, sizeof(UniLegacyHeader2)) == sizeof(UniLegacyHeader2), "can't read file, no header present");
		assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match, "<< Vec3(head.dimX,head.dimY,head.dimZ)<<" vs "<< grid->getSize() );
		assertMsg (head.gridType == grid->getType(), "grid type doesn't match "<< head.gridType<<" vs "<< grid->getType() );
		assertMsg (head.bytesPerElement == sizeof(T), "grid element size doesn't match "<< head.bytesPerElement <<" vs "<< sizeof(T) );
		gzread(gzf, &((*grid)[0]), sizeof(T)*head.dimX*head.dimY*head.dimZ);
	}
	else if (!strcmp(ID, "MNT2")) {
		// a bit ugly, almost identical to MNT3
		UniLegacyHeader3 head;
		assertMsg (gzread(gzf, &head, sizeof(UniLegacyHeader3)) == sizeof(UniLegacyHeader3), "can't read file, no header present");
		assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match, "<< Vec3(head.dimX,head.dimY,head.dimZ)<<" vs "<< grid->getSize() );
		assertMsg (unifyGridType(head.gridType)==unifyGridType(grid->getType()) , "grid type doesn't match "<< head.gridType<<" vs "<< grid->getType() );
#		if FLOATINGPOINT_PRECISION!=1
		Grid<T> temp(grid->getParent());
		void*  ptr  = &(temp[0]);
		gridReadConvert<T>(gzf, *grid, ptr, head.bytesPerElement);
#		else
		assertMsg (head.bytesPerElement == sizeof(T), "grid element size doesn't match "<< head.bytesPerElement <<" vs "<< sizeof(T) );
		gzread(gzf, &((*grid)[0]), sizeof(T)*head.dimX*head.dimY*head.dimZ);
#		endif
	} 
	else if (!strcmp(ID, "MNT3")) {
		// current file format
		UniHeader head;
		assertMsg (gzread(gzf, &head, sizeof(UniHeader)) == sizeof(UniHeader), "can't read file, no header present");
		assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match, "<< Vec3(head.dimX,head.dimY,head.dimZ)<<" vs "<< grid->getSize() );
		assertMsg (unifyGridType(head.gridType)==unifyGridType(grid->getType()) , "grid type doesn't match "<< head.gridType<<" vs "<< grid->getType() );
#		if FLOATINGPOINT_PRECISION!=1
		// convert float to double
		Grid<T> temp(grid->getParent());
		void*  ptr  = &(temp[0]);
		gridReadConvert<T>(gzf, *grid, ptr, head.bytesPerElement);
#		else
		assertMsg (head.bytesPerElement == sizeof(T), "grid element size doesn't match "<< head.bytesPerElement <<" vs "<< sizeof(T) );
		gzread(gzf, &((*grid)[0]), sizeof(T)*head.dimX*head.dimY*head.dimZ);
#		endif
	} else {
		debMsg( "Unknown header!" ,1);
	}
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};

template <class T>
void writeGridVol(const string& name, Grid<T>* grid) {
	debMsg( "writing grid " << grid->getName() << " to vol file " << name ,1);
	errMsg("Type not yet supported!");
}

struct volHeader { 
	char 	ID[3];
	char 	version;
	int		encoding;
	int		dimX, dimY, dimZ;
	int		channels;
	Vec3	bboxMin, bboxMax;
};

template <>
void writeGridVol<Real>(const string& name, Grid<Real>* grid) {
	debMsg( "writing real grid " << grid->getName() << " to vol file " << name ,1);
	
	volHeader header;
	header.ID[0] = 'V';
	header.ID[1] = 'O';
	header.ID[2] = 'L';
	header.version  = 3;
	header.encoding = 1; // float32 precision
	header.dimX = grid->getSizeX();
	header.dimY = grid->getSizeY();
	header.dimZ = grid->getSizeZ();
	header.channels = 1; // only 1 channel
	header.bboxMin = Vec3(-0.5);
	header.bboxMax = Vec3( 0.5);

	FILE* fp = fopen( name.c_str(), "wb" );
	if (fp == NULL) {
		errMsg("Cannot open '" << name << "'");
		return;
	}

	fwrite( &header, sizeof(volHeader), 1, fp );

#	if FLOATINGPOINT_PRECISION==1
	// for float, write one big chunk
	fwrite( &(*grid)[0], sizeof(float), grid->getSizeX()*grid->getSizeY()*grid->getSizeZ(), fp );
#	else
	// explicitly convert each entry to float - we might have double precision in mantaflow
	FOR_IDX(*grid) {
		float value = (*grid)[idx];
		fwrite( &value, sizeof(float), 1, fp );
	}
#	endif

	fclose(fp);
};


template <class T>
void readGridVol(const string& name, Grid<T>* grid) {
	debMsg( "writing grid " << grid->getName() << " to vol file " << name ,1);
	errMsg("Type not yet supported!");
}

template <>
void readGridVol<Real>(const string& name, Grid<Real>* grid) {
	debMsg( "reading real grid " << grid->getName() << " from vol file " << name ,1);
	
	volHeader header; 
	FILE* fp = fopen( name.c_str(), "rb" );
	if (fp == NULL) {
		errMsg("Cannot open '" << name << "'");
		return;
	}

	// note, only very basic file format checks here!
	assertMsg( fread( &header, 1, sizeof(volHeader) ,  fp) == sizeof(volHeader), "can't read file, no header present");
	if( header.dimX != grid->getSizeX() || header.dimY != grid->getSizeY() || header.dimZ != grid->getSizeZ()) errMsg( "grid dim doesn't match, "<< Vec3(header.dimX,header.dimY,header.dimZ)<<" vs "<< grid->getSize() );
#	if FLOATINGPOINT_PRECISION!=1
	errMsg("Not yet supported");
#	else
	fread( &((*grid)[0]), 1, sizeof(float)*header.dimX*header.dimY*header.dimZ , fp);
#	endif

	fclose(fp);
};

// 4d grids IO

template <class T>
void writeGrid4dUni(const string& name, Grid4d<T>* grid) {
	debMsg( "writing grid4d " << grid->getName() << " to uni file " << name ,1);
	
#	if NO_ZLIB!=1
	char ID[5] = "M4T3";
	UniHeader head;
	head.dimX = grid->getSizeX();
	head.dimY = grid->getSizeY();
	head.dimZ = grid->getSizeZ();
	head.dimT = grid->getSizeT();
	head.gridType = grid->getType();
	head.bytesPerElement = sizeof(T);
	snprintf( head.info, STR_LEN_GRID, "%s", buildInfoString().c_str() );	
	MuTime stamp; stamp.get();
	head.timestamp = stamp.time;
	
	if (grid->getType() & Grid4dBase::TypeInt)
		head.elementType = 0;
	else if (grid->getType() & Grid4dBase::TypeReal)
		head.elementType = 1;
	else if (grid->getType() & Grid4dBase::TypeVec3)
		head.elementType = 2;
	else if (grid->getType() & Grid4dBase::TypeVec4)
		head.elementType = 2;
	else 
		errMsg("unknown element type");
	
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf) errMsg("can't open file " << name);
	
	gzwrite(gzf, ID, 4);
#	if FLOATINGPOINT_PRECISION!=1
	Grid4d<T> temp(grid->getParent());
	gridConvertWrite< Grid4d<T> >( gzf, *grid, &(temp[0]), head);
#	else
	gzwrite(gzf, &head, sizeof(UniHeader));

	// can be too large - write in chunks
	for(int t=0; t<head.dimT; ++t) { 
		void* ptr = &((*grid)[           head.dimX*head.dimY*head.dimZ* t ]);
		gzwrite(gzf, ptr,      sizeof(T)*head.dimX*head.dimY*head.dimZ* 1);
	}
#	endif
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};

template <class T>
void readGrid4dUni(const string& name, Grid4d<T>* grid, int readTslice, Grid4d<T>* slice, void** fileHandle ) 
{
	if(grid)  debMsg( "reading grid "  << grid->getName()  << " from uni file " << name ,1);
	if(slice) debMsg( "reading slice " << slice->getName() << ",t="<<readTslice<<" from uni file " << name ,1);

#	if NO_ZLIB!=1
	gzFile gzf = NULL;
	char ID[5]={0,0,0,0,0};

	// optionally - reuse file handle, if valid one is passed in fileHandle pointer...
	if( (!fileHandle) || (fileHandle && (*fileHandle == NULL)) ) {
		gzf = gzopen(name.c_str(), "rb");
		if (!gzf) errMsg("can't open file "<<name);

		gzread(gzf, ID, 4);
		if( fileHandle) { *fileHandle = gzf; }
	} else {
		// optimized read - reduced sanity checks
		gzf = (gzFile)(*fileHandle); 
		void* ptr = &( (*slice)[ 0 ] );
		gzread(gzf, ptr, sizeof(T)* slice->getStrideT()* 1);  // quick and dirty...
		return;
	}
	
	if( (!strcmp(ID, "M4T2")) || (!strcmp(ID, "M4T3")) ) {
		int headerSize = -1;

		// current file format
		UniHeader head;
		if(!strcmp(ID, "M4T3")) {
			headerSize = sizeof(UniHeader);
			assertMsg (gzread(gzf, &head, sizeof(UniHeader)) == sizeof(UniHeader), "can't read file, no 4d header present");
			if(FLOATINGPOINT_PRECISION==1) assertMsg (head.bytesPerElement == sizeof(T), "4d grid element size doesn't match "<< head.bytesPerElement <<" vs "<< sizeof(T) );
		}
		// old header
		if(!strcmp(ID, "M4T2")) {
			UniLegacyHeader3 lhead;
			headerSize = sizeof(UniLegacyHeader3) + sizeof(int);
			assertMsg (gzread(gzf, &lhead, sizeof(UniLegacyHeader3)) == sizeof(UniLegacyHeader3), "can't read file, no 4dl header present");
			if(FLOATINGPOINT_PRECISION==1) assertMsg (lhead.bytesPerElement == sizeof(T), "4d grid element size doesn't match "<< lhead.bytesPerElement <<" vs "<< sizeof(T) ); 

			int fourthDim  = 0;
			gzread(gzf, &fourthDim, sizeof(fourthDim));

			head.dimX = lhead.dimX;
			head.dimY = lhead.dimY;
			head.dimZ = lhead.dimZ;
			head.dimT = fourthDim;
			head.gridType = lhead.gridType;
		}

		if(readTslice<0) {
			assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match, "<< Vec3(head.dimX,head.dimY,head.dimZ)<<" vs "<< grid->getSize() );
			assertMsg ( unifyGridType(head.gridType)==unifyGridType(grid->getType()) , "grid type doesn't match "<< head.gridType<<" vs "<< grid->getType() );

			// read full 4d grid
			assertMsg (head.dimT == grid->getSizeT(), "grid dim4 doesn't match, "<< head.dimT <<" vs "<< grid->getSize() );

			// can be too large - read in chunks
#			if FLOATINGPOINT_PRECISION!=1
			Grid4d<T> temp(grid->getParent());
			void* ptr = &(temp[0]); 
			for(int t=0; t<head.dimT; ++t) { 
				gridReadConvert4d<T>(gzf, *grid, ptr, head.bytesPerElement, t);
			}
#			else
			for(int t=0; t<head.dimT; ++t) { 
				void* ptr = &((*grid)[            head.dimX*head.dimY*head.dimZ* t ] );
				gzread(gzf, ptr,       sizeof(T)* head.dimX*head.dimY*head.dimZ*    1); 
			}
#			endif
		} else {
			// read chosen slice only
			assertMsg (head.dimX == slice->getSizeX() && head.dimY == slice->getSizeY() && head.dimZ == slice->getSizeZ(), "grid dim doesn't match, "<< Vec3(head.dimX,head.dimY,head.dimZ)<<" vs "<< slice->getSize() );
			assertMsg ( unifyGridType(head.gridType)==unifyGridType(slice->getType()) , "grid type doesn't match "<< head.gridType<<" vs "<< slice->getType() );

#			if FLOATINGPOINT_PRECISION!=1
			errMsg( "NYI (2)" ); // slice read not yet supported for double
#			else
			assertMsg( slice, "No 3d slice grid data given" );
			assertMsg (readTslice < head.dimT, "grid dim4 slice too large "<< readTslice <<" vs "<< head.dimT );
			void* ptr = &( (*slice)[ 0 ] );
			gzseek(gzf, sizeof(T)*head.dimX*head.dimY*head.dimZ* readTslice + headerSize + 4 , SEEK_SET );
			gzread(gzf, ptr,         sizeof(T)*head.dimX*head.dimY*head.dimZ* 1);
#			endif
		}
	} else {
		debMsg( "Unknown header!" ,1);
	}

	if( !fileHandle) { 
		gzclose(gzf);
	}
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};
void readGrid4dUniCleanup(void** fileHandle) {
	gzFile gzf = NULL; 
	if( fileHandle) {
		gzf = (gzFile)(*fileHandle);
		gzclose(gzf);
		*fileHandle = NULL;
	}
}

template<class T>
void writeGrid4dRaw(const string& name, Grid4d<T>* grid) {
	debMsg( "writing grid4d " << grid->getName() << " to raw file " << name ,1);
	
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf) errMsg("can't open file " << name);
	gzwrite(gzf, &((*grid)[0]), sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ()*grid->getSizeT());
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

template<class T>
void readGrid4dRaw(const string& name, Grid4d<T>* grid) {
	debMsg( "reading grid4d " << grid->getName() << " from raw file " << name ,1);
	
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (!gzf) errMsg("can't open file " << name);
	
	IndexInt bytes = sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ()*grid->getSizeT();
	IndexInt readBytes = gzread(gzf, &((*grid)[0]), bytes);
	assertMsg(bytes==readBytes, "can't read raw file, stream length does not match, "<<bytes<<" vs "<<readBytes);
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}


//*****************************************************************************
// particle data
//*****************************************************************************

static const int PartSysSize = sizeof(Vector3D<float>)+sizeof(int);

void writeParticlesUni(const std::string& name, BasicParticleSystem* parts ) {
	debMsg( "writing particles " << parts->getName() << " to uni file " << name ,1);
	
#	if NO_ZLIB!=1
	char ID[5] = "PB02";
	UniPartHeader head;
	head.dim      = parts->size();
	Vec3i         gridSize = parts->getParent()->getGridSize();
	head.dimX     = gridSize.x;
	head.dimY     = gridSize.y;
	head.dimZ     = gridSize.z;
	head.bytesPerElement = PartSysSize;
	head.elementType = 0; // 0 for base data
	snprintf( head.info, STR_LEN_PDATA, "%s", buildInfoString().c_str() );	
	MuTime stamp;
	head.timestamp = stamp.time;
	
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf) errMsg("can't open file " << name);
	
	gzwrite(gzf, ID, 4);
#	if FLOATINGPOINT_PRECISION!=1
	// warning - hard coded conversion of byte size here...
	gzwrite(gzf, &head, sizeof(UniPartHeader));
	for(int i=0; i<parts->size(); ++i) {
		Vector3D<float> pos  = toVec3f( (*parts)[i].pos );
		int             flag = (*parts)[i].flag;
		gzwrite(gzf, &pos , sizeof(Vector3D<float>) );
		gzwrite(gzf, &flag, sizeof(int)             );
	}
#	else
	assertMsg( sizeof(BasicParticleData) == PartSysSize, "particle data size doesn't match" );
	gzwrite(gzf, &head, sizeof(UniPartHeader));
	gzwrite(gzf, &(parts->getData()[0]), PartSysSize*head.dim);
#	endif
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};

void readParticlesUni(const std::string& name, BasicParticleSystem* parts ) {
	debMsg( "reading particles " << parts->getName() << " from uni file " << name ,1);
	
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (!gzf) errMsg("can't open file " << name);

	char ID[5]={0,0,0,0,0};
	gzread(gzf, ID, 4);
	
	if (!strcmp(ID, "PB01")) {
		errMsg("particle uni file format v01 not supported anymore");
	} else if (!strcmp(ID, "PB02")) {
		// current file format
		UniPartHeader head;
		assertMsg (gzread(gzf, &head, sizeof(UniPartHeader)) == sizeof(UniPartHeader), "can't read file, no header present");
		assertMsg ( ((head.bytesPerElement == PartSysSize) && (head.elementType==0) ), "particle type doesn't match");

		// re-allocate all data
		parts->resizeAll( head.dim );

		assertMsg (head.dim == parts->size() , "particle size doesn't match");
#		if FLOATINGPOINT_PRECISION!=1
		for(int i=0; i<parts->size(); ++i) {
			Vector3D<float> pos; int flag;
			gzread(gzf, &pos , sizeof(Vector3D<float>) );
			gzread(gzf, &flag, sizeof(int)             );
			(*parts)[i].pos  = toVec3d(pos);
			(*parts)[i].flag = flag;
		}
#		else
		assertMsg( sizeof(BasicParticleData) == PartSysSize, "particle data size doesn't match" );
		IndexInt bytes     = PartSysSize*head.dim;
		IndexInt readBytes = gzread(gzf, &(parts->getData()[0]), bytes);
		assertMsg( bytes==readBytes, "can't read uni file, stream length does not match, "<<bytes<<" vs "<<readBytes );
#		endif

		parts->transformPositions( Vec3i(head.dimX,head.dimY,head.dimZ), parts->getParent()->getGridSize() );
	}
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};

template <class T>
void writePdataUni(const std::string& name, ParticleDataImpl<T>* pdata ) {
	debMsg( "writing particle data " << pdata->getName() << " to uni file " << name ,1);
	
#	if NO_ZLIB!=1
	char ID[5] = "PD01";
	UniPartHeader head;
	head.dim      = pdata->size();
	head.bytesPerElement = sizeof(T);
	head.elementType = 1; // 1 for particle data, todo - add sub types?
	snprintf( head.info, STR_LEN_PDATA, "%s", buildInfoString().c_str() );	
	MuTime stamp;
	head.timestamp = stamp.time;
	
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf) errMsg("can't open file " << name);
	gzwrite(gzf, ID, 4);

#	if FLOATINGPOINT_PRECISION!=1
	// always write float values, even if compiled with double precision (as for grids)
	ParticleDataImpl<T> temp(pdata->getParent());
	temp.resize( pdata->size() );
	pdataConvertWrite( gzf, *pdata, &(temp[0]), head);
#	else
	gzwrite(gzf, &head, sizeof(UniPartHeader));
	gzwrite(gzf, &(pdata->get(0)), sizeof(T)*head.dim);
#	endif
	gzclose(gzf);

#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};

template <class T>
void readPdataUni(const std::string& name, ParticleDataImpl<T>* pdata ) {
	debMsg( "reading particle data " << pdata->getName() << " from uni file " << name ,1);
	
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (!gzf) errMsg("can't open file " << name );

	char ID[5]={0,0,0,0,0};
	gzread(gzf, ID, 4);
	
	if (!strcmp(ID, "PD01")) {
		UniPartHeader head;
		assertMsg (gzread(gzf, &head, sizeof(UniPartHeader)) == sizeof(UniPartHeader), "can't read file, no header present");
		assertMsg (head.dim == pdata->size() , "pdata size doesn't match");
#		if FLOATINGPOINT_PRECISION!=1
		ParticleDataImpl<T> temp(pdata->getParent());
		temp.resize( pdata->size() );
		pdataReadConvert<T>(gzf, *pdata, &(temp[0]), head.bytesPerElement);
#		else
		assertMsg ( ((head.bytesPerElement == sizeof(T)) && (head.elementType==1) ), "pdata type doesn't match");
		IndexInt bytes = sizeof(T)*head.dim;
		IndexInt readBytes = gzread(gzf, &(pdata->get(0)), sizeof(T)*head.dim);
		assertMsg(bytes==readBytes, "can't read uni file, stream length does not match, "<<bytes<<" vs "<<readBytes );
#		endif
	}
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}


//*****************************************************************************
// helper functions


KERNEL(idx) void knQuantize(Grid<Real>& grid, Real step)
{
	int    q  = int(grid(idx) / step + step*0.5);
	double qd = q * (double)step;
	grid(idx) = (Real)qd;
} 
PYTHON() void quantizeGrid(Grid<Real>& grid, Real step) { knQuantize(grid,step); }



// explicit instantiation
template void writeGridRaw<int> (const string& name, Grid<int>*  grid);
template void writeGridRaw<Real>(const string& name, Grid<Real>* grid);
template void writeGridRaw<Vec3>(const string& name, Grid<Vec3>* grid);
template void writeGridUni<int> (const string& name, Grid<int>*  grid);
template void writeGridUni<Real>(const string& name, Grid<Real>* grid);
template void writeGridUni<Vec3>(const string& name, Grid<Vec3>* grid);
template void writeGridVol<int> (const string& name, Grid<int>*  grid);
template void writeGridVol<Vec3>(const string& name, Grid<Vec3>* grid);
template void readGridVol<int> (const string& name, Grid<int>*  grid);
template void readGridVol<Vec3>(const string& name, Grid<Vec3>* grid);
template void writeGridTxt<int> (const string& name, Grid<int>*  grid);
template void writeGridTxt<Real>(const string& name, Grid<Real>* grid);
template void writeGridTxt<Vec3>(const string& name, Grid<Vec3>* grid);
template void readGridRaw<int>  (const string& name, Grid<int>*  grid);
template void readGridRaw<Real> (const string& name, Grid<Real>* grid);
template void readGridRaw<Vec3> (const string& name, Grid<Vec3>* grid);
template void readGridUni<int>  (const string& name, Grid<int>*  grid);
template void readGridUni<Real> (const string& name, Grid<Real>* grid);
template void readGridUni<Vec3> (const string& name, Grid<Vec3>* grid);

template void writePdataUni<int> (const std::string& name, ParticleDataImpl<int>* pdata );
template void writePdataUni<Real>(const std::string& name, ParticleDataImpl<Real>* pdata );
template void writePdataUni<Vec3>(const std::string& name, ParticleDataImpl<Vec3>* pdata );
template void readPdataUni<int>  (const std::string& name, ParticleDataImpl<int>* pdata );
template void readPdataUni<Real> (const std::string& name, ParticleDataImpl<Real>* pdata );
template void readPdataUni<Vec3> (const std::string& name, ParticleDataImpl<Vec3>* pdata );

template void readGrid4dUni<int>  (const string& name, Grid4d<int>*  grid, int readTslice, Grid4d<int>*  slice, void** fileHandle);
template void readGrid4dUni<Real> (const string& name, Grid4d<Real>* grid, int readTslice, Grid4d<Real>* slice, void** fileHandle);
template void readGrid4dUni<Vec3> (const string& name, Grid4d<Vec3>* grid, int readTslice, Grid4d<Vec3>* slice, void** fileHandle);
template void readGrid4dUni<Vec4> (const string& name, Grid4d<Vec4>* grid, int readTslice, Grid4d<Vec4>* slice, void** fileHandle);
template void writeGrid4dUni<int> (const string& name, Grid4d<int>*  grid);
template void writeGrid4dUni<Real>(const string& name, Grid4d<Real>* grid);
template void writeGrid4dUni<Vec3>(const string& name, Grid4d<Vec3>* grid);
template void writeGrid4dUni<Vec4>(const string& name, Grid4d<Vec4>* grid);
template void readGrid4dRaw<int>  (const string& name, Grid4d<int>*  grid);
template void readGrid4dRaw<Real> (const string& name, Grid4d<Real>* grid);
template void readGrid4dRaw<Vec3> (const string& name, Grid4d<Vec3>* grid);
template void readGrid4dRaw<Vec4> (const string& name, Grid4d<Vec4>* grid);
template void writeGrid4dRaw<int> (const string& name, Grid4d<int>*  grid);
template void writeGrid4dRaw<Real>(const string& name, Grid4d<Real>* grid);
template void writeGrid4dRaw<Vec3>(const string& name, Grid4d<Vec3>* grid);
template void writeGrid4dRaw<Vec4>(const string& name, Grid4d<Vec4>* grid);


} //namespace
