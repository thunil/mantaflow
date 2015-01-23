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
#include  <cstring>

using namespace std;

namespace Manta {

//! uni file header 
typedef struct {
	int dimX, dimY, dimZ; // grid size
	int gridType, elementType, bytesPerElement; // data type info
	char info[256]; // mantaflow build information
	unsigned long long timestamp; // creation time
} UniHeader;

//! in line with grid uni header
typedef struct {
	int dim; // number of partilces
	int dimX, dimY, dimZ; // underlying solver resolution (all data in local coordinates!)
	int elementType, bytesPerElement; // type id and byte size
	char info[256]; // mantaflow build information
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

void writeObjFile(const string& name, Mesh* mesh) {
	errMsg("obj exporter not yet implemented");
}

//*****************************************************************************
// conversion functions for double precision
//*****************************************************************************

#if NO_ZLIB!=1
template <class T>
void gridConvertWrite(gzFile& gzf, Grid<T>& grid, void* ptr, UniHeader& head) {
	errMsg("unknown type, not yet supported");
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

template <class T>
void pdataConvertWrite( gzFile& gzf, ParticleDataImpl<T>& pdata, void* ptr, UniPartHeader& head) {
	errMsg("unknown type, not yet supported");
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

#endif // NO_ZLIB!=1

template <class T>
void gridReadConvert(gzFile& gzf, Grid<T>& grid, void* ptr, int bytesPerElement) {
	errMsg("unknown type, not yet supported");
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
	errMsg("unknown pdata type, not yet supported");
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
	
	int bytes = sizeof(T)*grid->getSizeX()*grid->getSizeY()*grid->getSizeZ();
	int readBytes = gzread(gzf, &((*grid)[0]), bytes);
	assertMsg(bytes==readBytes, "can't read raw file, stream length does not match"<<bytes<<" vs "<<readBytes);
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

//! for test run debugging
PYTHON void printUniFileInfoString(const string& name) {
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (gzf) { 
		char ID[5]={0,0,0,0,0};
		gzread(gzf, ID, 4); 
		if (!strcmp(ID, "MNT2")) {
			UniHeader head;
			assertMsg (gzread(gzf, &head, sizeof(UniHeader)) == sizeof(UniHeader), "can't read file, no header present"); 
			gzclose(gzf);
			debMsg("File '"<<name<<"' info: "<< head.info ,1);
			return; // all good!
		}
		gzclose(gzf);
	}
#	endif
	debMsg("File '"<<name<<"', no valid info string found",1);
}

//! for auto-init & check of results of test runs
void getUniFileSize(const string& name, int& x, int& y, int& z, int* t=NULL) {
	x = y = z = 0;
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "rb");
	if (gzf) { 
		char ID[5]={0,0,0,0,0};
		gzread(gzf, ID, 4); 
		if ( (!strcmp(ID, "MNT2")) || (!strcmp(ID, "M4T2")) ) {
			UniHeader head;
			assertMsg (gzread(gzf, &head, sizeof(UniHeader)) == sizeof(UniHeader), "can't read file, no header present"); 
			x = head.dimX;
			y = head.dimY;
			z = head.dimZ;
		}
		// optionally , read fourth dim
		if ((!strcmp(ID, "M4T2")) && t) {
			int dimT = 0;
			gzread(gzf, &dimT, sizeof(int) );
			(*t) = dimT;
		}
		gzclose(gzf);
	}
#	endif
}
PYTHON Vec3 getUniFileSize(const string& name) {
	int x,y,z;
	getUniFileSize(name, x,y,z);
	return Vec3( Real(x), Real(y), Real(z) );
}


// actual read/write functions

template <class T>
void writeGridUni(const string& name, Grid<T>* grid) {
	debMsg( "writing grid " << grid->getName() << " to uni file " << name ,1);
	
#	if NO_ZLIB!=1
	char ID[5] = "MNT2";
	UniHeader head;
	head.dimX = grid->getSizeX();
	head.dimY = grid->getSizeY();
	head.dimZ = grid->getSizeZ();
	head.gridType = grid->getType();
	head.bytesPerElement = sizeof(T);
	snprintf( head.info, 256, "%s", buildInfoString().c_str() );	
	MuTime stamp; stamp.get();
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
	void* ptr = &((*grid)[0]);
#	if FLOATINGPOINT_PRECISION!=1
	// always write float values, even if compiled with double precision...
	Grid<T> temp(grid->getParent());
	// "misuse" temp grid as storage for floating point values (we have double, so it will always fit)
	gridConvertWrite( gzf, *grid, &(temp[0]), head);
#	endif
	gzwrite(gzf, &head, sizeof(UniHeader));
	gzwrite(gzf, ptr, sizeof(T)*head.dimX*head.dimY*head.dimZ);
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
};


template <class T>
void readGridUni(const string& name, Grid<T>* grid) {
	debMsg( "reading grid " << grid->getName() << " from uni file " << name ,1);

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
		// current file format
		UniHeader head;
		assertMsg (gzread(gzf, &head, sizeof(UniHeader)) == sizeof(UniHeader), "can't read file, no header present");
		assertMsg (head.dimX == grid->getSizeX() && head.dimY == grid->getSizeY() && head.dimZ == grid->getSizeZ(), "grid dim doesn't match, "<< Vec3(head.dimX,head.dimY,head.dimZ)<<" vs "<< grid->getSize() );
		assertMsg ( unifyGridType(head.gridType)==unifyGridType(grid->getType()) , "grid type doesn't match "<< head.gridType<<" vs "<< grid->getType() );
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
	snprintf( head.info, 256, "%s", buildInfoString().c_str() );	
	MuTime stamp; stamp.get();
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
		int bytes     = PartSysSize*head.dim;
		int readBytes = gzread(gzf, &(parts->getData()[0]), bytes);
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
	snprintf( head.info, 256, "%s", buildInfoString().c_str() );	
	MuTime stamp; stamp.get();
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
	if (!gzf) errMsg("can't open file " << name);

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
		int bytes = sizeof(T)*head.dim;
		int readBytes = gzread(gzf, &(pdata->get(0)), sizeof(T)*head.dim);
		assertMsg(bytes==readBytes, "can't read uni file, stream length does not match, "<<bytes<<" vs "<<readBytes );
#		endif
	}
	gzclose(gzf);
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

// explicit instantiation
template void writeGridRaw<int> (const string& name, Grid<int>*  grid);
template void writeGridRaw<Real>(const string& name, Grid<Real>* grid);
template void writeGridRaw<Vec3>(const string& name, Grid<Vec3>* grid);
template void writeGridUni<int> (const string& name, Grid<int>*  grid);
template void writeGridUni<Real>(const string& name, Grid<Real>* grid);
template void writeGridUni<Vec3>(const string& name, Grid<Vec3>* grid);
template void writeGridVol<int> (const string& name, Grid<int>*  grid);
template void writeGridVol<Vec3>(const string& name, Grid<Vec3>* grid);
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

#if ENABLE_GRID_TEST_DATATYPE==1
// dummy functions for test datatype - not really supported right now!
// but we need some function body for linking
template<> void writeGridRaw<nbVector>(const string& name, Grid<nbVector>* grid) {assertMsg(false,"Not supported right now.");};
template<> void writeGridUni<nbVector>(const string& name, Grid<nbVector>* grid) {assertMsg(false,"Not supported right now.");};
template<> void writeGridVol<nbVector>(const string& name, Grid<nbVector>* grid) {assertMsg(false,"Not supported right now.");};
template<> void writeGridTxt<nbVector>(const string& name, Grid<nbVector>* grid) {assertMsg(false,"Not supported right now.");};
template<> void readGridRaw<nbVector> (const string& name, Grid<nbVector>* grid) {assertMsg(false,"Not supported right now.");};
template<> void readGridUni<nbVector> (const string& name, Grid<nbVector>* grid) {assertMsg(false,"Not supported right now.");};
#endif // ENABLE_GRID_TEST_DATATYPE


} //namespace
