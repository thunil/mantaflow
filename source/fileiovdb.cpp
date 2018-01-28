
using namespace std;
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

#include "openvdb/openvdb.h"

namespace Manta {

	template <class T>
	void writeGridVDB(const string& name, Grid<T>* grid) {
	
	}
	template <>
	void writeGridVDB(const string& name, Grid<Real>* grid) {
		debMsg("Writing real grid " << grid->getName() << " to vdb file " << name, 1);

		openvdb::initialize();
		
		// Create an empty floating-point grid with background value 0.
		openvdb::FloatGrid::Ptr gridVDB = openvdb::FloatGrid::create();

		gridVDB->setTransform(
			openvdb::math::Transform::createLinearTransform(/*voxel size=*/0.08));

		// Get an accessor for coordinate-based access to voxels.
		openvdb::FloatGrid::Accessor accessor = gridVDB->getAccessor();

		// Identify the grid as a level set.
		gridVDB->setGridClass(openvdb::GRID_FOG_VOLUME);

		// Name the grid "density".
		gridVDB->setName("density");

		openvdb::io::File file(name);

		FOR_IJK(*grid) {
			
			openvdb::Coord xyz(i, j, k);
			accessor.setValue(xyz, (*grid)(i, j, k));	
		}

		// Add the grid pointer to a container.
		openvdb::GridPtrVec gridsVDB;
		gridsVDB.push_back(gridVDB);

		// Write out the contents of the container.
		file.write(gridsVDB);
		file.close();
	};

	template void writeGridVDB<int>(const string& name, Grid<int>*  grid);
	template void writeGridVDB<Vec3>(const string& name, Grid<Vec3>* grid);
	template void writeGridVDB<Real>(const string& name, Grid<Real>* grid);

}