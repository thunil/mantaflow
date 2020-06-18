/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2020 Sebastian Barschkis, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * General functions that make use of functions from other io files.
 *
 ******************************************************************************/

#include "mantaio.h"

using namespace std;

namespace Manta {

PYTHON() int load(const string& name, std::vector<PbClass*>& objects, float worldSize=1.0) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));

	if (ext == ".raw")
		return readGridsRaw(name, &objects);
	else if (ext == ".uni")
		return readGridsUni(name, &objects);
	else if (ext == ".vol")
		return readGridsVol(name, &objects);
	if (ext == ".vdb")
		return readObjectsVDB(name, &objects, worldSize);
	else if (ext == ".npz")
		return readGridsNumpy(name, &objects);
	else if (ext == ".txt")
		return readGridsTxt(name, &objects);
	else
		errMsg("file '" + name +"' filetype not supported");
	return 0;
}

PYTHON() int save(const string& name, std::vector<PbClass*>& objects, float worldSize=1.0, bool skipDeletedParts=false) {
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));

	if (ext == ".raw")
		return writeGridsRaw(name, &objects);
	else if (ext == ".uni")
		return writeGridsUni(name, &objects);
	else if (ext == ".vol")
		return writeGridsVol(name, &objects);
	if (ext == ".vdb")
		return writeObjectsVDB(name, &objects, worldSize, skipDeletedParts);
	else if (ext == ".npz")
		return writeGridsNumpy(name, &objects);
	else if (ext == ".txt")
		return writeGridsTxt(name, &objects);
	else
		errMsg("file '" + name +"' filetype not supported");
	return 0;
}

} //namespace
