/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Base class for all Python-exposed classes
 *
 ******************************************************************************/

// -----------------------------------------------------------------
// NOTE:
// Do not include this file in user code, include "manta.h" instead
// -----------------------------------------------------------------

#ifdef _MANTA_H
#ifndef _PTYPE_H
#define _PTYPE_H

#include <string>
#include <vector>
#include <map>

namespace Manta {
class FluidSolver;

struct PyObject {
	int dummy;
};

struct PbType {
	std::string S;
	std::string str() const;
};

struct PbTypeVec {
	std::vector<PbType> T;
	std::string str() const;
};


//! Base class for all classes
class PbClass {
public:
	PbClass(FluidSolver* parent, const std::string& name="");
	PbClass(const PbClass& a);
	virtual ~PbClass();

	// basic property setter/getters
	void setName(const std::string& name) { mName = name; }
	std::string getName() const { return mName; }
	FluidSolver* getParent() const { return mParent; }
	void setParent(FluidSolver* v) { mParent = v; }
	void checkParent();

	// hidden flag for GUI, debug output
	inline bool isHidden() { return mHidden; }
	inline void setHidden(bool v) { mHidden = v; }

protected:
	FluidSolver* mParent;
	std::string mName;
	bool mHidden;

};

} // namespace

#endif
#endif
