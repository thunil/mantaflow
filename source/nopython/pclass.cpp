/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL)
 * http://www.gnu.org/licenses
 *
 * Functions for property setting/getting via python
 *
 ******************************************************************************/

#include "manta.h"
#include "general.h"

using namespace std;
namespace Manta {

//******************************************************************************
// Helpers

string PbTypeVec::str() const {
	if (T.empty()) return "";
	string s="<";
	for (int i=0; i<(int)T.size(); i++) {
		s += T[i].str();
		s += (i!=(int)T.size()-1) ? ',' : '>';
	}
	return s;
}
string PbType::str() const {
	if (S=="float") return "Real";
	if (S=="manta.vec3") return "Vec3";
	return S;
}

//******************************************************************************
// PbClass

PbClass::PbClass(FluidSolver* parent, const string& name/*, PyObject* obj*/)
	:  mParent(parent), mName(name), mHidden(false)
{
}

PbClass::PbClass(const PbClass& a) : mParent(a.mParent), mName("_unnamed"), mHidden(false)
{
}

PbClass::~PbClass() {
}
/*
void PbClass::lock() {
	mMutex.lock();
}
void PbClass::unlock() {
	mMutex.unlock();
}
bool PbClass::tryLock() {
	return mMutex.tryLock();
}

PbClass* PbClass::getInstance(int idx) {
	if (idx<0 || idx > (int)mInstances.size())
		errMsg("PbClass::getInstance(): invalid index");
	return mInstances[idx];
}

int PbClass::getNumInstances() {
	return mInstances.size();
}*/
/*
bool PbClass::isNullRef(PyObject* obj) {
	return PyLong_Check(obj) && PyLong_AsDouble(obj)==0;
}

void PbClass::registerObject(PyObject* obj, PbArgs* args) {
	// cross link
	Pb::setReference(this, obj);
	mPyObject = obj;

	mInstances.push_back(this);

	if (args) {
		string _name = args->getOpt<std::string>("name",-1,"");
		if (!_name.empty()) setName(_name);
	}
}

PbClass* PbClass::createPyObject(const string& classname, const string& name, PbArgs& args, PbClass* parent) {
	return Pb::createPy(classname,name,args,parent);
}
*/
void PbClass::checkParent() {
	if (getParent() == NULL) {
		errMsg("New class " + mName + ": no parent given -- specify using parent=xxx !");
	}
}

} // namespace
