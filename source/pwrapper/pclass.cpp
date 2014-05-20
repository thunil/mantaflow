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

#include "pythonInclude.h"
#include "structmember.h"
#include "manta.h"
#include "general.h"

using namespace std;
namespace Manta {

//******************************************************************************
// Free functions

#ifdef GUI
    extern void updateQtGui(bool full, int frame, const std::string& curPlugin);
#else
    inline void updateQtGui(bool full, int frame, const std::string& curPlugin) {}
#endif

void pbPreparePlugin(FluidSolver* parent, const string& name) {
    if (parent)
        parent->pluginStart(name);
}

void pbFinalizePlugin(FluidSolver *parent, const string& name) {
    if (parent) {
        parent->pluginStop(name);
    }
    
    // GUI update, also print name of parent if there's more than one
    std::ostringstream msg;
    if (name != "FluidSolver::step") 
	{
        if(parent && (parent->getNumInstances()>0) )  msg << parent->getName() << string(".");
        msg << name;
    }
    updateQtGui(false, 0, msg.str() );
    
    // name unnamed PbClass Objects from var name
    PbClass::renameObjects();
}

void pbSetError(const string& fn, const string& ex) {
    cout << "Error in " << fn << endl;
    if (!ex.empty())
        PyErr_SetString(PyExc_RuntimeError, ex.c_str());
}
 
//******************************************************************************
// PbClass

vector<PbClass*> PbClass::mInstances;

PbClass::PbClass(FluidSolver* parent, const string& name, PyObject* obj)
    : mMutex(), mParent(parent), mPyObject(obj), mName(name), mHidden(false)
{
}

PbClass::PbClass(const PbClass& a) : mMutex(), mParent(a.mParent), mPyObject(0), mName("_unnamed"), mHidden(false)
{
}
    

PbClass::~PbClass() {   
    for(vector<PbClass*>::iterator it = mInstances.begin(); it != mInstances.end(); ++it) {
        if (*it == this) {
            mInstances.erase(it);
            break;
        }
    }
}
    
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
}
    
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

void PbClass::checkParent() {
    if (getParent() == NULL) {
        errMsg("New class " + mName + ": no parent given -- specify using parent=xxx !");
    }
}
//! Assign unnamed PbClass objects their Python variable name
void PbClass::renameObjects() {
    PyObject* sys_mod_dict = PyImport_GetModuleDict();
    PyObject* loc_mod = PyMapping_GetItemString(sys_mod_dict, (char*)"__main__");
    if (!loc_mod) return;
    PyObject* locdict = PyObject_GetAttrString(loc_mod, "__dict__");
    if (!locdict) return;
        
    // iterate all PbClass instances 
    for (size_t i=0; i<mInstances.size(); i++) {
        PbClass* obj = mInstances[i];
        if (obj->getName().empty()) {
            // empty, try to find instance in module local dictionary
            
            PyObject *lkey, *lvalue;
            Py_ssize_t lpos = 0;
            while (PyDict_Next(locdict, &lpos, &lkey, &lvalue)) {
                if (lvalue == obj->mPyObject) {
                    string varName = fromPy<string>(PyObject_Str(lkey));
                    obj->setName(varName);
                    //cout << "assigning variable name '" << varName << "' to unnamed instance" << endl;
                    break;
                }
            }
        }
    }
    Py_DECREF(locdict);
    Py_DECREF(loc_mod);    
}


} // namespace
