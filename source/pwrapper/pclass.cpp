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
#include "pclass.h"
#include "general.h"

using namespace std;
namespace Manta {

const string gDefaultModuleName = "manta";

//******************************************************************************
// Custom object definition

struct PbMethod {
    PbMethod(const string& n, const string& d, PbGenericFunction f) : name(n), doc(d), func(f) {}
    string name, doc;
    PbGenericFunction func;
    
    PyMethodDef def() {
        PyMethodDef def = {&name[0], (PyCFunction)func, METH_VARARGS | METH_KEYWORDS, &doc[0]};
        return def;
    }
};
struct PbGetSet {
    PbGetSet() {}
    PbGetSet(const string& n, const string& d, PbGetter g, PbSetter s) : name(n), doc(d), getter(g), setter(s) {}
    string name, doc;
    PbGetter getter;
    PbSetter setter;

    PyGetSetDef def() {
        PyGetSetDef def = {&name[0], getter, setter, &doc[0], NULL};
        return def;
    }
};

struct PbClassData {
    string pythonName;
    PbInitFunc init;
    PyTypeObject typeInfo;
    //PySequenceMethods seqInfo;
    vector<PbMethod> methods;
    map<string,PbGetSet> getset;
    PbClassData* baseclass;
    PbConstructor constructor;

    vector<PyMethodDef> genMethods;
    vector<PyGetSetDef> genGetSet;
};

struct PbObject {
    PyObject_HEAD
    PbClass *instance;
    PbClassData *classdef;
};

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
    if (name != "FluidSolver::step") {
        std::ostringstream msg;
        if(parent && (parent->getNumInstances()>0) )  msg << parent->getName() << string(".");
        msg << name;
        updateQtGui(false, 0, msg.str() );
    }
    
    // name unnamed PbClass Objects from var name
    PbWrapperRegistry::instance().renameObjects();
}

void pbSetError(const string& fn, const string& ex) {
    cout << "Error in " << fn << endl;
    if (!ex.empty())
        PyErr_SetString(PyExc_RuntimeError, ex.c_str());
}


//******************************************************************************
// Wrapper module def

PbWrapperRegistry::PbWrapperRegistry() {
    addClass("__modclass__", "__modclass__" , "");
    addClass("PbClass", "PbClass", "");
}

static PyObject* zeroGetter(PyObject* self, void* closure) {
    Py_INCREF(Py_None);
    return Py_None;
}

void PbWrapperRegistry::addClass(const string& pythonName, const string& internalName, const string& baseclass) {
    // store class info
    PbClassData* data = new PbClassData;
    data->pythonName = pythonName;
    data->baseclass = baseclass.empty() ? NULL : lookup(baseclass);    
    if (!baseclass.empty() && !data->baseclass)
        errMsg("Registering class '" + internalName + "' : Base class '" + baseclass + "' not found");
    mClasses[internalName] = data;
    mClassList.push_back(data);
    
    // register info object
    addGetSet(internalName, "__mantaclass__", zeroGetter, NULL);
    
    // register all methods of base classes
    if (data->baseclass) {
        for (vector<PbMethod>::iterator it = data->baseclass->methods.begin(); it != data->baseclass->methods.end(); ++it) {
            addMethod(internalName, it->name, it->func);
        }
    }        
}

void PbWrapperRegistry::addExternalInitializer(PbInitFunc func) {
    mExtInitializers.push_back(func);
}

void PbWrapperRegistry::addPythonPath(const string& path) {
    mPaths.push_back(path);
}

void PbWrapperRegistry::addPythonCode(const string& file, const string& code) {
    mCode += code + "\n";
}

void PbWrapperRegistry::addGetSet(const string& classname, const string& property, PbGetter getfunc, PbSetter setfunc) {
    if (mClasses.find(classname) == mClasses.end())
        errMsg("Register class " + classname + " before registering get/setter " + property);
    
    PbClassData* classdef = mClasses[classname];
    PbGetSet& def = classdef->getset[property];
    if (def.name.empty()) {
        def.name = property;
        def.doc = property;
    }
    if (getfunc) def.getter = getfunc;
    if (setfunc) def.setter = setfunc;
}

void PbWrapperRegistry::addAlias(const string& classname, const string& alias) {
    if (mClasses.find(classname) == mClasses.end())
        errMsg("Trying to register alias '" + alias +"' for non-existing class '" + classname + "'");
    mClasses[alias] = mClasses[classname];
}

void PbWrapperRegistry::addMethod(const string& classname, const string& methodname, PbGenericFunction func) {
    string aclass = classname;
    if (aclass.empty())
        aclass = "__modclass__";
    else if (mClasses.find(aclass) == mClasses.end())
        errMsg("Register class '" + aclass + "' before registering methods.");
    
    PbClassData* classdef = mClasses[aclass];
    classdef->methods.push_back(PbMethod(methodname,methodname,func));        
}

void PbWrapperRegistry::addConstructor(const string& classname, PbConstructor func) {
    if (mClasses.find(classname) == mClasses.end())
        errMsg("Register class " + classname + " before registering constructor.");
    
    PbClassData* classdef = mClasses[classname];
    classdef->constructor = func;    
}

PbClassData* PbWrapperRegistry::lookup(const string& name) {
    if (mClasses.find(name) != mClasses.end())
        return mClasses[name];
    
    for(map<string, PbClassData*>::iterator it = mClasses.begin(); it != mClasses.end(); ++it) {
        if (it->second->pythonName == name)
            return it->second;
    }
    return NULL;    
}

void PbWrapperRegistry::cleanup() {
    for(vector<PbClassData*>::iterator it = mClassList.begin(); it != mClassList.end(); ++it) {
        delete *it;
    }
    mClasses.clear();
    mClassList.clear();
}

PbWrapperRegistry& PbWrapperRegistry::instance() {
    static PbWrapperRegistry inst;
    return inst;
}

bool PbWrapperRegistry::canConvert(PbClassData* from, PbClassData* to) {
    if (from == to) return true;
    if (from->baseclass)
        return canConvert(from->baseclass, to);
    return false;
}

//! Assign unnamed PbClass objects their Python variable name
void PbWrapperRegistry::renameObjects() {
    PyObject* sys_mod_dict = PyImport_GetModuleDict();
    PyObject* loc_mod = PyMapping_GetItemString(sys_mod_dict, (char*)"__main__");
    if (!loc_mod) return;
    PyObject* locdict = PyObject_GetAttrString(loc_mod, "__dict__");
    if (!locdict) return;
        
    // iterate all PbClass instances 
    for (size_t i=0; i<PbClass::mInstances.size(); i++) {
        PbClass* obj = PbClass::mInstances[i];
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

//******************************************************************************
// Callback functions



PyMODINIT_FUNC PyInit_Main(void) {
#if PY_MAJOR_VERSION >= 3
    return PbWrapperRegistry::instance().initModule();   
#else
    PbWrapperRegistry::instance().initModule();   
#endif
}

void cbDealloc(PbObject* self) {
    if (self->instance) {
        // don't delete top-level objects
        if (self->instance->getParent() != self->instance)
            delete self->instance;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

PyObject* cbNew(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    PbObject *self = (PbObject*) type->tp_alloc(type, 0);
    if (self != NULL) {
        // lookup and link classdef
        self->classdef = PbWrapperRegistry::instance().lookup(type->tp_name);
        self->instance = NULL;
    } else
        errMsg("can't allocate new python class object");
    return (PyObject*) self;
}

PyObject* PbWrapperRegistry::initModule() {    
    
    // generate and terminate all method lists    
    PyMethodDef sentinelFunc = { NULL, NULL, 0, NULL };
    PyGetSetDef sentinelGetSet = { NULL, NULL, NULL, NULL, NULL };
    for (map<string, PbClassData*>::iterator it = mClasses.begin(); it != mClasses.end(); ++it) {
        it->second->genMethods.clear();
        it->second->genGetSet.clear();
        for (vector<PbMethod>::iterator i2 = it->second->methods.begin(); i2 != it->second->methods.end(); ++i2)
            it->second->genMethods.push_back(i2->def());
        for (map<string,PbGetSet>::iterator i2 = it->second->getset.begin(); i2 != it->second->getset.end(); ++i2)
            it->second->genGetSet.push_back(i2->second.def());
        
        it->second->genMethods.push_back(sentinelFunc);
        it->second->genGetSet.push_back(sentinelGetSet);
    }
    
    // prepare module info
#if PY_MAJOR_VERSION >= 3
    static PyModuleDef MainModule = {
        PyModuleDef_HEAD_INIT,
        gDefaultModuleName.c_str(),
        "Bridge module to the C++ solver",
        -1,
        NULL, NULL, NULL, NULL, NULL
    };
    // get generic methods (plugin functions)
    MainModule.m_methods = &mClasses["__modclass__"]->genMethods[0];
    
    // create module
    PyObject* module = PyModule_Create(&MainModule);
#else
    PyObject* module = Py_InitModule(gDefaultModuleName.c_str(), &mClasses["__modclass__"]->genMethods[0]);
#endif
    
    if (module == NULL)
        return NULL;
    

    // load classes
    for(map<string, PbClassData*>::iterator it = mClasses.begin(); it != mClasses.end(); ++it) {        
        char* nameptr = (char*)it->first.c_str();
        PbClassData& data = *it->second;
        
        // define python classinfo
        PyTypeObject t = {
            PyVarObject_HEAD_INIT(NULL, 0)
            (char*)data.pythonName.c_str(),// tp_name 
            sizeof(PbObject),          // tp_basicsize 
            0,                         // tp_itemsize 
            (destructor)cbDealloc,     // tp_dealloc 
            0,                         // tp_print 
            0,                         // tp_getattr 
            0,                         // tp_setattr 
            0,                         // tp_reserved 
            0,                         // tp_repr 
            0,                         // tp_as_number 
            0,                         // tp_as_sequence 
            0,                         // tp_as_mapping 
            0,                         // tp_hash  
            0,                         // tp_call 
            0,                         // tp_str 
            0,                         // tp_getattro 
            0,                         // tp_setattro 
            0,                         // tp_as_buffer 
            Py_TPFLAGS_DEFAULT | 
            Py_TPFLAGS_BASETYPE,       // tp_flags 
            nameptr,                   // tp_doc 
            0,                         // tp_traverse
            0,                         // tp_clear 
            0,                         // tp_richcompare 
            0,                         // tp_weaklistoffset 
            0,                         // tp_iter 
            0,                         // tp_iternext 
            &data.genMethods[0],   // tp_methods 
            0,                         // tp_members 
            &data.genGetSet[0],    // tp_getset 
            0,                         // tp_base 
            0,                         // tp_dict 
            0,                         // tp_descr_get 
            0,                         // tp_descr_set 
            0,                         // tp_dictoffset 
            (initproc)(data.constructor),// tp_init 
            0,                         // tp_alloc 
            cbNew                     // tp_new 
        };
        data.typeInfo = t;
        
        if (PyType_Ready(&data.typeInfo) < 0)
            continue;
    
        Py_INCREF(&data.typeInfo);
        PyModule_AddObject(module, (char*)data.pythonName.c_str(), (PyObject*) &data.typeInfo);
    }
    
    // externals
    for(vector<PbInitFunc>::iterator it = mExtInitializers.begin(); it != mExtInitializers.end(); ++it) {
        (*it)(module);
    }    
    return module;
}

PyObject* PbWrapperRegistry::createPyObject(const string& classname, const string& name, PbArgs& args, PbClass *parent) {
    PbClassData* classdef = lookup(classname);
    if (!classdef)
        errMsg("Class " + classname + " doesn't exist.");    
    
    // create object
    PyObject* obj = cbNew(&classdef->typeInfo, NULL, NULL);    
    PbObject* self = (PbObject*)obj;
    PyObject* nkw = 0;
    
    if (args.kwds())
        nkw = PyDict_Copy(args.kwds());
    else
        nkw = PyDict_New();
    
    PyObject* nocheck = Py_BuildValue("s","yes");
    PyDict_SetItemString(nkw, "nocheck", nocheck);
    if (parent) PyDict_SetItemString(nkw, "parent", parent->getPyObject());
    
    // create instance
    if (self->classdef->constructor(obj, args.linArgs(), nkw) < 0)
        errMsg("error raised in constructor"); // assume condition is already set
    
    Py_DECREF(nkw);
    Py_DECREF(nocheck);
    self->instance->setName(name);
        
    return obj;    
}
 
void PbWrapperRegistry::addInstance(PbClass* obj) {
    PbClass::mInstances.push_back(obj);
}

// prepare typeinfo and register python module
void PbWrapperRegistry::construct(const string& scriptname) {
    mScriptName = scriptname;
    
    // load main extension module
    PyImport_AppendInittab(gDefaultModuleName.c_str(), PyInit_Main);
}

void PbWrapperRegistry::runPreInit(vector<string>& args) {
    // add python directories to path
    PyObject *sys_path = PySys_GetObject((char*)"path");
    for (size_t i=0; i<mPaths.size(); i++) {
        PyObject *path = toPy(mPaths[i]);
        if (sys_path == NULL || path == NULL || PyList_Append(sys_path, path) < 0) {
            errMsg("unable to set python path");
        }
        Py_DECREF(path);
    }
    
    // run preinit code
    
    // provide arguments
    string iargs = "args = ['";
    for (size_t i=1; i<args.size(); i++) {
        if (i>1) iargs+="', '";
        iargs+=args[i];
    }
    iargs += "']\n";
    PyRun_SimpleString(iargs.c_str());
    
    // provide compile flags
    string cfl = "";
#ifdef CUDA
    cfl+="CUDA=True\n";
#else
    cfl+="CUDA=False\n";
#endif
#ifdef DEBUG
    cfl+="DEBUG=True\n";
#else
    cfl+="DEBUG=False\n";
#endif
#ifdef MT
    cfl+="MT=True\n";
#else
    cfl+="MT=False\n";
#endif
#ifdef GUI
    cfl+="GUI=True\n";
#else
    cfl+="GUI=False\n";
#endif
    cfl+="SCENEFILE='"+mScriptName+"'\n";
    PyRun_SimpleString(cfl.c_str());
    
    if (!mCode.empty())
        PyRun_SimpleString(mCode.c_str());
}

//******************************************************************************
// PbClass

PbClass::PbClass(FluidSolver* parent, const string& name, PyObject* obj)
    : mName(name), mPyObject(obj), mParent(parent), mHidden(false), mMutex()
{
}

PbClass::PbClass(const PbClass& a) : mName("_unnamed"), mPyObject(0), mParent(a.mParent), mHidden(false), mMutex()
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
    
PbClass* PbClass::fromPyObject(PyObject* o) {
    if(!PyObject_HasAttrString(o, "__mantaclass__"))
        return NULL;
        
    return ((PbObject*) o)->instance;
}

bool PbClass::isNullRef(PyObject* obj) {
    return PyLong_Check(obj) && PyLong_AsDouble(obj)==0;
}

void PbClass::setPyObject(PyObject* obj) {
    if(!PyObject_HasAttrString(obj, "__mantaclass__"))
        return;
    
    // cross link
    ((PbObject*) obj)->instance = this;
    mPyObject = obj;
    
    mInstances.push_back(this);
}

PbClass* PbClass::createPyObject(const string& classname, const string& name, PbArgs& args, PbClass* parent) {
    PyObject* obj = PbWrapperRegistry::instance().createPyObject(classname, name, args, parent);
    return ((PbObject*)obj)->instance;
}

bool PbClass::canConvertTo(const string& classname) {
    PbClassData* from = ((PbObject*)mPyObject)->classdef;
    PbClassData* dest = PbWrapperRegistry::instance().lookup(classname);
    if (!dest)
        errMsg("Classname '" + classname + "' is not registered.");
    return PbWrapperRegistry::instance().canConvert(from, dest);
}

void PbClass::checkParent() {
    if (getParent() == NULL) {
        errMsg("New class " + mName + ": no parent given -- specify using parent=xxx !");
    }
}

PyObject* PbClass::assignNewPyObject(const string& classname) {
    PbClassData* classdef = PbWrapperRegistry::instance().lookup(classname);
    assertMsg(classdef,"python class " + classname + " does not exist.");
    
    // allocate new object
    PbObject *obj = (PbObject*) classdef->typeInfo.tp_alloc(&(classdef->typeInfo), 0);
    assertMsg(obj, "cannot allocate new python object");
    
    obj->classdef = classdef;
    setPyObject((PyObject*)obj);
    
    return mPyObject;
}


vector<PbClass*> PbClass::mInstances;

} // namespace