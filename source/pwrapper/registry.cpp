/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Auto python registry
 *
 ******************************************************************************/

#include "pythonInclude.h"
#include "structmember.h"
#include "manta.h"

using namespace std;

const string gDefaultModuleName = "manta";

namespace Pb {

//******************************************************************************
// Custom object definition

struct Method {
    Method(const string& n, const string& d, GenericFunction f) : name(n), doc(d), func(f) {}
    string name, doc;
    GenericFunction func;
    
    PyMethodDef def() {
        PyMethodDef def = {&name[0], (PyCFunction)func, METH_VARARGS | METH_KEYWORDS, &doc[0]};
        return def;
    }
};
struct GetSet {
    GetSet() {}
    GetSet(const string& n, const string& d, Getter g, Setter s) : name(n), doc(d), getter(g), setter(s) {}
    string name, doc;
    Getter getter;
    Setter setter;

    PyGetSetDef def() {
        PyGetSetDef def = {&name[0], getter, setter, &doc[0], NULL};
        return def;
    }
};

struct ClassData {
    string cName, pyName;
    InitFunc init;
    PyTypeObject typeInfo;
    //PySequenceMethods seqInfo;
    vector<Method> methods;
    map<string,GetSet> getset;
    ClassData* baseclass;
    string baseclassName;
    Constructor constructor;

    vector<PyMethodDef> genMethods;
    vector<PyGetSetDef> genGetSet;
};

struct PbObject {
    PyObject_HEAD
    Manta::PbClass *instance;
    ClassData *classdef;
};

//******************************************************
// Internal wrapper class

//! Registers all classes and methods exposed to Python. 
/*! This class is only used internally by Pb:: framwork. 
 *  Please use the functionality of PbClass to lookup and translate pointers. */
class WrapperRegistry {
public:
    static WrapperRegistry& instance();
    void addClass(const std::string& name, const std::string& internalName, const std::string& baseclass);
    void addExternalInitializer(InitFunc func);
    void addMethod(const std::string& classname, const std::string& methodname, GenericFunction method);
    void addConstructor(const std::string& classname, Constructor method);
    void addGetSet(const std::string& classname, const std::string& property, Getter getfunc, Setter setfunc);
    void addPythonPath(const std::string& path);
    void addPythonCode(const std::string& file, const std::string& code);
    PyObject* createPyObject(const std::string& classname, const std::string& name, Manta::PbArgs& args, Manta::PbClass *parent);
    void construct(const std::string& scriptname);
    void cleanup();
    void renameObjects();
    void runPreInit(const std::vector<std::string>& args);
    PyObject* initModule();
    ClassData* lookup(const std::string& name);
    bool canConvert(ClassData* from, ClassData* to);
private:
    ClassData* getOrConstructClass(const string& name);
    void registerBaseclasses();
    void registerAliases();
    void addParentMethods(ClassData* cls, ClassData* base);
    WrapperRegistry();
    std::map<std::string, ClassData*> mClasses;
    std::vector<ClassData*> mClassList;
    std::vector<InitFunc> mExtInitializers;
    std::vector<std::string> mPaths;
    std::string mCode, mScriptName;
};

//******************************************************************************
// Callback functions

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
        self->classdef = WrapperRegistry::instance().lookup(type->tp_name);
        self->instance = NULL;
    } else
        errMsg("can't allocate new python class object");
    return (PyObject*) self;
}

PyMODINIT_FUNC PyInit_Main(void) {
#if PY_MAJOR_VERSION >= 3
    return WrapperRegistry::instance().initModule();   
#else
    WrapperRegistry::instance().initModule();   
#endif
}

//******************************************************
// WrapperRegistry

WrapperRegistry::WrapperRegistry() {
    addClass("__modclass__", "__modclass__" , "");
    addClass("PbClass", "PbClass", "");
}

ClassData* WrapperRegistry::getOrConstructClass(const string& classname) {
    map<string,ClassData*>::iterator it = mClasses.find(classname);

    if (it != mClasses.end())
        return it->second;
    ClassData* data = new ClassData;
    data->cName = classname;
    data->baseclass = NULL;
    mClasses[classname] = data;
    mClassList.push_back(data);
    return data;
}

void replaceAll(string& source, string const& find, string const& replace) {
    for(string::size_type i = 0; (i = source.find(find, i)) != std::string::npos;)
    {
        source.replace(i, find.length(), replace);
        i += replace.length() - find.length() + 1;
    }
}

void WrapperRegistry::addClass(const string& pyName, const string& internalName, const string& baseclass) {
    ClassData* data = getOrConstructClass(internalName);
    
    // regularize python name
    string pythonName = pyName;
    replaceAll(pythonName, "<", "_");
    replaceAll(pythonName, ">", "");
    replaceAll(pythonName, ",", "_");
    //cout << pythonName << endl;
    
    if (data->pyName.empty()) 
        data->pyName = pythonName;
    mClasses[pythonName] = data;
    if (!baseclass.empty())
        data->baseclassName = baseclass;    
}

void WrapperRegistry::addExternalInitializer(InitFunc func) {
    mExtInitializers.push_back(func);
}

void WrapperRegistry::addPythonPath(const string& path) {
    mPaths.push_back(path);
}

void WrapperRegistry::addPythonCode(const string& file, const string& code) {
    mCode += code + "\n";
}

void WrapperRegistry::addGetSet(const string& classname, const string& property, Getter getfunc, Setter setfunc) {
    ClassData* classdef = getOrConstructClass(classname);
    GetSet& def = classdef->getset[property];
    if (def.name.empty()) {
        def.name = property;
        def.doc = property;
    }
    if (getfunc) def.getter = getfunc;
    if (setfunc) def.setter = setfunc;
}

void WrapperRegistry::addMethod(const string& classname, const string& methodname, GenericFunction func) {
    string aclass = classname;
    if (aclass.empty())
        aclass = "__modclass__";
    
    ClassData* classdef = getOrConstructClass(aclass);
    for(int i=0; i<classdef->methods.size(); i++)
        if (classdef->methods[i].name == methodname) return; // avoid duplicates
    classdef->methods.push_back(Method(methodname,methodname,func)); 
}

void WrapperRegistry::addConstructor(const string& classname, Constructor func) {
    ClassData* classdef = getOrConstructClass(classname);
    classdef->constructor = func;
}

void WrapperRegistry::addParentMethods(ClassData* cur, ClassData* base) {
    if (base == 0) return;

    for (vector<Method>::iterator it = base->methods.begin(); it != base->methods.end(); ++it)
        addMethod(cur->cName, it->name, it->func);

    for (map<string,GetSet>::iterator it = base->getset.begin(); it != base->getset.end(); ++it)
        addGetSet(cur->cName, it->first, it->second.getter, it->second.setter);

    addParentMethods(cur, base->baseclass);
}

void WrapperRegistry::registerBaseclasses() {
    for (int i=0; i<mClassList.size(); i++) {
        string bname = mClassList[i]->baseclassName;
        if(!bname.empty()) {
            mClassList[i]->baseclass = lookup(bname);
            if (!mClassList[i]->baseclass)
                errMsg("Registering class '" + mClassList[i]->cName + "' : Base class '" + bname + "' not found");
        }
    }

    for (int i=0; i<mClassList.size(); i++) {
        addParentMethods(mClassList[i], mClassList[i]->baseclass);
    }
}

void WrapperRegistry::registerAliases() {
    for (map<string, ClassData*>::iterator it = mClasses.begin(); it != mClasses.end(); ++it) {
        if (it->first != it->second->pyName && it->first.find('<') == string::npos) {
            mCode += it->first + " = " + it->second->pyName + "\n";
        }
    }
}

ClassData* WrapperRegistry::lookup(const string& name) {
    for(map<string, ClassData*>::iterator it = mClasses.begin(); it != mClasses.end(); ++it) {
        if (it->first == name || it->second->cName == name)
            return it->second;
    }
    return NULL;    
}

void WrapperRegistry::cleanup() {
    for(vector<ClassData*>::iterator it = mClassList.begin(); it != mClassList.end(); ++it) {
        delete *it;
    }
    mClasses.clear();
    mClassList.clear();
}

WrapperRegistry& WrapperRegistry::instance() {
    static WrapperRegistry inst;
    return inst;
}

bool WrapperRegistry::canConvert(ClassData* from, ClassData* to) {
    if (from == to) return true;
    if (from->baseclass)
        return canConvert(from->baseclass, to);
    return false;
}

void WrapperRegistry::runPreInit(const vector<string>& args) {
    // add python directories to path
    PyObject *sys_path = PySys_GetObject((char*)"path");
    for (size_t i=0; i<mPaths.size(); i++) {
        PyObject *path = Manta::toPy(mPaths[i]);
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
    
    if (!mCode.empty()) {
        mCode = "from manta import *\n" + mCode;
        PyRun_SimpleString(mCode.c_str());
    }
}

PyObject* WrapperRegistry::createPyObject(const string& classname, const string& name, Manta::PbArgs& args, Manta::PbClass *parent) {
    ClassData* classdef = lookup(classname);
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

// prepare typeinfo and register python module
void WrapperRegistry::construct(const string& scriptname) {
    mScriptName = scriptname;

    registerBaseclasses();
    registerAliases();
    
    // load main extension module
    PyImport_AppendInittab(gDefaultModuleName.c_str(), PyInit_Main);
}

PyObject* WrapperRegistry::initModule() {
    // generate and terminate all method lists    
    PyMethodDef sentinelFunc = { NULL, NULL, 0, NULL };
    PyGetSetDef sentinelGetSet = { NULL, NULL, NULL, NULL, NULL };
    for (int i=0; i<mClassList.size(); i++) {
        ClassData* cls = mClassList[i];
        cls->genMethods.clear();
        cls->genGetSet.clear();
        for (vector<Method>::iterator i2 = cls->methods.begin(); i2 != cls->methods.end(); ++i2)
            cls->genMethods.push_back(i2->def());
        for (map<string,GetSet>::iterator i2 = cls->getset.begin(); i2 != cls->getset.end(); ++i2)
            cls->genGetSet.push_back(i2->second.def());
        
        cls->genMethods.push_back(sentinelFunc);
        cls->genGetSet.push_back(sentinelGetSet);
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
    for(map<string, ClassData*>::iterator it = mClasses.begin(); it != mClasses.end(); ++it) {        
        char* nameptr = (char*)it->first.c_str();
        ClassData& data = *it->second;
        
        // define python classinfo
        PyTypeObject t = {
            PyVarObject_HEAD_INIT(NULL, 0)
            (char*)data.pyName.c_str(),// tp_name 
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
        PyModule_AddObject(module, (char*)data.pyName.c_str(), (PyObject*) &data.typeInfo);
    }
    
    // externals
    for(vector<InitFunc>::iterator it = mExtInitializers.begin(); it != mExtInitializers.end(); ++it) {
        (*it)(module);
    }    
    return module;
}


//******************************************************
// Register members and exposed functions

void setup(const std::string& filename, const std::vector<std::string>& args) {
    WrapperRegistry::instance().construct(filename);
    Py_Initialize();
    WrapperRegistry::instance().runPreInit(args);
}

void finalize() {
    Py_Finalize();
    WrapperRegistry::instance().cleanup();
}

bool canConvert(PyObject* obj, const string& classname) {
    ClassData* from = ((PbObject*)obj)->classdef;
    ClassData* dest = WrapperRegistry::instance().lookup(classname);
    if (!dest)
        errMsg("Classname '" + classname + "' is not registered.");
    return WrapperRegistry::instance().canConvert(from, dest);
}

Manta::PbClass* objFromPy(PyObject* obj) {
    if (obj->ob_type->tp_dealloc != (destructor)cbDealloc) // not a manta object
        return NULL;
        
    return ((PbObject*) obj)->instance;
}

PyObject* copyObject(Manta::PbClass* cls, const string& classname) {
    ClassData* classdef = WrapperRegistry::instance().lookup(classname);
    assertMsg(classdef,"python class " + classname + " does not exist.");
    
    // allocate new object
    PbObject *obj = (PbObject*) classdef->typeInfo.tp_alloc(&(classdef->typeInfo), 0);
    assertMsg(obj, "cannot allocate new python object");
    
    obj->classdef = classdef;
    cls->registerObject((PyObject*)obj, 0);
    
    return cls->getPyObject();
}

Manta::PbClass* createPy(const std::string& classname, const std::string& name, Manta::PbArgs& args, Manta::PbClass* parent) {
    PyObject* obj = WrapperRegistry::instance().createPyObject(classname, name, args, parent);
    return ((PbObject*)obj)->instance;
}

void setReference(Manta::PbClass* cls, PyObject* obj) {
    ((PbObject*) obj)->instance = cls;
}

Register::Register(const string& className, const string& funcName, GenericFunction func) {
    WrapperRegistry::instance().addMethod(className, funcName, func);
}
Register::Register(const string& className, const string& funcName, Constructor func) {
    WrapperRegistry::instance().addConstructor(className, func);
}
Register::Register(const string& className, const string& property, Getter getter, Setter setter) {
    WrapperRegistry::instance().addGetSet(className, property, getter, setter);
}
Register::Register(const string& className, const string& pyName, const string& baseClass) {
    WrapperRegistry::instance().addClass(pyName, className, baseClass);
}
Register::Register(const string& file, const string& pythonCode) {
    WrapperRegistry::instance().addPythonCode(file, pythonCode);
}
Register::Register(InitFunc func) { 
    WrapperRegistry::instance().addExternalInitializer(func);        
}

} // namespace