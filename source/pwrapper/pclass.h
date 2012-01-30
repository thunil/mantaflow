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

#ifndef _PTYPE_H
#define _PTYPE_H

#include "pconvert.h"
#include <string>
#include <vector>
#include <map>

#ifdef GUI
#   include <QMutex>
#else
struct QMutex {
    void lock() {};
    void unlock() {};
    bool tryLock() {return true;};
};
#endif

// forward declaration to minimize Python.h includes
#ifndef PyObject_HEAD
#ifndef PyObject_Fake
struct _object;
typedef _object PyObject;
#define PyObject_Fake
#endif
#endif

namespace Manta {
class PbClassData;
class FluidSolver;

//! Base class for all classes exposed to Python
class PbClass {
public:
    PbClass(FluidSolver* parent, const std::string& name="", PyObject* obj=NULL);
    virtual ~PbClass();
    
    // basic property setter/getters
    void setName(const std::string& name) { mName = name; }
    std::string getName() const { return mName; }
    PyObject* getPyObject() const { return mPyObject; }
    void setPyObject(PyObject* obj);
    FluidSolver* getParent() { return mParent; }
    void setParent(FluidSolver* v) { mParent = v; }
    void checkParent();
    
    // hidden flag for GUI, debug output
    inline bool isHidden() { return mHidden; }
    inline void setHidden(bool v) { mHidden = v; }
    
    void lock();
    void unlock();
    bool tryLock();
    
    // PbClass instance registry
    static int getNumInstances();
    static PbClass* getInstance(int index);
    
    // converters
    static PbClass* fromPyObject(PyObject* o);
    static bool isNullRef(PyObject* o);
    static PbClass* createPyObject(const std::string& classname, const std::string& name, PbArgs& args, PbClass *parent);
    bool canConvertTo(const std::string& classname);
    
private:
    PbClass& operator=(const PbClass&); // disable copy constructor
    
protected:
    QMutex mMutex;
    FluidSolver* mParent;
    PyObject* mPyObject;
    std::string mName;
    bool mHidden;
        
    friend class PbWrapperRegistry;    
    static std::vector<PbClass*> mInstances;
};

typedef void (*PbInitFunc)(PyObject*);
typedef PyObject* (*PbGenericFunction)(PyObject* self, PyObject* args, PyObject* kwds);
typedef int (*PbConstructor)(PyObject* self, PyObject* args, PyObject* kwds);
typedef PbClass* (*PbClassFactoryFunc)(PbArgs&, PbClass*);
typedef PyObject* (*PbGetter)(PyObject* self, void* closure);
typedef int (*PbSetter)(PyObject* self, PyObject* value, void* closure);

//! Singeton. Registers all classes and methods exposed to Python. 
/*! This class is only used internally by Pb* framwork. 
 *  Please use the functionality of PbClass to lookup and translate pointers. */
class PbWrapperRegistry {
public:
    static PbWrapperRegistry& instance();
    void addClass(const std::string& name, const std::string& internalName, const std::string& baseclass);
    void addExternalInitializer(PbInitFunc func);
    void addMethod(const std::string& classname, const std::string& methodname, PbGenericFunction method);
    void addConstructor(const std::string& classname, PbConstructor method);
    void addAlias(const std::string& classname, const std::string& alias);
    void addGenericFunction(const std::string& funcname, PbGenericFunction func);
    void addInstance(PbClass* obj);
    void addGetSet(const std::string& classname, const std::string& property, PbGetter getfunc, PbSetter setfunc);
    void addPythonPath(const std::string& path);
    void addPythonCode(const std::string& file, const std::string& code);
    PyObject* createPyObject(const std::string& classname, const std::string& name, PbArgs& args, PbClass *parent);
    void construct(const std::string& scriptname);
    void cleanup();
    void renameObjects();
    void runPreInit(std::vector<std::string>& args);
    PyObject* initModule();
    PbClassData* lookup(const std::string& name);
    bool canConvert(PbClassData* from, PbClassData* to);
private:
    PbWrapperRegistry();
    std::map<std::string, PbClassData*> mClasses;
    std::vector<PbClassData*> mClassList;
    std::vector<PbInitFunc> mExtInitializers;
    std::vector<std::string> mPaths;
    std::string mCode, mScriptName;
};

//!\cond Register

void pbFinalizePlugin(FluidSolver* parent, const std::string& name);
void pbPreparePlugin(FluidSolver* parent, const std::string& name);
void pbSetError(const std::string& fn, const std::string& ex);

struct PbRegisterExtInit {
    PbRegisterExtInit(PbInitFunc func) { 
        PbWrapperRegistry::instance().addExternalInitializer(func);        
    }
};

//!\endcond

#define PB_REGISTER_EXTERNAL(name) \
        static const Manta::PbRegisterExtInit _regobj_ext_##name(name);

// Define plugin documentation group
// all plugin functions and classes will automatically be added to this group
//! \defgroup Plugins Plugin functions and classes
//! \defgroup PyClasses Classes exposed to Python
       
} // namespace        

// to avoid the constant incomplete class nagging
#ifndef _PCLASS_NOFLUIDSOLVER
#   include "fluidsolver.h"
#endif

#endif
