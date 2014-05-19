/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor: Process replacement text of PYTHON keywords
 *
 ******************************************************************************/

#include "prep.h"
#include <cstdlib>
#include <set>
#include <sstream>
#include <iostream>
using namespace std;

#define STR(x) #x

//******************************************************
// Templates for code generation

const string TmpFunction = STR(
static PyObject* _P_$FUNCNAME$ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
    try {
        PbArgs _args(_linargs, _kwds);
        FluidSolver *parent = _args.obtainParent();
        pbPreparePlugin(parent, "$FUNCNAME$" );
        PyObject *_retval = 0;
        { 
            ArgLocker _lock;
            $ARGLOADER$
            @IF(RET_VOID)
                _retval = getPyNone();
                $FUNCNAME$($CALLSTRING$);
            @ELSE
                _retval = toPy($FUNCNAME$($CALLSTRING$));
            @END
            _args.check();
        }
        pbFinalizePlugin(parent,"$FUNCNAME$" );
        return _retval;
    } catch(std::exception& e) {
        pbSetError("$FUNCNAME$",e.what());
        return 0;
    }
}
static const Pb::Register _RP_$FUNCNAME$ ("","$FUNCNAME$",_P_$FUNCNAME$);
);

const string TmpMemberFunction = STR(
static PyObject* _$FUNCNAME$ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
    try {
        PbArgs _args(_linargs, _kwds);
        $CLASS$* pbo = dynamic_cast<$CLASS$*>(Pb::objFromPy(_self));
        pbPreparePlugin(pbo->getParent(), "$CLASS$::$FUNCNAME$");
        PyObject *_retval = 0;
        { 
            ArgLocker _lock;
            $ARGLOADER$
            pbo->_args.copy(_args);
            @IF(RET_VOID)
                _retval = getPyNone();
                pbo->$FUNCNAME$($CALLSTRING$);
            @ELSE
                _retval = toPy(pbo->$FUNCNAME$($CALLSTRING$));
            @END
            pbo->_args.check();
        }
        pbFinalizePlugin(pbo->getParent(),"$CLASS$::$FUNCNAME$");
        return _retval;
    } catch(std::exception& e) {
        pbSetError("$CLASS$::$FUNCNAME$",e.what());
        return 0;
    }
});

const string TmpConstructor = STR(
static int _$CLASS$ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
    PbClass* obj = Pb::objFromPy(_self); 
    if (obj) delete obj; 
    try {
        PbArgs _args(_linargs, _kwds);
        pbPreparePlugin(0, "$CLASS$::$FUNCNAME$" );
        { 
            ArgLocker _lock;
            $ARGLOADER$
            obj = new $CLASS$($CALLSTRING$);
            obj->registerObject(_self, &_args); 
            _args.check();
        }
        pbFinalizePlugin(obj->getParent(),"$CLASS$::$FUNCNAME$" );
        return 0;
    } catch(std::exception& e) {
        pbSetError("$CLASS$::$FUNCNAME$",e.what());
        return -1;
    }
});

const string TmpGetSet = STR(
static PyObject* _GET_$NAME$(PyObject* self, void* cl) {
    $CLASS$* pbo = dynamic_cast<$CLASS$*>(Pb::objFromPy(self));
    return toPy(pbo->$NAME$);
}
static int _SET_$NAME$(PyObject* self, PyObject* val, void* cl) {
    $CLASS$* pbo = dynamic_cast<$CLASS$*>(Pb::objFromPy(self));
    pbo->$NAME$ = fromPy<$TYPE$ >(val); 
    return 0;
});

const string TmpRegisterMethod = STR(
@IF(CTPL)
    static const Pb::Register _R_$IDX$ ("$CLASS$<$CT$>","$FUNCNAME$",$CLASS$<$CT$>::_$FUNCNAME$);
@ELSE
    static const Pb::Register _R_$IDX$ ("$CLASS$","$FUNCNAME$",$CLASS$::_$FUNCNAME$);
@END
);

const string TmpRegisterGetSet = STR(
@IF(CTPL)
    static const Pb::Register _R_$IDX$ ("$CLASS$<$CT$>$","$PYNAME$",$CLASS$<$CT$>::_GET_$NAME$,$CLASS$<$CT$>::_SET_$NAME$);
@ELSE
    static const Pb::Register _R_$IDX$ ("$CLASS$","$PYNAME$",$CLASS$::_GET_$NAME$,$CLASS$::_SET_$NAME$);
@END
);

const string TmpRegisterClass = STR(
@IF(CTPL)
    static const Pb::Register _R_$IDX$ ("$CLASS$<$CT$>","$PYNAME$<$CT$>","$BASE$$BTPL$");
    template<> const char* Namify<$CLASS$<$CT$> >::S = "$CLASS$<$CT$>";
@ELSE
    static const Pb::Register _R_$IDX$ ("$CLASS$","$PYNAME$","$BASE$$BTPL$");
    template<> const char* Namify<$CLASS$ >::S = "$CLASS$";
@END
);

const string TmpAlias = STR(
static const Pb::Register _R_$IDX$ ("$CLASS$","$PYNAME$","");
);

//******************************************************
// Code generation functions

string generateLoader(const Argument& arg) {
    bool integral = isIntegral(arg.type.name);
    Type ptrType = arg.type;
    if (integral) {
        ptrType.isPointer = false;
        ptrType.isRef = false;
        ptrType.isConst = false;
    } else if (arg.type.isRef) {
        ptrType.isPointer = true;
        ptrType.isRef = false;
        ptrType.isConst = false;
    }
    
    stringstream loader;
    loader << arg.type.minimal << " " << arg.name << " = ";
    if (arg.type.isRef && !integral)
        loader << "*";

    loader << "_args." << (arg.value.empty() ? "get" : "getOpt");
    loader << "<" << ptrType.build() << " >";
    loader << "(" << arg.index << ",\"" << arg.name << "\",";
    if (!arg.value.empty())
        loader << arg.value << ",";
    loader << "&_lock); ";

    return loader.str();
}

// global for tracking state between python class and python function registrations
bool gFoundConstructor = false;

void processPythonFunction(const Block& block, const string& code, Sink& sink) {
    const Function& func = block.func;

    // PYTHON(...) keyword options
    for (size_t i=0; i<block.options.size(); i++) {
        errMsg(block.line0, "unknown keyword " + block.options[i].name);    
    }

    bool isConstructor = func.returnType.minimal.empty();
    bool isPlugin = !block.parent;
    if (isConstructor) gFoundConstructor = true;
    if (isPlugin && sink.isHeader)
        errMsg(block.line0,"plugin python functions can't be defined in headers.");

    // replicate function
    const string signature = func.minimal + block.initList;
    if (gDocMode) {
        // document free plugins
        if (isPlugin)
            sink.inplace << "//! \\ingroup Plugins\n";
        sink.inplace << "PYTHON " << signature << "{}\n";
        return;
    }
    sink.inplace << block.linebreaks() << signature << code;

    // generate variable loader
    string loader = "";
    for (int i=0; i<func.arguments.size(); i++)
        loader += generateLoader(func.arguments[i]);

    // generate glue layer function
    const string table[] = { "FUNCNAME", func.name, 
                             "ARGLOADER", loader, 
                             "CLASS", isPlugin ? "" : block.parent->name,
                             "CTPL", (isPlugin || !block.parent->isTemplated()) ? "" : "$CT$",
                             "CALLSTRING", func.callString(),
                             "RET_VOID", (func.returnType.name=="void") ? "Y" : "",
                             "@end" };
    string callerTempl = isConstructor ? TmpConstructor : (isPlugin ? TmpFunction : TmpMemberFunction);
    sink.inplace << replaceSet(callerTempl, table);

    // register member functions
    if (!isPlugin) {
        const string reg = replaceSet(TmpRegisterMethod, table);
        sink.link << '+' << block.parent->name << '^' << reg << '\n';
    }
}

void processPythonVariable(const Block& block, Sink& sink) {
    const Function& var = block.func;

    if (!block.parent)
        errMsg(block.line0, "python variables can only be used inside classes");

    // process options
    string pythonName = var.name;
    for (size_t i=0; i<block.options.size(); i++) {
        if (block.options[i].name == "name") 
            pythonName = block.options[i].value;
        else
            errMsg(block.line0, "PYTHON(opt): illegal option. Supported options are: 'name'");
    }

    // generate glue layer function
    const string table[] = { "NAME", var.name, 
                             "CLASS", block.parent->name,
                             "CTPL", !block.parent->isTemplated() ? "" : "$CT$",
                             "PYNAME", pythonName,
                             "TYPE", var.returnType.minimal,
                             "@end" };

    // output function and accessors
    sink.inplace << block.linebreaks() << var.minimal << ";";
    sink.inplace << replaceSet(TmpGetSet, table);

    // register accessors
    const string reg = replaceSet(TmpRegisterGetSet, table);
    sink.link << '+' << block.parent->name << '^' << reg << '\n';
}

void processPythonClass(const Block& block, const string& code, Sink& sink) {
    const Class& cls = block.cls;
    string pythonName = cls.name;

    if (!sink.isHeader)
        errMsg(block.line0, "PYTHON classes can only be defined in header files.");
    
    // PYTHON(...) keyword options
    for (size_t i=0; i<block.options.size(); i++) {
        if (block.options[i].name == "name") 
            pythonName = block.options[i].value;
        else
            errMsg(block.line0, "PYTHON(opt): illegal kernel option. Supported options are: 'name'");
    }
    
    if (gDocMode) {
        sink.inplace << "//! \\ingroup PyClasses\nPYTHON " << cls.minimal;
        return;
    }

    // register class
    const string table[] = { "CLASS", cls.name,
                             "BASE", cls.baseClass.name,
                             "BTPL", cls.baseClass.isTemplated() ? "<$BT$>" : "",
                             "PYNAME", pythonName,
                             "CTPL", cls.isTemplated() ? "CT" : "",
                             "@end" };

    // register class
    string reg = replaceSet(TmpRegisterClass, table);
    sink.link << '+' << cls.name << '^' << reg << '\n';
    // instantiate directly if not templated
    if (!cls.isTemplated())
        sink.link << '>' << cls.name << "^\n";
    // chain the baseclass instantiation
    if (cls.baseClass.isTemplated())
        sink.link << '@' << cls.name << '^' << cls.templateTypes.names() << '^' 
                  << cls.baseClass.name << '^' << cls.baseClass.templateTypes.names() << '\n';

    // write signature
    sink.inplace << block.linebreaks() << cls.minimal << "{";

    // remove first {, and steal two linebreaks so we can add a #define later
    string ncode = code.substr(1);
    stealLinebreaks(ncode, 2);
    
    // scan code for member functions
    gFoundConstructor = false;
    processText(ncode.substr(1), block.line0, sink, &cls);
    if (!gFoundConstructor)
        errMsg(block.line0, "no PYTHON constructor found in class '" + cls.name + "'");
    
    // add secret bonus members to class and close
    sink.inplace << "public: PbArgs _args;";
    sink.inplace << "}\n";
    // add a define to make commenting out classes, and #ifdefs work correctly
    sink.inplace << "#define _C_" << cls.name << '\n';
}

void processPythonInstantiation(const Block& block, const Type& aliasType, const string& aliasName, Sink& sink) {
    if (!sink.isHeader)
        errMsg(block.line0, "instantiate allowed in headers only");

    const string table[] = {"CLASS", strip(aliasType.build()), "PYNAME", aliasName, "@end"};
    sink.link << '>' << aliasType.name << '^' << aliasType.templateTypes.listText << '\n';
    if (!aliasName.empty())
        sink.link << '&' << replaceSet(TmpAlias,table) << '\n';
    sink.link << block.linebreaks();
}
