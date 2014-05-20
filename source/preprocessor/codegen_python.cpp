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
$TEMPLATE$ static PyObject* _$WRAPPER$ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
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
);

const string TmpMemberFunction = STR(
$TEMPLATE$ static PyObject* _$WRAPPER$ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
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
$TEMPLATE$ static int _$WRAPPER$ (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
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
@ELIF(CLASS)
    static const Pb::Register _R_$IDX$ ("$CLASS$","$FUNCNAME$",$CLASS$::_$FUNCNAME$);
@ELSE
    static const Pb::Register _RP_$FUNCNAME$ ("","$FUNCNAME$",_$FUNCNAME$);
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

const string TmpTemplateWrapper = STR(
static $RET$ _$FUNCNAME$ (PyObject* s, PyObject* l, PyObject* kw) {
    PbArgs args(l, kw);
    int hits=0;
    $RET$ (*call)(PyObject*,PyObject*,PyObject*);

    $TEMPLATE_CHECK$
    
    if (hits == 1)
        return call(s,l,kw);
    if (hits == 0)
        pbSetError("$FUNCNAME$", "Can't deduce template parameters");
    else
        pbSetError("$FUNCNAME$", "Argument matches multiple templates");
    return @IF(CONSTRUCTOR) -1 @ELSE 0 @END ; 
}
);

const string TmpTemplateChecker = STR(
template $TEMPLATE$
static bool $NAME$ (PbArgs& A) {
    return $CHK$;
}
);

//******************************************************
// Code generation functions

string generateLoader(const Argument& arg) {
    bool integral = isIntegral(arg.type.name);
    Type ptrType = arg.type;
    string optCall =  arg.value.empty() ? "get" : "getOpt";

    ptrType.isConst = false;
    if (integral) {
        ptrType.isPointer = false;
        ptrType.isRef = false;
    } else if (arg.type.isPointer) {
        ptrType.isPointer = false;
        ptrType.isRef = false;
        optCall = arg.value.empty() ? "getPtr" : "getPtrOpt";
    }
    
    stringstream loader;
    loader << arg.type.build() << " " << arg.name << " = ";
    loader << "_args." << optCall;
        
    loader << "<" << ptrType.build() << " >";
    loader << "(\"" << arg.name << "\"," << arg.index << ",";
    if (!arg.value.empty())
        loader << arg.value << ",";
    loader << "&_lock); ";

    return loader.str();
}

// global for tracking state between python class and python function registrations
bool gFoundConstructor = false;

void processPythonFunction(const Block& block, const string& code, Sink& sink, vector<Instantiation>& inst) {
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
    for (int i=0; i<(int)func.arguments.size(); i++)
        loader += generateLoader(func.arguments[i]);

    // generate glue layer function
    const string table[] = { "FUNCNAME", func.name, 
                             "ARGLOADER", loader, 
                             "TEMPLATE", func.isTemplated() ? "template "+func.templateTypes.minimal : "",
                             "WRAPPER", (func.isTemplated() ? "T_":"") + func.name,
                             "CLASS", isPlugin ? "" : block.parent->name,
                             "CTPL", (isPlugin || !block.parent->isTemplated()) ? "" : "$CT$",
                             "CALLSTRING", func.callString(),
                             "RET_VOID", (func.returnType.name=="void") ? "Y" : "",
                             "@end" };
    string callerTempl = isConstructor ? TmpConstructor : (isPlugin ? TmpFunction : TmpMemberFunction);
    sink.inplace << replaceSet(callerTempl, table);


    // drop a marker for function template wrapper
    if (func.isTemplated()) {
        stringstream num; num << inst.size();
        sink.inplace << "$" << num.str() << "$";
        inst.push_back(Instantiation(isPlugin ? "" : block.parent->name, func));
    }

    // register functions
    const string reg = replaceSet(TmpRegisterMethod, table);
    if (isPlugin) {
        sink.inplace << reg;
    } else {
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

void processPythonClass(const Block& block, const string& code, Sink& sink, vector<Instantiation>& inst) {
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
    processText(ncode.substr(1), block.line0, sink, &cls, inst);
    if (!gFoundConstructor)
        errMsg(block.line0, "no PYTHON constructor found in class '" + cls.name + "'");
    
    // add secret bonus members to class and close
    sink.inplace << "public: PbArgs _args;";
    sink.inplace << "}\n";
    // add a define to make commenting out classes, and #ifdefs work correctly
    sink.inplace << "#define _C_" << cls.name << '\n';
}

void processPythonInstantiation(const Block& block, const Type& aliasType, Sink& sink, vector<Instantiation>& inst) {
    string parent = block.parent ? block.parent->name : "";
    // for template functions, add to instantiation list
    bool isFunction = false;
    for (int i=0; i<(int)inst.size(); i++) {
        if (inst[i].cls == parent && inst[i].func.name == aliasType.name) {
            inst[i].templates.push_back(aliasType.templateTypes.listText);
            isFunction = true;
            break;
        }
    }
    // otherwise, assume it's a class, and put down a link-time instantiation request
    if (!isFunction) {
        sink.link << '>' << aliasType.name << '^' << aliasType.templateTypes.listText << '\n';    
    }
}

void processPythonAlias(const Block& block, const Type& aliasType, const string& aliasName, Sink& sink) {
    const string table[] = {"CLASS", strip(aliasType.build()), "PYNAME", aliasName, "@end"};
    if (!aliasName.empty())
        sink.link << '&' << replaceSet(TmpAlias,table) << '\n';
}

// build the template argument checker needed for template deduction in the wrapper
string buildTemplateChecker(string& out, const Function& func) {
    stringstream chk;
    for (int k=0; k<func.arguments.size(); k++) {
        stringstream num; num << k;    
        Type type = func.arguments[k].type;
        type.isPointer = false;
        type.isRef = false;
        type.isConst = false;
        chk << "A.typeCheck<" << type.build() << " >(" 
            << num.str() << ",\"" << func.arguments[k].name << "\")";

        if (k != func.arguments.size()-1)
            chk << " && ";
    }

    static int idx = 0;
    stringstream name; 
    name << "_K_" << idx++;
    const string table[] = { "TEMPLATE", func.templateTypes.minimal,
                             "NAME", name.str(),
                             "CHK", chk.str(),
                             "@end" };
    out+= replaceSet(TmpTemplateChecker,table);
    return name.str();
}

// add a wrapper for all templated function 
void postProcessInstantiations(Sink& sink, vector<Instantiation>& inst) {
    string out = sink.inplace.str();
    for (int i=0; i<(int)inst.size(); i++) {
        Instantiation& cur = inst[i];
        if (cur.templates.size() == 0)
            errMsg(0, cur.cls + "::" + cur.func.name + " : templated function without instantiation detected.");

        string wrapper = "";
        string chkFunc = buildTemplateChecker(wrapper, cur.func);
        stringstream chkCall;

        // build argument checker        
        for (int j=0; j<cur.templates.size(); j++) {
            stringstream num; num << j;
            chkCall << "if (" << chkFunc << "<" << cur.templates[j] << ">(args)) {";
            chkCall << "hits++; call = _T_" << cur.func.name << "<" << cur.templates[j] <<">; }";
        }
        const string table[] = { "CONSTRUCTOR", cur.func.name == cur.cls ? "Y":"",
                                 "FUNCNAME", cur.func.name,
                                 "TEMPLATE_CHECK", chkCall.str(),
                                 "RET", cur.func.name == cur.cls ? "int" : "PyObject*",
                                 "@end" };
        wrapper += replaceSet(TmpTemplateWrapper, table);
        
        stringstream num; 
        num << "$" << i << "$";
        replaceAll(out, num.str(), wrapper);
    }
    sink.inplace.str(out);
}