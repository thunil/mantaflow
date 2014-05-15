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
static PyObject* _P_$FUNCNAME (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
    try {
        PbArgs _args(_linargs, _kwds);
        FluidSolver *parent = _args.obtainParent();
        pbPreparePlugin(parent, "$FUNCNAME" );
        PyObject *_retval = 0;
        { 
            ArgLocker _lock;
            $ARGLOADER
            @IF($RET_VOID)
                _retval = getPyNone();
                $FUNCNAME($CALLSTRING);
            @ELSE
                _retval = d_toPy($FUNCNAME($CALLSTRING));
            @END
            _args.check();
        }
        pbFinalizePlugin(parent,"$FUNCNAME" );
        return _retval;
    } catch(std::exception& e) {
        pbSetError("$FUNCNAME",e.what());
        return 0;
    }
}
static const PbRegFunc _RP_$FUNCNAME ("","$FUNCNAME",_P_$FUNCNAME);
);

const string TmpMemberFunction = STR(
static PyObject* _$FUNCNAME (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
    try {
        PbArgs _args(_linargs, _kwds);
        $CLASS* pbo = dynamic_cast<$CLASS*>(PbClass::fromPyObject(_self));
        pbPreparePlugin(pbo->getParent(), "$CLASS::$FUNCNAME");
        PyObject *_retval = 0;
        { 
            ArgLocker _lock;
            $ARGLOADER
            pbo->_args.copy(_args);
            @IF($RET_VOID)
                _retval = getPyNone();
                pbo->$FUNCNAME($CALLSTRING);
            @ELSE
                _retval = d_toPy(pbo->$FUNCNAME($CALLSTRING));
            @END
            _args.check();
        }
        pbFinalizePlugin(pbo->getParent(),"$CLASS::$FUNCNAME");
        return _retval;
    } catch(std::exception& e) {
        pbSetError("$CLASS::$FUNCNAME",e.what());
        return 0;
    }
});

const string TmpConstructor = STR(
static int _$CLASS (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {
    PbClass* obj = PbClass::fromPyObject(_self); 
    if (obj) delete obj; 
    try {
        PbArgs _args(_linargs, _kwds);
        pbPreparePlugin(0, "$CLASS::$FUNCNAME" );
        { 
            ArgLocker _lock;
            $ARGLOADER
            obj = new $CLASS($CALLSTRING);
            std::string _name = _args.getOpt<std::string>("name",""); 
            obj->setPyObject(_self); 
            if (!_name.empty()) obj->setName(_name); 
            _args.check();
        }
        pbFinalizePlugin(obj->getParent(),"$CLASS::$FUNCNAME" );
        return 0;
    } catch(std::exception& e) {
        pbSetError("$CLASS::$FUNCNAME",e.what());
        return -1;
    }
});

const string TmpRegisterMethod = STR(
@IF($CLASSTPL)
    static const PbRegFunc _R_$CLASS_$CLASSTPL_$FUNCNAME ("$CLASS<$CLASSTPL>$","$FUNCNAME",$CLASS<$CLASSTPL>::_$FUNCNAME);
@ELSE
    static const PbRegFunc _R_$CLASS_$FUNCNAME ("$CLASS","$FUNCNAME",$CLASS::_$FUNCNAME);
@END
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

string createConverters(const string& name, const string& tb, const string& nl, const string& nlr) {
    return "template<> " + name + "* fromPy<" + name + "*>(PyObject* obj) {" + nl +
           tb+ "if (PbClass::isNullRef(obj)) return 0;" + nl +
           tb+ "PbClass* pbo = PbClass::fromPyObject(obj);" + nl +
           tb+ "if (!pbo || !(pbo->canConvertTo(\"" + name + "\")))" + nl +
           tb+tb+ "throw Error(\"can't convert argument to type '" + name + "'\");" + nl +
           tb+ "return dynamic_cast<" + name + "*>(pbo);" + nl +
           "}" + nlr +
           "template<> PyObject* toPy< " + name + " >( " + name + "& v) {" + nl +
           tb+ "if (v.getPyObject()) return v.getPyObject();" + nl +
           tb+ name + "* co = new " + name + " (v); " +
           tb+ "return co->assignNewPyObject(\""+name+"\");" + nl +
           "}" + nlr;
}

// globals for tracking state between python class and python function registrations
string gLocalReg, gParent;
bool gFoundConstructor = false, gIsTemplated=false;

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
                             "CLASSTPL", (isPlugin || !block.parent->isTemplated()) ? "" : "$CT",
                             "CALLSTRING", func.callString(),
                             "RET_VOID", (func.returnType.name=="void") ? "Y" : "",
                             "@end" };
    string callerTempl = isConstructor ? TmpConstructor : (isPlugin ? TmpFunction : TmpMemberFunction);
    sink.inplace << replaceSet(callerTempl, table);

    // register member functions
    if (!isPlugin) {
        const string reg = replaceSet(TmpRegisterMethod, table);
        if (block.parent->isTemplated())
            sink.tplReg.push_back(RegCall(block.parent->name,reg));
        else 
            sink.ext << reg << "\n";
    }
}

string processPythonVariable(const Block& block, Sink& sink) {
    // beautify code
    string nl = gDebugMode ? "\n" : "";
    string tb = gDebugMode ? "\t" : "";
    const int line = block.line0;
    
    if (gParent.empty())
        errMsg(line, "PYTHON variables con only be used inside classes");
    
    string name=block.func.name;
    string pname=name, type = block.func.returnType.name;
    
    // process options
    for (size_t i=0; i<block.options.size(); i++) {
        if (block.options[i].name == "name") 
            pname = block.options[i].value;
        else
            errMsg(line, "PYTHON(opt): illegal option. Supported options are: 'name'");
    }
    
    // define getter / setter
    string nn = gParent + "_" +name;
    string gethdr = "PyObject* _get_" + nn + "(PyObject* self, void* cl)";
    string sethdr = "int _set_" + nn + "(PyObject* self, PyObject* val, void* cl)";
    
    // replicate original code plus accessor
    string code = "";
    code += type + " " + name + ";" + nl;
    code += tb + "friend " + gethdr + ";" + nl;
    code += tb + "friend " + sethdr + ";" + nl;
    
    // add get/setter
    string getter = gethdr+" { return d_toPy(fromPy<" + gParent+"*>(self)->" + name + "); }";
    string setter = sethdr+" { fromPy<" + gParent+"*>(self)->" + name + "=fromPy<" + type + " >(val); return 0;}";
        
    // register
    gRegText += "PbWrapperRegistry::instance().addGetSet(\""+gParent+"\",\""+pname+"\",_get_"+nn+ ",_set_"+nn+");\n";
    if (gIsHeader) {
        gRegText += getter+"\n"+setter+"\n";
    } else {
        gLocalReg += getter + nl + setter + nl;
        gRegText += "extern " + gethdr + ";\nextern " + sethdr + ";\n";
    }
    
    return block.linebreaks() + code;
}

void processPythonClass(const Block& block, const string& code, Sink& sink) {
    const Class& cur = block.cls;
    


    // beautify code
    string nl = gDebugMode ? "\n" : "";
    string tb = gDebugMode ? "\t" : "";
    const List<Type>& templArgs = block.cls.templateTypes;
    const int line = block.line0;

    if (!gIsHeader && !templArgs.empty())
        errMsg(line, "PYTHON template classes can only be defined in header files.");
    
    string name = block.cls.name;
    string pname = name; // for now
    string baseclassName = block.cls.baseClass.name;
    const List<Type>& baseclassTempl = block.cls.baseClass.templateTypes;
    
    // process options
    for (size_t i=0; i<block.options.size(); i++) {
        if (block.options[i].name == "name") 
            pname = block.options[i].value;
        else
            errMsg(line, "PYTHON(opt): illegal kernel option. Supported options are: 'name'");
    }
    
    // class registry
    string baseclass = baseclassName, modBase = baseclassName, registry = "", implInst = "";
    if (!baseclassTempl.empty()) {
        // baseclass known ? try to implicitly instantiate base class
        string targ="", tcarg="" ,bclist="";
        bool chain=false;
        for (size_t i=0; i<baseclassTempl.size(); i++) {
            // check if template arg            
            int index = -1;
            for (size_t j=0; j<templArgs.size(); j++) {
                if (templArgs[j].name == baseclassTempl[i].name) {
                    index = j;
                    chain=true;
                    break;
                }
            }
            if (index>=0) {
                targ += "@"+baseclassTempl[i].name+"@";
                stringstream s;
                s << "$" << index << "$";
                bclist += s.str();
            } else {
                targ += baseclassTempl[i].name;
                bclist += baseclassTempl[i].name;
            }
            if (i!=baseclassTempl.size()-1) { targ += ","; bclist += ","; }
        }
        for (size_t i=0; i<templArgs.size(); i++) {
            tcarg += "@" + templArgs[i].name + "@";
            if (i!=templArgs.size()-1) tcarg += ",";
        }
        string aliasn = "_" + baseclassName + "_" + targ;
        replaceAll(aliasn,",","_");
        
        // need defer chain ?
        if (chain){
            gRegText += "@chain " + name + " " + tcarg + " " + baseclassName + " " + targ + "\n";
        } else {
            gRegText += "@instance " + baseclassName + " " + targ + " " + aliasn + "\n";    
        }
        modBase += "<" + bclist + ">";
        baseclass += baseclassTempl.minimal;
    }
    if (templArgs.empty())
        registry = "PbWrapperRegistry::instance().addClass(\"" + pname + "\", \"" + name + "\", \"" + modBase + "\");";
    else {
        registry = "@template " + name + " PbWrapperRegistry::instance().addClass(\"@\", \"" + name + "<$$>\", \"" + modBase + "\");";
    }
    
    // register class
    if (gIsHeader) {
        // make sure we
        string fn = gFilename;
        const size_t p = fn.find_last_of('/');
        if (p != string::npos) fn=fn.substr(p+1);
        gRegText += "#include \"" + fn + "\"\n";        
    }
    gRegText += registry + "\n";
    
    // register converters
    gLocalReg = ""; 
    if (templArgs.empty()) {
        if (gIsHeader)
            gRegText += createConverters(name, "", " ", "\n");
        else
            gLocalReg += createConverters(name, tb, nl, nl);
    }
    
    // tokenize and parse contained python functions
    gParent = name;
    gIsTemplated = !templArgs.empty();
    gFoundConstructor = false;

    // create class
    string pclass = "";
    if (gDocMode) {
        pclass += "//! \\ingroup PyClasses\nPYTHON ";
    }
    if (gIsTemplated)
        pclass += "template" + templArgs.minimal + nl;
    pclass += "class " + name + " : public " + baseclass + " {" + nl;
    
    sink.inplace << (block.linebreaks() + implInst);
    sink.inplace << pclass;
    processText(code.substr(1), line, sink, &cur);
    gParent = "";    
    if (!gFoundConstructor)
        errMsg(line, "no PYTHON constructor found in class '" + name + "'");
    if (!gIsHeader && gIsTemplated)
        errMsg(line, "PYTHON class template can only be used in header files.");
    
    pclass = nl;
    if (!gDocMode) {
        pclass += "public:" + nl;
        pclass += tb+"PbArgs _args;" + nl;
    }
    pclass += "};" + nl;
    
    sink.inplace << pclass << gLocalReg;
}

set<string> gAliasRegister;
string processPythonInstantiation(const Block& block, const Type& aliasType, const string& aliasName) {
    gRegText += "@instance " + aliasType.name + " " + aliasType.templateTypes.listText + " " + aliasName + "\n";
    
    if (gAliasRegister.find(aliasName) == gAliasRegister.end()) {
        gAliasRegister.insert(aliasName);
        return block.linebreaks() + "typedef " + aliasType.minimal + " " + aliasName + "; ";
    }
    return block.linebreaks();
}
