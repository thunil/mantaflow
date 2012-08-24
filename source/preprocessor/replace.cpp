/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor: Process replacement text of individual keywords
 *
 ******************************************************************************/

#include "prep.h"
#include <cstdlib>
#include <set>
#include <sstream>
#include <iostream>
using namespace std;

string buildline(int lb) {
    string lbs="";
    for (int i=0; i<lb; i++) lbs+="\n";
    return lbs;
}

enum FunctionType { FtPlugin, FtMember, FtConstructor };
void createPythonWrapper(const vector<Argument>& args, const string& fname, const string& totalname, FunctionType ftype, string& header, string& footer, string& callList, bool isClass) {
    // beautify code for debug mode
    string nl = (gDebugMode && ftype != FtConstructor) ? "\n" : " ";
    string tb = (gDebugMode && ftype != FtConstructor) ? (isClass ? "\t\t" : "\t") : "";
    string tb1 = (gDebugMode && ftype != FtConstructor) ? (isClass ? "\t" : "") : "";
    string tb2 = tb1+tb;
    
    const string argstr = (ftype == FtMember) ? "__args" : "_args";
    // load arguments from python
    callList = "";
    string loader = "", argList = "";
    for (size_t i=0; i<args.size(); i++) {
        string type = args[i].type;        
        bool itype = isIntegral(args[i].type);
        string name = args[i].name;
        stringstream num; num << args[i].number;
        if (!args[i].templ.empty()) type+= "<" + args[i].templ + ">";
        if (args[i].isPointer) type += "*";        
        string completeType = (args[i].isConst) ? ("const " + type) : type;
                
        if (args[i].isRef) {
            if (args[i].value.empty()) {                
                // for strings etc. direct get, for PbClass use pointers
                if (itype)
                    loader += tb2+ completeType + "& " + name + " = " + argstr + ".get< " + type + " > (" + num.str() + ",\"" + name + "\", &_lock);" + nl;
                else
                    loader += tb2+ completeType + "& " + name + " = *" + argstr + ".get< " + type + "* > (" + num.str() + ",\"" + name + "\", &_lock);" + nl;
            } else {                
                // for strings etc. direct get, for PbClass use pointers
                if (itype) {
                    loader += tb2+ completeType + "& " + name + " = " + argstr + ".getOpt< " + type + " > (" + num.str() + ",\"" + name + "\", " + args[i].value + ", &_lock);" + nl;
                } else {
                    loader += tb2+ completeType + "* _ptr_" + name + " = " + argstr + ".getOpt< " + type + "* > (" + num.str() + ",\"" + name + "\", 0, &_lock);" + nl;
                    loader += tb2+ completeType + "& " + name + " = (_ptr_" + name + ") ? (*_ptr_" + name + ") : (" + args[i].value + ");" + nl; 
                }
            }
            type += "&";
        }
        else {
            loader += tb2+ completeType + " " + name + " = ";
            if (args[i].value.empty()) {
                loader += argstr + ".get< " + type + " > (" + num.str() + ",\"" + name + "\", &_lock);" + nl;
            } else {
                loader += argstr + ".getOpt< " + type + " > (" + num.str() + ",\"" + name + "\", " + args[i].value + ", &_lock);" + nl;
            }
        }
        if (i!=0) callList +=", ";
        callList += name;
    }
    
    // generate header & footer
    if (ftype == FtConstructor) {
        header = "int " + fname + " (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {" + nl;
        header += tb2+ "PbClass* obj = PbClass::fromPyObject(_self);" + nl;
        header += tb2+ "if (obj) delete obj;" + nl;                
    } else
        header = "PyObject* " + fname + " (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {" + nl;
    header += tb+"try {" + nl;
    header += tb2+ "PbArgs " + argstr + "(_linargs, _kwds);" + nl;
    if (ftype == FtMember)
        header += tb2+ "PyObject *_retval = NULL;" + nl;
    if (ftype == FtPlugin) {
        header += tb2+ "FluidSolver *parent = _args.obtainParent();" + nl;
        header += tb2+ "pbPreparePlugin(parent, \""+totalname+"\");" + nl;
        header += tb2+ "PyObject *_retval = NULL;" + nl;
    } else if (ftype == FtMember)
        header += tb2+ "pbPreparePlugin(this->mParent, \""+totalname+"\");" + nl;        
    else 
        header += tb2+ "pbPreparePlugin(0, \""+totalname+"\");" + nl;    
    header += tb2+ "{ ArgLocker _lock;" + nl;
    header += loader;
    if (ftype == FtMember)
        header += tb2+ "this->_args.copy(__args);" + nl;
            
    if (ftype == FtConstructor) {
        footer =  tb2+ "std::string _name = _args.getOpt<std::string>(\"name\",\"\");" + nl;
        footer += tb2+ "obj->setPyObject(_self);" + nl;
        footer += tb2+ "if (!_name.empty()) obj->setName(_name);" + nl;
        footer += tb2+ "_args.check(); }" + nl;
        footer += tb2+ "pbFinalizePlugin(obj->getParent(),\"" +totalname+"\");" + nl;
        footer += tb2+ "return 0;" + nl;
    } else if (ftype == FtMember) {
        footer =  tb2+ "this->_args.check(); }" + nl;
        footer += tb2+ "pbFinalizePlugin(this->mParent,\"" +totalname+"\");" + nl;
        footer += tb2+ "return _retval;" + nl;    
    } else if (ftype == FtPlugin) {
        footer =  tb2+ "_args.check(); }" + nl;
        footer += tb2+ "pbFinalizePlugin(parent,\"" +totalname+"\");" + nl;
        footer += tb2+ "return (_retval) ? _retval : getPyNone();" + nl;    
    }
    footer += tb+ "} catch(std::exception& e) {" + nl;
    footer += tb2+ "pbSetError(\"" + totalname + "\",e.what());" + nl;
    footer += tb2+ "return " + ((ftype==FtConstructor) ? "-1;" : "0;") + nl;
    footer += tb+ "}" + nl;
    footer += tb1+ "}" + nl;
}

string createConverters(const string& name, const string& tb, const string& nl, const string& nlr) {
    return "template<> " + name + "* fromPy<" + name + "*>(PyObject* obj) {" + nl +
           tb+ "if (PbClass::isNullRef(obj)) return 0;" + nl +
           tb+ "PbClass* pbo = PbClass::fromPyObject(obj);" + nl +
           tb+ "if (!pbo || !(pbo->canConvertTo(\"" + name + "\")))" + nl +
           tb+tb+ "throw Error(\"can't convert argument to type '" + name + "'\");" + nl +
           tb+ "return dynamic_cast<" + name + "*>(pbo);" + nl +
           "}" + nlr +
           "template<> PyObject* toPy<" + name + "*>(" + name + "* v) {" + nl +
           tb+ "return v->getPyObject();" + nl +
           "}" + nlr;
}

string processKernel(int lb, const string& kname, const vector<Argument>& opts, const string& templ, const vector<Argument>& args, const string& codebase, int line, bool isClass) {    
    // beautify code
    string nl = gDebugMode ? "\n" : "";
    string tb = gDebugMode ? "\t" : "";    
    string tb2 = tb+tb, tb3=tb2+tb, tb4=tb3+tb, tb5=tb4+tb;
    
    if (gDocMode) {
        string ds = "//! \\ingroup Kernels\nKERNEL<";
        for(size_t i=0; i<opts.size(); i++) { if (i!=0) ds+=", "; ds+=opts[i].complete; }
        ds += "> " + kname + " (";
        for(size_t i=0; i<args.size(); i++) { if (i!=0) ds+=", "; ds+=args[i].complete;}
        return ds + " ) {}\n";
    }
    
    // process options
    bool idxMode = false, reduce = false, mt = gUseMT, pts = false;
    string bnd = "0";
    for (size_t i=0; i<opts.size(); i++) {
        if (opts[i].name == "ijk") 
            idxMode = false;
        else if (opts[i].name == "index" || opts[i].name == "idx")
            idxMode = true;
        else if (opts[i].name == "st" || opts[i].name == "single")
            mt = false;
        else if (opts[i].name == "pts" || opts[i].name == "particle" || opts[i].name == "points")
            pts = true;
        else if (opts[i].name == "bnd")
            bnd = opts[i].value;
        else if (opts[i].name == "bnd")
            bnd = opts[i].value;
        else if (opts[i].name == "reduce")
            reduce = true;
        else
            errMsg(line, "KERNEL(opt): illegal kernel option. Supported options are: 'ijk', 'idx', 'bnd=x', 'reduce', 'st', 'pts'");
    }
    string qualifier = reduce ? "" : "const";
    string tbbcall = reduce ? "tbb::parallel_reduce" : "tbb::parallel_for";

	bool mtTbb    = false;
	bool mtOpenMp = false; // testing... NT_DEBUG
	if(mt) {
		mtTbb = true;
	}
    
    // point out illegal paramter combinations
    if (bnd != "0" && idxMode)
        errMsg(line, "KERNEL(opt): can't combine index mode with bounds iteration.");    
    if (reduce && !isClass)
        errMsg(line, "KERNEL(opt): reduce can only be used in kernel classes.");
    if (pts && (idxMode || bnd != "0" ))
        errMsg(line, "KERNEL(opt): Modes 'ijk', 'idx' and 'bnd' can't be applied to particle kernels.");

    // depnding on loop/function, adapt return statement
    string code = codebase;
    if (!isClass) {
        replaceAll(code, "return;", "continue;");
        replaceAll(code, "return ", "continue ");
    }

    // parse arguments
    string loader = "", initList = "", argList = "", copier = "", members = "", basegrid="";
    for (size_t i=0; i<args.size(); i++) {
        string type = args[i].type;
        string name = args[i].name;
        if (!args[i].templ.empty()) type+= "<" + args[i].templ + ">";                
        if (args[i].isPointer) type += "*";                
        if (args[i].isRef) type += "&";
        if (args[i].isConst) type = "const "+type;
        
        initList += (i==0) ? "" : ", ";
        argList += (i==0) ? "" : ", ";
        copier += (i==0) ? "" : ", ";
        
        string aname = "_" + name;
        string sname = (isClass ? "" : "m_" )+ name;
        loader += tb2+ type + " " + name + " = " + sname + ";" + nl;
        copier += sname + "(o." + sname + ")";
        members += tb+ type + " " + sname + ";" + nl;
        initList += sname + "(" + aname + ")";
        argList += type + " " + aname;
        if (!args[i].value.empty())
            argList += " = " + args[i].value;
        
        // figure out basegrid
        if (basegrid.empty()) {
            if (!pts && !type.find("Grid") != string::npos) {
                if (!args[i].isPointer) basegrid += "&";
                basegrid += "_" + name;
            } else if (pts) {
                if (args[i].isPointer) basegrid += "*";
                basegrid += "_" + name;
            }
        }
    }
    if (basegrid.empty()) 
        errMsg(line, "KERNEL: use at least one grid to call the kernel.");
    if (isClass)
        loader = "";
    
    // create kernel class
    string kclass = "", kclassname = kname;
    if (!templ.empty()) kclass += "template <" + templ + ">" + nl;
    kclass += "struct " + kclassname + " : public " + (pts ? "Particle" : "") + "KernelBase { " + nl;
    kclass += "public:" + nl;
    // init constructor
    kclass += tb+ kclassname + " (" + argList + ") :" + nl;
    if (pts)
        kclass += tb2+ "ParticleKernelBase((" + basegrid + ").size()), " + initList + nl;
    else
        kclass += tb2+ "KernelBase(" + basegrid + ", " + bnd + "), " + initList + nl;
    kclass += tb + "{" + nl;
    if (reduce) kclass += tb2+ "setup();" + nl;
    kclass += tb2+ "run();" + nl;
    kclass += tb+ "}" + nl;    
    if (mtTbb) {
        // multithreading using intel tbb
        kclass += tb+ "void operator() (const tbb::blocked_range<size_t>& r) " + qualifier + "{" + nl;
        kclass += loader;
        if (pts) {
            kclass += tb2+ "for (int i=r.begin(); i!=(int)r.end(); i++)" + nl;
            kclass += isClass ? (tb3+"(*this)(i);") : code;                    
        } else if (idxMode) {
            kclass += tb2+ "for (int idx=r.begin(); idx!=(int)r.end(); idx++)" + nl;
            kclass += isClass ? (tb3+"(*this)(idx);") : code;
        } else {
            kclass += tb2+ "const int _maxX = maxX, _maxY=maxY;" + nl;            
            kclass += tb2+ "for (int k=r.begin(); k!=(int)r.end(); k++)" + nl;
            kclass += tb3+ "for (int j=" + bnd + "; j<_maxY; j++)" + nl;
            kclass += tb4+ "for (int i=" + bnd + "; i<_maxX; i++)" + nl;
            kclass += isClass ? (tb5+"(*this)(i,j,k);") : code;
        }            
        kclass += nl + tb+"}" + nl;
        kclass += tb+ "void run() {" + nl;
        kclass += tb2+ tbbcall + "(tbb::blocked_range<size_t>("+bnd+", " + (pts ? "size" : (idxMode ? "maxCells" : "maxZ")) + "), *this);"+ nl;
        kclass += tb+ "}" + nl;
	} else if(mtOpenMp) {
	} else {
        // simple loop
        kclass += tb+ "void run() {" + nl;
        kclass += loader;
        if (pts) {
            kclass += tb2+ "const int _sz = size;"+ nl;
            kclass += tb2+ "for (int i=0; i < _sz; i++)" + nl;
            kclass += isClass ? (tb3+"(*this)(i);") : code;
        } else if (idxMode) {
            kclass += tb2+ "const int _maxCells = maxCells;"+ nl;
            kclass += tb2+ "for (int idx=0; idx < _maxCells; idx++)" + nl;
            kclass += isClass ? (tb3+"(*this)(idx);") : code;
        } else {
            kclass += tb2+ "const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ;" + nl;
            kclass += tb2+ "for (int k="+bnd+"; k < _maxZ; k++)" + nl;
            kclass += tb3+ "for (int j="+bnd+"; j < _maxY; j++)" + nl;
            kclass += tb4+ "for (int i="+bnd+"; i < _maxX; i++)" + nl;
            kclass += isClass ? (tb5+"(*this)(i,j,k);") : code;
        }
        kclass += nl + tb+"}" + nl;
    }
    if (reduce && mtTbb) {
        // split constructor
        kclass += tb+ kclassname + " (" + kclassname + "& o, tbb::split) : " + nl;
        if (!pts)
            kclass += tb2+ "KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.X, o.Y, o.Z)," + nl;
        else
            kclass += tb2+ "KernelBase(o.size)," + nl;
        kclass += tb2+ copier + nl;
        kclass += tb + "{" + nl; 
        kclass += tb2+ "setup();" + nl;
        kclass += tb+ "}" + nl;        
    }    
    kclass += nl + members;
    kclass += isClass ? code.substr(1) : "}";
    kclass += ";" + nl;
   
    return buildline(lb) + kclass;
}

string processPlugin(int lb, const string& fname, const vector<Argument>& opts, const vector<Argument>& args, const string& code, int line, const string& type) {    
    // beautify code for debug mode
    string nl = gDebugMode ? "\n" : " ";
    string tb1 = gDebugMode ? "\t" : "";
    bool isClass = (type=="class") || (type=="struct");
    
    if (!opts.empty())
        errMsg(line,"Keyword PLUGIN does not support options");
    
    if (gDocMode) {
        string ds = "//! \\ingroup Plugins\nPLUGIN " + fname + "( ";
        for(size_t i=0; i<args.size(); i++) {  if (i!=0) ds+=", "; ds+=args[i].complete;}
        return ds + " ) {}\n";
    }
    
    string wrapperName = "plugin_" + fname;
    
    // register function
    gRegText += "extern PyObject* " + wrapperName + " (PyObject* _self, PyObject* _linargs, PyObject* _kwds);\n" + // declare
                "PbWrapperRegistry::instance().addGenericFunction(\"" + fname + "\", " + wrapperName + ");\n"; // register                    
    
    // generate wrapper
    string codeInline = "";
    codeInline += code.substr(1,code.size()-2);    
    string header, footer, callList;
    createPythonWrapper(args, isClass ? "step" : wrapperName, fname, FtPlugin, header, footer, callList, isClass);
        
    if (isClass) {        
        string caller = "", classcode="";
        // create caller
        caller += "PyObject* " + wrapperName + " (PyObject* _self, PyObject* _linargs, PyObject* _kwds) {" + nl;
        caller += tb1+"static " + fname + " inst;" + nl;
        caller += tb1+"return inst.step(_self, _linargs, _kwds);" + nl;
        caller += "}" + nl;
        
        // create class
        classcode += "class " + fname + " {" + nl;
        classcode += "public:" + nl;
        classcode += replaceFunctionHeader(codeInline, "step", header, footer, line) + nl;
        classcode += "};" + nl;        
        return buildline(lb) + classcode + nl + caller;
    }
    else {
        // replicate original function
        string callstr = "";
        string func = type + " " + fname + "(";
        for (size_t i=0; i<args.size(); i++) {
            if (i!=0) { func += ", "; callstr += ","; }
            func += args[i].complete;
            callstr += args[i].name;
        }
        // add args, parent arguments
        if (args.size()>0) func += ",";
        func += " FluidSolver* parent = NULL, PbArgs& _args = PbArgs::EMPTY";
        func += ") {" + nl + codeInline + nl + "}" + nl;
        if (args.size()>0) callstr += ",";
        callstr += " parent, _args";
        
        string caller = header;
        if (type=="void")
            caller += tb1+"_retval = getPyNone(); " + fname + "(" + callstr + ");" + nl;
        else
            caller += tb1+"_retval = toPy<" + type + ">(" + fname + "(" + callstr + ") );" + nl;
        caller += footer;
        
        // inject code directly into caller function        
        return buildline(lb) + func + caller;
    }
}

// globals for tracking state between python class and python function registrations
string gLocalReg, gParent;
bool gFoundConstructor = false, gIsTemplated=false;

string processPythonFunction(int lb, const string& name, const string& type, const vector<Argument>& args, const string& initlist, const string& code, int) {
    // beautify code
    string nl = gDebugMode ? "\n" : "";
    string tb = gDebugMode ? "\t" : "";
    
    // is header file ?
    bool isHeader = gFilename[gFilename.size()-2] == '.' && gFilename[gFilename.size()-1] == 'h';
    bool isConstructor = type.empty();
    if (isConstructor) gFoundConstructor = true;
        
    // generate caller
    string clname = "_" + gParent + (gIsTemplated ? "@" : "");
    string fname = "_" + name;
    string codeInline = code;
    if (code[0] != ';') {
        codeInline = "{" + nl;
        codeInline += code.substr(1,code.size()-1) + nl;        
    }
        
    string header, footer, callList;
    createPythonWrapper(args, isConstructor ? (clname+clname) : fname, gParent+"::"+name, isConstructor ? FtConstructor : FtMember, header, footer, callList, true);    
    
    string caller = tb + header + nl;
    if (type == "void") {
        caller += tb+tb+tb+ "_retval = getPyNone();" + nl;
        caller += tb+tb+tb+ name + "(" + callList + ");" + nl + nl;
    } else
        caller += tb+tb+tb+ "_retval = toPy(" + name + "(" + callList + ") );" + nl + nl;
    caller += footer;
    
    // replicate original function
    string func = type + " " + name + "(";
    for (size_t i=0; i<args.size(); i++) {
        if (i!=0) func += ", ";
        func += args[i].complete;
    }
    func += ") " + initlist + codeInline + nl;
    
    // register
    string regname = clname + (isConstructor ? clname : fname);
    string regHeader = (isConstructor ? "int " : "PyObject* ") + regname + " (PyObject* _self, PyObject* _linargs, PyObject* _kwds)";
    string regDecl = "", regCall = "";
    if (isConstructor) {
        caller = "";
        if (gIsTemplated) {
            regDecl = "@template " + gParent + " " + header + "obj = new " + gParent + "<$> (" + callList + ");" + footer;
            regCall = "@template " + gParent + " PbWrapperRegistry::instance().addConstructor(\"" + gParent + "<$>\", " + regname + ");";        
        } else {
            regDecl = header + "obj = new " + gParent + "(" + callList + ");" + footer;
            regCall = "PbWrapperRegistry::instance().addConstructor(\"" + gParent + "\", " + regname + ");";
        }
    } else {
        if (gIsTemplated) {
            regDecl = "@template " + gParent + " " + regHeader + " { return dynamic_cast<" + gParent +"<$>*>(PbClass::fromPyObject(_self))->_" + name + "(_self, _linargs, _kwds); }";
            regCall = "@template " + gParent + " PbWrapperRegistry::instance().addMethod(\"" + gParent + "<$>\", \"" + name + "\", " + regname + ");";
        } else {
            regDecl = regHeader + " { return fromPy<" + gParent +"*>(_self)->_" + name + "(_self, _linargs, _kwds); }";
            regCall = "PbWrapperRegistry::instance().addMethod(\"" + gParent + "\", \"" + name + "\", " + regname + ");";
        }
    }
    if (isHeader) {
        gRegText += regDecl + "\n" + regCall + "\n";
    } else {
        gLocalReg += regDecl + nl;
        gRegText += "extern " + regHeader + ";\n" + regCall + "\n";
    }
    
    if (gDocMode) {
        caller = "";
        func = "PYTHON " + func;
    }    
    return buildline(lb) + func + nl + caller;
}

string processPythonVariable(int lb, const string& name, const vector<Argument>& opts, const string& type, int line) {
    // beautify code
    string nl = gDebugMode ? "\n" : "";
    string tb = gDebugMode ? "\t" : "";
    
    // is header file ?
    bool isHeader = gFilename[gFilename.size()-2] == '.' && gFilename[gFilename.size()-1] == 'h';
    
    if (gParent.empty())
        errMsg(line, "PYTHON variables con only be used inside classes");
    
    string pname = name;
    
    // process options
    for (size_t i=0; i<opts.size(); i++) {
        if (opts[i].name == "name") 
            pname = opts[i].value;
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
    string getter = gethdr+" { return toPy(fromPy<" + gParent+"*>(self)->" + name + "); }";
    string setter = sethdr+" { fromPy<" + gParent+"*>(self)->" + name + "=fromPy<" + type + " >(val); return 0;}";
        
    // register
    gRegText += "PbWrapperRegistry::instance().addGetSet(\""+gParent+"\",\""+pname+"\",_get_"+nn+ ",_set_"+nn+");\n";
    if (isHeader) {
        gRegText += getter+"\n"+setter+"\n";
    } else {
        gLocalReg += getter + nl + setter + nl;
        gRegText += "extern " + gethdr + ";\nextern " + sethdr + ";\n";
    }
    
    return buildline(lb) + code;
}

string processPythonClass(int lb, const string& name, const vector<Argument>& opts, const std::vector<Argument>& templArgs, const string& baseclass, const string& code, int line) {
    // beautify code
    string nl = gDebugMode ? "\n" : "";
    string tb = gDebugMode ? "\t" : "";
    
    // is header file ?
    bool isHeader = gFilename[gFilename.size()-2] == '.' && gFilename[gFilename.size()-1] == 'h';
    
    if (!isHeader && !templArgs.empty())
        errMsg(line, "PYTHON template classes can only be defined in header files.");
    
    string pname = name; // for now
    
    // process options
    for (size_t i=0; i<opts.size(); i++) {
        if (opts[i].name == "name") 
            pname = opts[i].value;
        else
            errMsg(line, "PYTHON(opt): illegal kernel option. Supported options are: 'name'");
    }
    
    // class registry
    string registry = "";
    if (templArgs.empty())
        registry = "PbWrapperRegistry::instance().addClass(\"" + pname + "\", \"" + name + "\", \"" + baseclass + "\");";
    else
        registry = "@template " + name + " PbWrapperRegistry::instance().addClass(\"@\", \"" + name + "<$>\", \"" + baseclass + "\");";
    
    // register class
    if (isHeader) {
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
        if (isHeader)
            gRegText += createConverters(name, "", " ", "\n");
        else
            gLocalReg += createConverters(name, tb, nl, nl);
    }
    
    // tokenize and parse contained python functions
    gParent = name;
    gIsTemplated = !templArgs.empty();
    gFoundConstructor = false;
    string newText = processText(code.substr(1), line);
    gParent = "";    
    if (!gFoundConstructor)
        errMsg(line, "no PYTHON constructor found in class '" + name + "'");
    if (!isHeader && gIsTemplated)
        errMsg(line, "PYTHON class template can only be used in header files.");
    
    // create class
    string pclass = "";
    if (gDocMode) {
        pclass += "//! \\ingroup PyClasses\nPYTHON ";
    }
    if (gIsTemplated)
        pclass += "template<" + listArgs(templArgs) + "> " + nl;
    pclass += "class " + name + " : public " + baseclass + " {" + nl;
    pclass += newText + nl;
    if (!gDocMode) {
        pclass += "protected:" + nl;
        pclass += tb+"PbArgs _args;" + nl;
    }
    pclass += "};" + nl;
    
    
    return buildline(lb) + pclass + gLocalReg;
}

set<string> gAliasRegister;
string processPythonInstantiation(int lb, const string& name, const std::vector<Argument>& templArgs, const string& aliasname, int) {
    gRegText += "@instance " + name + " " + listArgs(templArgs) + " " + aliasname + "\n";
    gRegText += createConverters(name + "<" + listArgs(templArgs) + ">", "", " ", "\n");
    
    if (gAliasRegister.find(aliasname) == gAliasRegister.end()) {
        gAliasRegister.insert(aliasname);
        return buildline(lb) + "typedef " + name + "<" + listArgs(templArgs) + "> " + aliasname + "; ";
    }
    return buildline(lb);
}