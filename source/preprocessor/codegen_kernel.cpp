/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor: Process replacement text of KERNEL keyword
 *
 ******************************************************************************/

#include "prep.h"
#include <cstdlib>
#include <set>
#include <sstream>
#include <iostream>
using namespace std;

string processKernel(int lb, const string& kname, const ArgList& opts, Argument retType, const ArgList& templArgs, const ArgList& args, const string& codebase, int line, bool isClass) {
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
    bool idxMode = false, reduce = false, pts = false;
    string bnd = "0";
    for (size_t i=0; i<opts.size(); i++) {
        if (opts[i].name == "ijk") 
            idxMode = false;
        else if (opts[i].name == "index" || opts[i].name == "idx")
            idxMode = true;
        else if (opts[i].name == "st" || opts[i].name == "single")
            gMTType = None;
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
    string loader = "", initList = "", argList = "", copier = "", members = "", basegrid="", callList = "";
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
        callList += (i==0) ? "" : ", ";
        
        string aname = "_" + name;
        string sname = (isClass ? "" : "m_" )+ name;
        loader += tb2+ type + " " + name + " = " + sname + ";" + nl;
        copier += sname + "(o." + sname + ")";
        members += tb+ type + " " + sname + ";" + nl;
        initList += sname + "(" + aname + ")";
        callList += aname;
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
    
    string kclass = "", kclassname = kname, caller="";
    
    // create return type adapter
    if (!isClass && retType.type != "void") {
        if (retType.isPointer || retType.isRef)
            errMsg(line, "KERNEL: return type cannot be pointer/ref");
        
        // rename kernel 
        kclassname = "_kernel_" + kname;
        
        // build caller function with provided ret arg
        caller += "inline " + retType.complete + " " + kname + "(" + argList + ", " + retType.complete + "& __ret) {" + nl;
        caller += tb + kclassname + " __kn (" + callList + ",__ret);" + nl;
        caller += tb + "return __ret;" + nl;
        caller += "}" + nl;
        
        // build caller function with generated ret arg
        caller += "inline " + retType.complete + " " + kname + "(" + argList + ") {" + nl;
        caller += tb + retType.complete + " __ret((" + basegrid + ")->getParent());" + nl;
        caller += tb + kclassname + " __kn (" + callList + ",__ret);" + nl;
        caller += tb + "return __ret;" + nl;
        caller += "}" + nl;
        
        // additional argument
        argList += "," + retType.complete + "& _ret";
        initList += ", _m_ret(_ret) ";
        members += tb+ retType.complete + "& _m_ret;" + nl;
        loader += tb2 + retType.complete + "& _ret = _m_ret;"  + nl;        
    }
        
    // create kernel class
    if (!templArgs.empty()) kclass += "template <" + listArgs(templArgs) + ">" + nl;
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
    
    // create main kernel function 
    if (gMTType == TBB) {
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
        if (pts)
            kclass += tb2+ tbbcall + "(tbb::blocked_range<size_t>(0, size), *this);"+ nl;
        else if (idxMode)
            kclass += tb2+ tbbcall + "(tbb::blocked_range<size_t>(0, maxCells), *this);"+ nl;
        else
            kclass += tb2+ tbbcall + "(tbb::blocked_range<size_t>(minZ, maxZ), *this);"+ nl;
        kclass += tb+ "}" + nl;
	} else if(gMTType == OpenMP) {
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
            kclass += tb2+ "for (int k=minZ; k < _maxZ; k++)" + nl;
            kclass += tb3+ "for (int j="+bnd+"; j < _maxY; j++)" + nl;
            kclass += tb4+ "for (int i="+bnd+"; i < _maxX; i++)" + nl;
            kclass += isClass ? (tb5+"(*this)(i,j,k);") : code;
        }
        kclass += nl + tb+"}" + nl;
    }
    if (reduce) {
        if (gMTType == TBB) {
            // split constructor
            kclass += tb+ kclassname + " (" + kclassname + "& o, tbb::split) : " + nl;
            if (!pts)
                kclass += tb2+ "KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z)," + nl;
            else
                kclass += tb2+ "KernelBase(o.size)," + nl;
            kclass += tb2+ copier + nl;
            kclass += tb + "{" + nl; 
            kclass += tb2+ "setup();" + nl;
            kclass += tb+ "}" + nl;
        }
    }    
    kclass += nl + members;
    kclass += isClass ? code.substr(1) : "}";
    kclass += ";" + nl;
   
    return buildline(lb) + kclass + caller;
}
