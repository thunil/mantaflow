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

string processKernel(int lb, const string& kname, const ArgList& opts, Argument retType, const ArgList& returnArg, const ArgList& templArgs, const ArgList& args, const string& code, int line) {
    // beautify code
	string nlr = "\n";
    string nl = gDebugMode ? "\n" : " ";
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
    string bnd = "0", reduceOp="";
    for (size_t i=0; i<opts.size(); i++) {
        if (opts[i].name == "ijk") 
            idxMode = false;
        else if (opts[i].name == "index" || opts[i].name == "idx")
            idxMode = true;
        else if (opts[i].name == "st" || opts[i].name == "single")
            gMTType = MTNone;
        else if (opts[i].name == "pts" || opts[i].name == "particle" || opts[i].name == "points")
            pts = true;
        else if (opts[i].name == "bnd")
            bnd = opts[i].value;
        else if (opts[i].name == "bnd")
            bnd = opts[i].value;
        else if (opts[i].name == "reduce") {
            reduce = true;
            reduceOp = opts[i].value;
            assert(reduceOp == "+" || reduceOp == "-" || reduceOp == "*" || reduceOp == "/" || reduceOp == "min" || reduceOp == "max",
                   "invalid 'reduce' operator. Expected reduce= +|-|*|/|min|max");
        } else
            errMsg(line, "KERNEL(opt): illegal kernel option. Supported options are: 'ijk', 'idx', 'bnd=x', 'reduce', 'st', 'pts'");
    }
    string qualifier = (!reduce && gMTType==MTTBB) ? "const" : "";
    string tbbcall = reduce ? "tbb::parallel_reduce" : "tbb::parallel_for";
    const bool haveOuter = !reduce && !returnArg.empty() && gMTType == MTTBB;
    
	// point out illegal paramter combinations
    if (bnd != "0" && idxMode)
        errMsg(line, "KERNEL(opt): can't combine index mode with bounds iteration.");    
    if (pts && (idxMode || bnd != "0" ))
        errMsg(line, "KERNEL(opt): Modes 'ijk', 'idx' and 'bnd' can't be applied to particle kernels.");

    // check type consistency of first 'returns' with return type
    if (!returnArg.empty() && retType.type != "void") {
        assert(returnArg.size() == 1, "multiple returns statement only work for 'void' kernels");
        const Argument& rt = returnArg[0];
        assert(rt.type == retType.type && rt.templ == retType.templ &&
               rt.isRef == retType.isRef && rt.isPointer == retType.isPointer &&
               rt.isConst == retType.isConst, "return type does not match type in first 'returns' statement");
    }
    if (retType.type != "void" && returnArg.empty())
        errMsg(line, "return argument specified without matching 'returns' initializer");
    
    // parse arguments
    string initList = "", argList = "", copier = "", members = "", basegrid="", baseobj="", callList = "", orgArgList="", mCallList, outerCopier="", lineAcc="";
    for (size_t i=0; i<args.size(); i++) {
        string type = args[i].getType();
        string name = args[i].name;
        
        initList += (i==0) ? "" : ", ";
        argList += (i==0) ? "" : ", ";
        copier += (i==0) ? "" : ", ";
        callList += (i==0) ? "" : ", ";
        orgArgList += (i==0) ? "" : ", ";
        outerCopier += (i==0) ? "" : ", ";
        
        string aname = "_" + name;
        string sname = "m_" + name;
        copier += sname + "(o." + sname + ")";
        outerCopier += "o._inner." + sname;
        members += tb+ type + " " + sname + ";" + nl;
        initList += sname + "(" + aname + ")";
        callList += aname;
        mCallList += ", "+sname;
        orgArgList += type + " " + name;
        argList += type + " " + aname;
        if (!args[i].value.empty())
            argList += " = " + args[i].value;
        if (i<3) {
            Argument refa(args[i]), refb(args[i]);
            refa.isRef = true;
            refb.isRef = false;
            char num = '0' + i;
            string avar = (haveOuter ? "_inner." : "") + sname;
            lineAcc += tb + refa.getType() + " getArg" + num + "() { return " + avar + "; }" + nl;
            //lineAcc += tb + "void setArg" + num + "(" + refa.getType() + " _v) { " + avar + " = _v; }" + nl;
            lineAcc += tb + "typedef " + refb.getType() + " type" + num + ";" + nl;
        }
        
        // figure out basegrid and baseobj
        if (basegrid.empty()) {
            if (!pts && type.find("Grid") != string::npos) {
                if (!args[i].isPointer) basegrid += "&";
                basegrid += "_" + name;
            } else if (pts) {
                if (args[i].isPointer) basegrid += "*";
                basegrid += "_" + name;
            }
        }
        if (baseobj.empty() && (type.find("Grid") != string::npos || type.find("System") != string::npos)) {
            if (args[i].isPointer) baseobj += "*";
            baseobj += "_" + name;
        }
    }
    if (basegrid.empty()) 
        errMsg(line, "KERNEL: use at least one grid to call the kernel.");
    
    // include parent if possible
    string initParent = "", copyParent = "", parentMember = "";
    if (!baseobj.empty()) {
        initParent = "parent(("+baseobj+").getParent()), ";
        copyParent = "parent(o.parent), ";
        parentMember = tb + "FluidSolver* parent;" + nl;
    }
    
    string kclass = "", kclassname = kname, callerClass="";
    string outerMembers = parentMember, outerArgList = argList, initRetval = "", retList = "";
        
    // add return args as members
    for (size_t i=0; i<returnArg.size(); i++) {
        Argument arg = returnArg[i];
        
        outerMembers += tb + arg.getTypeName() + ";"+nl;
        callList += ", " + arg.name;
        copier +=", " + arg.name + "(" + arg.value + ")";
        initRetval += arg.name + "(" + arg.value + "), ";
        retList += (i==0) ? "" : ", ";
        retList += arg.name;
            
        if (haveOuter) {
            initList +=", m_" + arg.name + "(_"+arg.name+ ")";\
            arg.isRef = true;
            arg.name = "_" + arg.name;
            argList += ", " + arg.getTypeName();
            arg.name = "m" + arg.name;          
        } else {
            initList +=", " + arg.name + "(" + arg.value + ")";
        }
        members += tb + arg.getTypeName() + ";"+nl;        
        
        mCallList += ", " + arg.name;
        // ref it
        Argument arg2 = returnArg[i];
        arg2.isRef = true;
        orgArgList += ", " + arg2.getTypeName();         
    }

    // define return conversion operator
    string outOp = "";
    if (retType.type != "void") {
        outOp += tb + "operator " + retType.complete + " () {" + nl;
        outOp += tb2+ "return " + returnArg[0].name + ";" + nl;
        outOp += tb + "}" + nl;
        Argument crRet = retType;
        crRet.isRef = true; crRet.isConst = true;
        outOp += tb + crRet.getType() + " getRet() const { return " + returnArg[0].name + "; }" + nl;
        outOp += lineAcc;
    }
        
    // create outer class for non-reduce return values
    if (haveOuter) {
        kclassname = "_kernel_" + kname;
        if (!templArgs.empty()) callerClass += "template <" + listArgs(templArgs) + ">" + nl;
        callerClass += "struct " + kname + " : public " + (pts ? "Particle" : "") + "KernelBase { " + nl;
        callerClass += tb + kname + " ( " + outerArgList + ") : ";
        if (pts)
            callerClass += "ParticleKernelBase((" + basegrid + ").size()), ";
        else
            callerClass += "KernelBase(" + basegrid + ", " + bnd + "), ";
        callerClass += initParent + initRetval + " _inner(" + callList + ") { }" + nl;
        
        // copy constructor
        callerClass += tb+ kname + " (const " + kname + "& o) : ";
        if (!pts)
            callerClass += "KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), ";
        else
            callerClass += "ParticleKernelBase(o.size), ";
        callerClass += copyParent + initRetval + " _inner(" + outerCopier + "," + retList + ") {}" + nl;
    
        callerClass += tb + "void run() { _inner.run(); }" + nl;
        callerClass += outOp;
        callerClass += outerMembers;
        // declar inner class
        callerClass += tb + kclassname;
        if (!templArgs.empty()) {
            callerClass += "<" + templArgs[0].name;
            for (size_t i=1; i<templArgs.size(); i++) callerClass += ", "+templArgs[i].name;
            callerClass += ">";
        }
        callerClass += " _inner;" + nl;
        callerClass += "};" + nl;        
    }
    
    // create kernel class
    if (!templArgs.empty()) kclass += "template <" + listArgs(templArgs) + ">" + nl;
    kclass += "struct " + kclassname + " : public " + (pts ? "Particle" : "") + "KernelBase { " + nl;
    
    // init constructor
    kclass += tb+ kclassname + " (" + argList + ") :" + nl;
    if (pts)
        kclass += tb2+ "ParticleKernelBase((" + basegrid + ").size()), " + initParent + initList + nl;
    else
        kclass += tb2+ "KernelBase(" + basegrid + ", " + bnd + "), " + initParent + initList + nl;
    kclass += tb + "{" + nl;
    kclass += tb2+ "run();" + nl;
    kclass += tb+ "}" + nl;
    
    // copy constructor
    if (!haveOuter) {
        kclass += tb+ kclassname + " (const " + kclassname + "& o) : ";
        if (!pts)
            kclass += "KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z), ";
        else
            kclass += "ParticleKernelBase(o.size), ";
        kclass += copyParent + copier + " {}" + nl;
    }
    
    // code point
    if (pts)
        kclass+= tb+ "inline void op(int i,"+orgArgList+") "+qualifier+nl;
    else if (idxMode)
        kclass+= tb+ "inline void op(int idx,"+orgArgList+") "+qualifier+nl;
    else
        kclass+= tb+ "inline void op(int i,int j,int k,"+orgArgList+") "+qualifier+nl;
    kclass += code + nl;
    
    // create main kernel function 
    if (gMTType == MTTBB) {
        // multithreading using intel tbb
        kclass += tb+ "void operator() (const tbb::blocked_range<size_t>& r) " + qualifier + "{" + nl;
        if (pts) {
            kclass += tb2+ "for (int i=r.begin(); i!=(int)r.end(); i++)" + nl;
            kclass += tb3+ "op(i"+mCallList+");" + nl;
        } else if (idxMode) {
            kclass += tb2+ "for (int idx=r.begin(); idx!=(int)r.end(); idx++)" + nl;
            kclass += tb3+"op(idx"+mCallList+");" + nl;
        } else {
            kclass += tb2+ "const int _maxX = maxX, _maxY=maxY;" + nl;
            kclass += tb2+ "if (maxZ>1) {" +nl;
            kclass += tb3+ "for (int k=r.begin(); k!=(int)r.end(); k++)" + nl;
            kclass += tb3+ "for (int j=" + bnd + "; j<_maxY; j++)" + nl;
            kclass += tb3+ "for (int i=" + bnd + "; i<_maxX; i++)" + nl;
            kclass += tb4+"op(i,j,k"+mCallList+");" + nl;
            kclass += tb2+ "} else {" + nl;
            kclass += tb3+ "const int k=0;" + nl;
            kclass += tb3+ "for (int j=r.begin(); j!=(int)r.end(); j++)" + nl;
            kclass += tb3+ "for (int i=" + bnd + "; i<_maxX; i++)" + nl;
            kclass += tb4+"op(i,j,k"+mCallList+");" + nl;
            kclass += tb2+ "}";
        }            
        kclass += tb+"}" + nl;
        kclass += tb+ "void run() {" + nl;
        if (pts)
            kclass += tb2+ tbbcall + "(tbb::blocked_range<size_t>(0, size), *this);"+ nl;
        else if (idxMode)
            kclass += tb2+ tbbcall + "(tbb::blocked_range<size_t>(0, maxCells), *this);"+ nl;
        else {
            kclass += tb2+ "if (maxZ>1) " + nl;
            kclass += tb3+ tbbcall + "(tbb::blocked_range<size_t>(minZ, maxZ), *this);"+ nl;
            kclass += tb2+ "else " + nl;
            kclass += tb3+ tbbcall + "(tbb::blocked_range<size_t>("+ bnd +", maxY), *this);"+ nl;
        }
        kclass += tb+ "}" + nl;
        if (reduce) {
            // split constructor
            kclass += tb+ kclassname + " (" + kclassname + "& o, tbb::split) : " + nl;
            if (!pts)
                kclass += tb2+ "KernelBase(o.maxX, o.maxY, o.maxZ, o.maxCells, o.minZ, o.X, o.Y, o.Z)," + nl;
            else
                kclass += tb2+ "ParticleKernelBase(o.size)," + nl;
            kclass += tb2+ copyParent + copier + " {}" + nl;
            // join
            kclass += tb + "void join(const " + kclassname+"& _other) {" + nl;
            for (size_t i=0; i<returnArg.size(); i++) {
                if (reduceOp == "min" || reduceOp == "max")
                    kclass += tb2 + returnArg[i].name + " = " + reduceOp + "(" + returnArg[i].name + ", _other." + returnArg[i].name +");"+nl;
                else
                    kclass += tb2 + returnArg[i].name + " " + reduceOp + "= _other." + returnArg[i].name +";"+nl;
            }
            kclass += tb + "}" + nl;
        }
	}
    else if(gMTType == MTOpenMP)
    {
		// correct for the introduced linebreaks
        if (!gDebugMode) {
            ostringstream reline;
            reline << nlr << "#line " << (line) << nlr;
            callerClass += reline.str();
        }
        
        string directive = "#pragma omp for", postDir = "", preDir = "";
        if (reduce) {
            directive = "#pragma omp for nowait";
            for (size_t i=0; i<returnArg.size(); i++)
                preDir += tb3 + returnArg[i].getType() + " _loc_" + returnArg[i].name + " = " + returnArg[i].value + ";" + nl;
            postDir = nlr + tb3 + "#pragma omp critical" + nlr;
            postDir += tb3 + "{" + nl;
            for (size_t i=0; i<returnArg.size(); i++) {
                if (reduceOp == "*" || reduceOp == "/" || reduceOp == "+" || reduceOp == "-")
                    postDir += tb4 + returnArg[i].name + " " + reduceOp + "= _loc_" + returnArg[i].name + ";" + nl;                
                else {
                    postDir += tb4 + returnArg[i].name + " = " + reduceOp + "(" + returnArg[i].name + ", _loc_" + returnArg[i].name + ");" + nl;                
                }                
            }
            postDir += tb3 + "}" + nl;                
            
            // rebuild call list with local return args
            mCallList = "";
            for (size_t i=0; i<args.size(); i++) mCallList += ", m_" + args[i].name;
            for (size_t i=0; i<returnArg.size(); i++) mCallList += ", _loc_" + returnArg[i].name;
        }
        
        kclass += tb+ "void run() {" + nl;
        if (pts) {
            kclass += tb2+ "const int _sz = size;"+ nlr;
            kclass += tb2+ "#pragma omp parallel" + nlr;
            kclass += tb2+ "{" + nl + preDir + nlr;         
            kclass += tb3+ directive + nlr;
            kclass += tb3+ "for (int i=0; i < _sz; i++)" + nl;
            kclass += tb4+ "op(i"+mCallList+");" + nl;
            kclass += postDir + tb2 + "}" + nl;
        } else if (idxMode) {
            kclass += tb2+ "const int _maxCells = maxCells;"+ nlr;
            kclass += tb2+ "#pragma omp parallel" + nlr;
            kclass += tb2+ "{" + nl + preDir + nlr;                     
            kclass += tb3+ directive + nlr;                     
            kclass += tb3+ "for (int idx=0; idx < _maxCells; idx++)" + nl;
            kclass += tb4+ "op(idx"+mCallList+");" + nl;
            kclass += postDir + tb2 + "}" + nl;
        } else {
            kclass += tb2+ "const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ;" + nl;
            kclass += tb2+ "if (maxZ>1) {"+nlr;
            kclass += tb3+ "#pragma omp parallel" + nlr;
            kclass += tb3+ "{" + nl + preDir + nlr;                     
            kclass += tb4+ directive + nlr;
            kclass += tb4+ "for (int k=minZ; k < _maxZ; k++)" + nl;
            kclass += tb4+ "for (int j="+bnd+"; j < _maxY; j++)" + nl;
            kclass += tb4+ "for (int i="+bnd+"; i < _maxX; i++)" + nl;
            kclass += tb5+ "op(i,j,k"+mCallList+");" + nl;
            kclass += postDir + tb3 + "}";        
            kclass += tb2+"} else {"+nl;
            kclass += tb3+ "const int k=0;"+nlr;
            kclass += tb3+ "#pragma omp parallel" + nlr;
            kclass += tb3+ "{" + nl + preDir + nlr;
            kclass += tb4+ directive + nlr;
            kclass += tb4+ "for (int j="+bnd+"; j < _maxY; j++)" + nl;
            kclass += tb4+ "for (int i="+bnd+"; i < _maxX; i++)" + nl;
            kclass += tb5+ "op(i,j,k"+mCallList+");" + nl;            
            kclass += postDir + tb3 + "}";
            kclass += tb2 + "}";        
        }
        kclass += tb+"}"+nl;
    }
    else
    {
        kclass += tb+ "void run() {" + nl;
        if (pts) {
            kclass += tb2+ "const int _sz = size;"+ nl;
            kclass += tb2+ "for (int i=0; i < _sz; i++)" + nl;
            kclass += tb3+ "op(i"+mCallList+");" + nl;
        } else if (idxMode) {
            kclass += tb2+ "const int _maxCells = maxCells;"+ nl;
            kclass += tb2+ "for (int idx=0; idx < _maxCells; idx++)" + nl;
            kclass += tb3+ "op(idx"+mCallList+");" + nl;
        } else {
            kclass += tb2+ "const int _maxX = maxX, _maxY=maxY, _maxZ = maxZ;" + nl;
            kclass += tb2+ "for (int k=minZ; k < _maxZ; k++)" + nl;
            kclass += tb3+ "for (int j="+bnd+"; j < _maxY; j++)" + nl;
            kclass += tb4+ "for (int i="+bnd+"; i < _maxX; i++)" + nl;
            kclass += tb3+ "op(i,j,k"+mCallList+");" + nl;
        }
        kclass += tb+"}" + nl;
    }
    if (!haveOuter) kclass += outOp;
    
    kclass += parentMember + members + nl;
    kclass += "};" + nl;
   
	return buildline(lb) + kclass + callerClass;
}
