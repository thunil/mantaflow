/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
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

#define kernelAssert(x,msg) if(!(x)){errMsg(line,msg);}

string processKernel(const Block& block, const string& code) {
    // beautify code
	string nlr = "\n";
    string nl = gDebugMode ? "\n" : " ";
    string tb = gDebugMode ? "\t" : "";    
    string tb2 = tb+tb, tb3=tb2+tb, tb4=tb3+tb, tb5=tb4+tb;

    const List<Argument>& args = block.func.arguments;
    const string kernelName = block.func.name;
    const Type& retType = block.func.returnType;
    const List<Type>& templArgs = block.func.templateTypes;
    const List<Argument>& opts = block.options;
    const List<Argument>& returnArg = block.reduceArgs;
    const int line = block.line0;
    
    if (gDocMode) {
        string ds = "//! \\ingroup Kernels\nKERNEL<";
        for(size_t i=0; i<opts.size(); i++) { if (i!=0) ds+=", "; ds+=opts[i].minimalText; }
        ds += "> " + kernelName + " (";
        for(size_t i=0; i<args.size(); i++) { if (i!=0) ds+=", "; ds+=args[i].minimalText;}
        return ds + " ) {}\n";
    }
    
    // process options
    bool idxMode = false, reduce = false, pts = false;
    string bnd = "0", reduceOp="";
	MType mtType = gMTType;
    for (size_t i=0; i<opts.size(); i++) {
        if (opts[i].name == "ijk") 
            idxMode = false;
        else if (opts[i].name == "index" || opts[i].name == "idx")
            idxMode = true;
        else if (opts[i].name == "st" || opts[i].name == "single")
            mtType = MTNone;
        else if (opts[i].name == "pts" || opts[i].name == "particle" || opts[i].name == "points")
            pts = true;
        else if (opts[i].name == "bnd")
            bnd = opts[i].value;
        else if (opts[i].name == "bnd")
            bnd = opts[i].value;
        else if (opts[i].name == "reduce") {
            reduce = true;
            reduceOp = opts[i].value;
            kernelAssert(reduceOp == "+" || reduceOp == "-" || reduceOp == "*" || reduceOp == "/" || reduceOp == "min" || reduceOp == "max",
                   "invalid 'reduce' operator. Expected reduce= +|-|*|/|min|max");
        } else
            errMsg(line, "KERNEL(opt): illegal kernel option '"+ opts[i].name +"' Supported options are: 'ijk', 'idx', 'bnd=x', 'reduce', 'st', 'pts'");
    }
    string qualifier = (!reduce && mtType==MTTBB) ? "const" : "";
    string tbbcall = reduce ? "tbb::parallel_reduce" : "tbb::parallel_for";
    const bool haveOuter = !reduce && !returnArg.empty() && mtType == MTTBB;
    
	// point out illegal paramter combinations
    if (bnd != "0" && idxMode)
        errMsg(line, "KERNEL(opt): can't combine index mode with bounds iteration.");    
    if (pts && (idxMode || bnd != "0" ))
        errMsg(line, "KERNEL(opt): Modes 'ijk', 'idx' and 'bnd' can't be applied to particle kernels.");

    // check type consistency of first 'returns' with return type
    if (!returnArg.empty() && retType.name != "void") {
        kernelAssert(returnArg.size() == 1, "multiple returns statement only work for 'void' kernels");
        const Type& rt = returnArg[0].type;
        kernelAssert(rt == retType, "return type does not match type in first 'returns' statement");
    }
    if (retType.name != "void" && returnArg.empty())
        errMsg(line, "return argument specified without matching 'returns' initializer");
    
    // parse arguments
    string initList = "", argList = "", copier = "", members = "", basegrid="", baseobj="", callList = "", orgArgList="", mCallList, outerCopier="", lineAcc="";
    for (size_t i=0; i<args.size(); i++) {
        string type = args[i].type.build();
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
            Type refa(args[i].type), refb(args[i].type);
            refa.isRef = true;
            refb.isRef = false; refb.isConst = false;
            char num = '0' + i;
            string avar = (haveOuter ? "_inner." : "") + sname;
            lineAcc += tb + refa.build() + " getArg" + num + "() { return " + avar + "; }" + nl;
            //lineAcc += tb + "void setArg" + num + "(" + refa.type.build() + " _v) { " + avar + " = _v; }" + nl;
            lineAcc += tb + "typedef " + refb.build() + " type" + num + ";" + nl;
        }
        
        // figure out basegrid and baseobj
        if (basegrid.empty()) {
            if (!pts && type.find("Grid") != string::npos) {
                if (!args[i].type.isPointer) basegrid += "&";
                basegrid += "_" + name;
            } else if (pts) {
                if (args[i].type.isPointer) basegrid += "*";
                basegrid += "_" + name;
            }
        }
        if (baseobj.empty() && (type.find("Grid") != string::npos || type.find("System") != string::npos)) {
            if (args[i].type.isPointer) baseobj += "*";
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
    
    string kclass = "", kclassname = kernelName, callerClass="";
    string outerMembers = parentMember, outerArgList = argList, initRetval = "", retList = "";
        
    // add return args as members
    for (size_t i=0; i<returnArg.size(); i++) {
        Argument arg = returnArg[i];
        
        outerMembers += tb + arg.type.build() + " " + arg.name + ";"+nl;
        callList += ", " + arg.name;
        copier +=", " + arg.name + "(" + arg.value + ")";
        initRetval += arg.name + "(" + arg.value + "), ";
        retList += (i==0) ? "" : ", ";
        retList += arg.name;
            
        if (haveOuter) {
            initList +=", m_" + arg.name + "(_"+arg.name+ ")";\
            arg.type.isRef = true;
            arg.name = "_" + arg.name;
            argList += ", " + arg.type.build();
            arg.name = "m" + arg.name;          
        } else {
            initList +=", " + arg.name + "(" + arg.value + ")";
        }
        members += tb + arg.type.build() + " " + arg.name + ";"+nl;        
        
        mCallList += ", " + arg.name;
        // ref it
        Argument arg2 = returnArg[i];
        arg2.type.isRef = true;
        orgArgList += ", " + arg2.type.build();         
    }

    // define return conversion operator
    string outOp = "";
    if (retType.name != "void") {
        outOp += tb + "operator " + retType.minimal + " () {" + nl;
        outOp += tb2+ "return " + returnArg[0].name + ";" + nl;
        outOp += tb + "}" + nl;
        Type crRet = retType;
        crRet.isRef = true; crRet.isConst = true;
        outOp += tb + crRet.build() + " getRet() const { return " + returnArg[0].name + "; }" + nl;
        outOp += lineAcc;
    }
        
    // create outer class for non-reduce return values
    if (haveOuter) {
        kclassname = "_kernel_" + kernelName;
        if (!templArgs.empty()) callerClass += "template " + templArgs.minimal + nl;
        callerClass += "struct " + kernelName + " : public " + (pts ? "Particle" : "") + "KernelBase { " + nl;
        callerClass += tb + kernelName + " ( " + outerArgList + ") : ";
        if (pts)
            callerClass += "ParticleKernelBase((" + basegrid + ").size()), ";
        else
            callerClass += "KernelBase(" + basegrid + ", " + bnd + "), ";
        callerClass += initParent + initRetval + " _inner(" + callList + ") { }" + nl;
        
        // copy constructor
        callerClass += tb+ kernelName + " (const " + kernelName + "& o) : ";
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
    if (!templArgs.empty()) kclass += "template " + templArgs.minimal + nl;
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
    if (mtType == MTTBB) {
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
    else if(mtType == MTOpenMP)
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
                preDir += tb3 + returnArg[i].type.build() + " _loc_" + returnArg[i].name + " = " + returnArg[i].value + ";" + nl;
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
		// init kernel thread info fields
		preDir += tb3 + "this->threadId = omp_get_thread_num(); this->threadNum = omp_get_num_threads();" + nl;
        
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

	debMsg( line, "Kernel summary '"<< kernelName <<"'. Basegrid: "<< basegrid <<", baseobj: "<<baseobj<<", mt: "<< mtType );
   
	return block.linebreaks() + kclass + callerClass;
}
