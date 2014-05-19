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

#define STR(x) #x

//******************************************************
// Templates for code generation

// TP: why do we need getArg? just directly access the argument via its name...
const string TmpAccessor = STR(
inline $TYPE$ getArg$IDX$() { return $NAME$; }
typedef $TYPE$ type$IDX$;
);

// Single kernel, default
const string TmpSingleKernel = STR(
$TEMPLATE$ struct $KERNEL$ : public KernelBase {
    $KERNEL$($ARGS$) : 
@IF(PTS) 
        KernelBase($BASE$.size()) $INIT$ $LOCALSET$
@ELSE
        KernelBase($BASE$,$BND$) $INIT$ $LOCALSET$
@END
    {
        run();
    }
@IF(IJK)
    inline void op(int i, int j, int k, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@ELSE
    inline void op(int idx, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@END

@IF(RET_NAME)
    inline operator $RET_TYPE$() { return $RET_NAME$; }
    inline $RET_TYPE$ & getRet() { return $RET_NAME$; }
@END
    $ACCESSORS$

    $RUN$
    $MEMBERS$
    $LOCALS$
};
);

// Necesary for TBB with nontrivial return values 
const string TmpDoubleKernel = STR(
// inner kernel
$TEMPLATE$ struct _$KERNEL$ : public KernelBase {
    _$KERNEL$(const KernelBase& base, $ARGS$ $LOCALARG$) : 
        KernelBase(base) $INIT$ $LOCALINIT${}

@IF(IJK)
    inline void op(int i, int j, int k, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@ELSE
    inline void op(int idx, $ARGS$ $LOCALARG$) $CONST$ $CODE$
@END
    $RUN$
    $MEMBERS$
    $LOCALS_REF$
};

// outer kernel with accessors
$TEMPLATE$ struct $KERNEL$ : public KernelBase {
    $KERNEL$($ARGS$) :
@IF(PTS) 
        KernelBase($BASE$.size()) $COMMA$ _inner(KernelBase($BASE$.size()),$CALL$)
@ELSE
        KernelBase($BASE$,$BND$) $COMMA$ _inner(KernelBase($BASE$,$BND$),$CALL$)
@END
        $INIT$ $LOCALSET$
    {
        run();
    }

    void run() { _inner.run(); }

@IF(RET_NAME)
    inline operator $RET_TYPE$() { return $RET_NAME$; }
    inline $RET_TYPE$ & getRet() { return $RET_NAME$; }
@END
    $ACCESSORS$
    $MEMBERS$
    $LOCALS$
    _$KERNEL$$TPL$ _inner;
};
);

const string TmpRunSimple = STR(
void run() {
@IF(IJK)
    const int _maxX = maxX; 
    const int _maxY = maxY;
    for (int k=minZ; k< maxZ; k++)
    for (int j=$BND$; j< _maxY; j++)
    for (int i=$BND$; i< _maxX; i++)
        op(i,j,k, $CALL$);
@ELSE
    const int _sz = size;
    for (int i=0; i < _sz; i++)
        op(i, $CALL$);
@END
}
);

const string TmpRunTBB = STR(
void operator() (const tbb::blocked_range<size_t>& r) $CONST$ {
@IF(IJK)
    const int _maxX = maxX;
    const int _maxY = maxY;
    if (maxZ>1) {
        for (int k=r.begin(); k!=(int)r.end(); k++)
        for (int j=$BND$; j<_maxY; j++)
        for (int i=$BND$; i<_maxX; i++)
            op(i,j,k,$CALL$);
    } else {
        const int k=0;
        for (int j=r.begin(); j!=(int)r.end(); j++)
        for (int i=$BND$; i<_maxX; i++)
            op(i,j,k,$CALL$);
    }
@ELSE
    for (int idx=r.begin(); idx!=(int)r.end(); idx++)
        op(idx, $CALL$);
@END
}
void run() {
@IF(IJK)
    if (maxZ>1)
        tbb::parallel_$METHOD$ (tbb::blocked_range<size_t>(minZ, maxZ), *this);
    else
        tbb::parallel_$METHOD$ (tbb::blocked_range<size_t>($BND$, maxY), *this);
@ELSE
    tbb::parallel_$METHOD$ (tbb::blocked_range<size_t>(0, size), *this);
@END
}
@IF(REDUCE)
    $IKERNEL$ ($IKERNEL$& o, tbb::split) : KernelBase(o) $COPY$ $LOCALSET$ {}
    
    void join(const $IKERNEL$ & o) {
        $JOINER$
    }
@END
);

const string TmpRunOMP = STR(
void run() {
@IF(IJK)
    const int _maxX = maxX; 
    const int _maxY = maxY;
    if (maxZ > 1) {
        $PRAGMA$ omp parallel $NL$
        {
            $OMP_DIRECTIVE$
            for (int k=minZ; k < maxZ; k++)
            for (int j=$BND$; j < _maxY; j++)
            for (int i=$BND$; i < _maxX; i++)
               op(i,j,k,$CALL$);
           $OMP_POST$
        }
    } else {
        const int k=0;
        $PRAGMA$ omp parallel $NL$
        {
            $OMP_DIRECTIVE$
            for (int j=$BND$; j < _maxY; j++)
            for (int i=$BND$; i < _maxX; i++)
                op(i,j,k,$CALL$);
            $OMP_POST$
        }
    }
@ELSE
    const int _sz = size;
    $PRAGMA$ omp parallel $NL$
    { 
        $OMP_DIRECTIVE$ 
        for (int i=0; i < _sz; i++)
            op(i,$CALL$);
        $OMP_POST$
    }
@END
}
);

const string TmpOMPDirective = STR (
this->threadId = omp_get_thread_num(); 
this->threadNum = omp_get_num_threads();
@IF(REDUCE)
    $OMP_PRE$
    $PRAGMA$ omp for nowait $NL$
@ELSE
    $PRAGMA$ omp for $NL$
@END
);


#define kernelAssert(x,msg) if(!(x)){errMsg(block.line0,string("KERNEL: ") + msg);}

void processKernel(const Block& block, const string& code, Sink& sink) {
    const Function& kernel = block.func;
    
    if (gDocMode) {
        sink.inplace << "//! \\ingroup Kernels\n" << block.func.minimal << "{}\n";
        return;
    }

    // process options
    bool idxMode = false, reduce = false, pts = false;
    bool hasLocals = !block.locals.empty(), hasRet = kernel.returnType.name != "void";
    string bnd = "0", reduceOp="";

    MType mtType = gMTType;
    for (size_t i=0; i<block.options.size(); i++) {
        const string& opt = block.options[i].name;
        if (opt == "ijk") 
            idxMode = false;
        else if (opt == "index" || opt == "idx")
            idxMode = true;
        else if (opt == "st" || opt == "single")
            mtType = MTNone;
        else if (opt == "pts" || opt == "particle" || opt == "points")
            pts = true;
        else if (opt == "bnd")
            bnd = block.options[i].value;
        else if (opt == "reduce") {
            reduce = true;
            reduceOp = block.options[i].value;
            if (!(reduceOp == "+" || reduceOp == "-" || reduceOp == "*" ||
                  reduceOp == "/" || reduceOp == "min" || reduceOp == "max"))
                errMsg(block.line0, "invalid 'reduce' operator. Expected reduce= +|-|*|/|min|max");
        } else
            errMsg(block.line0, "illegal kernel option '"+ opt +
                                "' Supported options are: 'ijk', 'idx', 'bnd=x', 'reduce=x', 'st', 'pts'");
    }
    
    // point out illegal paramter combinations
    kernelAssert (bnd == "0" || !idxMode, "can't combine index mode with bounds iteration.");    
    kernelAssert (!pts || (!idxMode && bnd == "0" ), 
        "KERNEL(opt): Modes 'ijk', 'idx' and 'bnd' can't be applied to particle kernels.");

    // check type consistency of first 'returns' with return type
    if (hasRet && kernel.returnType.name != "void") {
        kernelAssert(block.locals.size() == 1, "multiple returns statement only work for 'void' kernels");
        const Type& rt = block.locals[0].type;
        kernelAssert(rt == kernel.returnType, "return type does not match type in first 'returns' statement");
    }
    kernelAssert(kernel.returnType.name == "void" || hasRet, 
        "return argument specified without matching 'returns' initializer");
    
    // figure out basegrid
    string baseGrid;
    for (int i=0; i<kernel.arguments.size(); i++) {
        const string& type = kernel.arguments[i].type.name;
        bool isGrid = type.find("Grid") != string::npos;
        if (isGrid || pts) { 
            baseGrid = kernel.arguments[i].name;
            if (isGrid && !kernel.arguments[i].type.isPointer)
                baseGrid = "&"+baseGrid;
            break;
        }
    }
    kernelAssert(!baseGrid.empty(), ": use at least one grid to call the kernel.");

    // build accesors
    stringstream accessors;
    for (int i=0; i<kernel.arguments.size(); i++) {
        stringstream num; num << i;
        const string table[] = { "TYPE", kernel.arguments[i].type.build(true),
                                 "NAME", kernel.arguments[i].name, 
                                 "IDX", num.str(), 
                                 "@end" };
        accessors << replaceSet(TmpAccessor, table);
    }

    // build locals, and reduce joiners
    stringstream joiner, preReduce, postReduce;
    for (int i=0; i<block.locals.size(); i++) {
        const string& name = block.locals[i].name;
        const string type = block.locals[i].type.build();
        const string& value = block.locals[i].value;

        preReduce << type << " _L_" << name << " = " << value << ";";
        if (reduceOp == "min" || reduceOp == "max") {
            joiner << name << " = " << reduceOp << "(" << name << ",o." << name << "); ";
            postReduce << name << " = " << reduceOp << "(" << name << ",_L_" << name << "); ";
        } else {
            joiner << name << " " << reduceOp << "= o." << name << "; ";
            postReduce << name << " " << reduceOp << "= _L_" << name << "; ";
        }         
    }
    const string ompPost = reduce ? "\n#pragma omp critical\n{"+postReduce.str()+"}":"";
    bool doubleKernel = mtType == MTTBB && hasRet && !reduce;
    
    const string table[] = { "IDX", idxMode ? "Y":"",
                             "PTS", pts ? "Y":"",
                             "IJK", (!pts && !idxMode) ? "Y":"",
                             "REDUCE", reduce ? "Y":"",
                             "TEMPLATE", kernel.isTemplated() ? "template "+kernel.templateTypes.minimal : "",
                             "TPL", kernel.isTemplated() ? "<"+kernel.templateTypes.names()+">" : "",
                             "KERNEL", kernel.name,
                             "IKERNEL", (doubleKernel ? "_":"") + kernel.name,
                             "ARGS", kernel.arguments.listText,
                             "LOCALARG", block.locals.full(true),
                             "BASE", baseGrid,
                             "INIT", kernel.arguments.copier("", false),
                             "LOCALINIT", block.locals.copier("", false),
                             "LOCALSET", block.locals.copier("", true),
                             "COPY", kernel.arguments.copier("o.", false),
                             "MEMBERS", kernel.arguments.createMembers(),
                             "LOCALS", block.locals.createMembers(false),
                             "LOCALS_REF", block.locals.createMembers(true),
                             "ACCESSORS", accessors.str(),
                             "CONST", (!reduce && mtType==MTTBB) ? "const" : "",
                             "CODE", code,
                             "RET_TYPE", hasRet ? block.locals[0].type.minimal : "",
                             "RET_NAME", hasRet ? block.locals[0].name : "",
                             "BND", bnd,
                             "CALL", kernel.callString() + (hasLocals ? ","+block.locals.names() : ""),
                             "METHOD", reduce ? "reduce" : "for",
                             "PRAGMA", "\n#pragma",
                             "NL", "\n",
                             "COMMA", ",",
                             "JOINER", joiner.str(),
                             "OMP_PRE", preReduce.str(),
                             "OMP_POST", ompPost,
                             "@end" };

    // generate kernel
    string templ = doubleKernel ? TmpDoubleKernel : TmpSingleKernel;
    if (mtType == MTNone)
        replaceAll(templ, "$RUN$", TmpRunSimple);
    else if (mtType == MTTBB)
        replaceAll(templ, "$RUN$", TmpRunTBB);
    else if (mtType == MTOpenMP) {
        string ompTempl = TmpRunOMP;
        replaceAll(ompTempl, "$OMP_DIRECTIVE$", TmpOMPDirective);
        replaceAll(templ, "$RUN$", ompTempl);
    }

    // synthesize code
    sink.inplace << block.linebreaks() << replaceSet(templ, table);

    // adjust lines after OMP block
    if (mtType == MTOpenMP)
        sink.inplace << "\n#line " << block.line1 << " \"" << sink.infile << "\"\n" << endl;
}
