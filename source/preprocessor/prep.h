/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor Declarations
 *
 ******************************************************************************/

#ifndef _PREP_H
#define _PREP_H

#include "util.h"
#include "code.h"
#include "tokenize.h"

// from main.cpp
enum MType { MTNone = 0, MTTBB, MTOpenMP};
extern std::string gFilename;
extern bool gDebugMode;
extern MType gMTType;
extern bool gDocMode;

// functions from merge.cpp
void generateMerge(int num, char* files[]);

// functions from tokenize.cpp
void processText(const std::string& text, int baseline, Sink& sink, const Class *parent);

// functions from parse.cpp
void parseBlock(const std::string& kw, const std::vector<Token>& tokens, const Class *parent, Sink& sink);

// functions from codegen_XXX.cpp
void processKernel(const Block& block, const std::string& code, Sink& sink);
void processPythonFunction(const Block& block, const std::string& code, Sink& sink);
void processPythonVariable(const Block& block, Sink& sink);
void processPythonClass(const Block& block, const std::string& code, Sink& sink);
void processPythonInstantiation(const Block& block, const Type& aliasType, const std::string& aliasName, Sink& sink);

#endif // _PREP_H
