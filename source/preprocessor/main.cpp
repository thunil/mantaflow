/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor Main
 *
 ******************************************************************************/

#include <string>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include "prep.h"

using namespace std;

string gFilename;
string gRegText;
bool gDebugMode = true;
bool gDocMode;
bool gIsHeader;
MType gMTType = MTNone;


void usage() {
    cerr << "preprocessor error: Unknown parameters." << endl;
    cerr << "  Usage : prep generate <dbg_mode> <mt_type> <inputdir> <inputfile> <outputfile>" << endl;
    cerr << "     or : prep docgen <dbg_mode> <mt_type> <inputdir> <inputfile> <outputfile>" << endl;
    cerr << "     or : prep merge <outfile> <infiles...>" << endl;
    exit(1);
}

void doMerge(int argc, char* argv[]) {    
    if (argc < 4) usage();
    
    string merged = "";
    const string outfile(argv[2]);
    
    for (int i=3; i<argc; i++) {
        merged += "@filename " + string(argv[i]) + "\n";        
        merged += readFile(argv[i]) + "\n";
    }
    writeFile(outfile, generateMerge(merged));    
}

void doGenerate(int argc, char* argv[], bool docs) {
    gDocMode = docs;
    gDebugMode = false;
    gMTType    = MTNone;
    if (argc != 7) usage();
    
    // set constants
    const string indir(argv[4]), infile(argv[5]), outfile(argv[6]);
    gFilename = infile;
    gIsHeader = gFilename[gFilename.size()-2] == '.' && gFilename[gFilename.size()-1] == 'h';
    gRegText = "";
    
    // TP : only enable in cmake's PREP_DEBUG mode (passed via cmd line option dbg_mode)
    gDebugMode = atoi(argv[2]) != 0;
    if (!strcmp(argv[3],"TBB")) gMTType = MTTBB;
    if (!strcmp(argv[3],"OPENMP")) gMTType = MTOpenMP;
    
    // load complete file into buffer    
    string text = readFile(indir+infile);
    if (text.empty()) {
        cerr << "preprocessor error: Can't read file '" << infile << "'" << endl;
        exit(1);
    }
    // pad text for easier lexing lookups
    text += "\n\n\n";
    
    // process text
    if (infile.size() > 3 && !infile.compare(infile.size()-3, 3, ".py")) {
        // python file, only registering
        replaceAll(text, "\n", "\\n");
        replaceAll(text, "\r", "\\r");
        replaceAll(text, "\t", "\\t");
        replaceAll(text, "\"", "<qtm>"); // split into two ops to avoid interference
        replaceAll(text, "<qtm>", "\\\"");
        gRegText += "PbWrapperRegistry::instance().addPythonCode(\"" + gFilename + "\", \"" + text + "\");\n";
        
        writeFile(outfile, "");
    } 
    else {    
        // C++ files, normal processing
        Sink sink(outfile);
        if (gDocMode) {
            sink.inplace << "/*! \\file " + infile + " */\n";
        } else {
            string header = "\n\n\n\n\n// DO NOT EDIT !\n";
			header += "// This file is generated using the MantaFlow preprocessor (prep generate).";
			header += "\n\n\n\n\n";
            sink.inplace << header;
            sink.headerExt = header + "#include \"" + infile + "\"\nnamespace Manta ";
			// Note - PREP_DEBUG mode is meant to debug problems in the generated files, so errors should not be redirected to the sources
			// so dont always add this line...
	        if (!gDebugMode) 
				sink.inplace << "#line 1 \"" << indir << infile << "\"\n";
        }
        processText(text, 1, sink, 0);
		sink.write();
    }
}


int main(int argc, char* argv[]) {
    // command line options
    if (argc < 2) usage();

    // use merger
    if (!strcmp(argv[1],"merge")) 
        doMerge(argc, argv);
    else if (!strcmp(argv[1],"generate")) 
        doGenerate(argc, argv, false);
    else if (!strcmp(argv[1],"docgen"))
        doGenerate(argc, argv, true);
    else 
        usage();
    
    return 0;
}
 
