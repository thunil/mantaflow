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
bool gUseMT;
bool gDocMode;

void errMsg(int line, const string& text) {
    cerr << gFilename << ":" << line << ": error: " << text << endl;
    exit(1);
}

string readFile(const string& name) {
    ifstream t(name.c_str());
    if (!t.good()) return "";
    
    string str;
    
    t.seekg(0, ios::end);   
    str.reserve((unsigned)t.tellg());
    t.seekg(0, ios::beg);

    str.assign((istreambuf_iterator<char>(t)),
                istreambuf_iterator<char>());
    
    t.close();
    return str;
}

void writeFile(const string& name, const string& text) {
    ofstream ofs(name.c_str(), ios::binary | ios::out);
    if (!ofs.good()) {
        cerr << "preprocessor error: Can't write to file '" << name << "'" << endl;
        exit(1);
    }
    ofs << text;
    ofs.close();
}

void replaceAll(string& text, const string& pattern, const string& repl) {
    for(;;) {
        const size_t pos = text.find(pattern);
        if (pos == string::npos) return;
        text.replace(pos, pattern.size(), repl);
    }
    return;
}

void usage() {
    cerr << "preprocessor error: Unknown parameters." << endl;
    cerr << "  Usage : prep generate <dbg_mode> <usemt> <inputfile> <outputfile>" << endl;
    cerr << "     or : prep docgen <inputfile> <outputfile>" << endl;
    cerr << "     or : prep merge <outfile> <infiles...>" << endl;
    exit(1);
}

void doMerge(int argc, char* argv[]) {    
    if (argc < 4) usage();
    
    string merged = "";
    const string outfile(argv[2]);
    
    for (int i=3; i<argc; i++) {
        merged += readFile(argv[i]) + "\n";
    }
    writeFile(outfile, generateMerge(merged));    
}

void doGenerate(int argc, char* argv[], bool docs) {
    gDocMode = docs;
    if (docs && argc != 4) usage();
    if (!docs && argc != 6) usage();
    
    // set constants
    const string infile (argv[docs ? 2 : 4]), outfile(argv[docs ? 3 : 5]);
    gFilename = infile;
    gRegText = "";
    gParent = "";
    // TP : only enable in cmake's PREP_DEBUG mode (passed via cmd line option dbg_mode)
	gDebugMode = false; if (!docs) gDebugMode = atoi(argv[2]) != 0;
    gUseMT     = false; if (!docs) gUseMT     = atoi(argv[3]) != 0;

	bool addLineWarnings = false;
#	ifdef _WIN32
	addLineWarnings = true;
	// NT_DEBUG ? debug on? 
	//gDebugMode = true;
#	endif
    
    // load complete file into buffer    
    string text = readFile(infile);
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
        string newText;
        if (docs) {
            string fn = outfile;
            //while (fn.find('/') != string::npos) fn = fn.substr(fn.find('/')+1);
            newText += "/*! \\file " + fn + " */\n";
        } else {
            newText += "\n\n\n\n\n// This file is generated using the MantaFlow preprocessor (prep generate). Do not edit.\n\n\n\n\n";
            
			// NT_DEBUG , always add? 
			// TP: No -- PREP_DEBUG mode is meant to debug problems in the generated files, so erros should not be redirected to the sources
	        if (!gDebugMode) 
				newText += "#line 1 \"" + infile + "\"\n";
        }
        newText += processText(text, 1);


		// windows special - mark generated files with warnings to prevent editing the wrong file in MSVC
		if( (!docs) && addLineWarnings) 
		{
			// position counters
			int lineCntOrg = 0;   // line count of original file
			int lineCnt    = -32; // interval line counter, skip the first few lines...
			int lastPos    = 0;   // get last string chunk
			int lastLength = 1;   // string chunk length
			// how often to add warning?
			const int lineInterval = 10;

			// construct output
			std::ostringstream modifiedNew;
			for(int c=0; c < (int)newText.length(); c++,lastLength++) {
				bool skip = false;
				if(newText[c] == '\n') {
					lineCntOrg++;

					// check for define coming up...
					if ( (c < (int)newText.length()-1) &&  ( newText[c+1] == '#' ) ) {
						skip = true;
					} else {
						// extended line?
						if( (c > 0) && ( newText[c-1] == '\\' ) ) {
							// dont increase
							skip = true;
						} else {
							// "regular" newline
							lineCnt++;
						}
					}

					if(!skip) {
						if( ( (lineCnt>0) && ((lineCnt%lineInterval) == (lineInterval-1)) )  ||
							(c==(int)newText.length()-1) ) 
						{
							modifiedNew << newText.substr( lastPos, lastLength );
							// add warning
							modifiedNew << "\n// AUTO-GENERATED FILE, DO NOT EDIT!!! \n\n";
							// make sure we find the original file & line, subtract length of added header at top (10 lines)
							modifiedNew << "#line " <<(lineCntOrg-10)<< " \"" + infile + "\" "; // no newline here! is in the file anyway...

							// start over
							lastPos = c;
							lastLength = 0;
						}
					}
				} // newline
			}
			newText = modifiedNew.str();
		}        

        // write down complete file from buffer
        writeFile(outfile, newText);        
    }
        
    // write reg functions into file
    if (!docs)
        writeFile(outfile + ".reg", gRegText);
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
 