/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Main file
 *
 ******************************************************************************/

#include <stdio.h>
#include "pythonInclude.h"
#include "pclass.h"
#include "pconvert.h"
#include "general.h"

namespace Manta {
    extern void guiMain(int argc, char* argv[]);
    extern void guiWaitFinish();
}

using namespace std;
using namespace Manta;

//*****************************************************************************
// main...

void runScript(vector<string>& args) {
    string filename = args[0];
    
    // Initialize extension classes and wrappers
    srand(0);
    PbWrapperRegistry::instance().construct(filename);
    Py_Initialize();
    PbWrapperRegistry::instance().runPreInit(args);
        
    // Try to load python script
    FILE* fp = fopen(filename.c_str(),"r");
    if (fp == NULL) {
        debMsg("Cannot open '" << filename << "'", 0);
        Py_Finalize();
        return;
    }
    
    // Run the python script file
    debMsg("Loading script '" << filename << "'", 0);
    PyRun_SimpleFileEx(fp, filename.c_str(), 1);
    
    debMsg("Script finished.", 0);
#ifdef GUI
    guiWaitFinish();
#endif

    // finalize
    Py_Finalize();
    PbWrapperRegistry::instance().cleanup();    
}

int main(int argc,char* argv[]) {
    if (argc<=1) {
        cerr << "Usage : Syntax is 'ddf <config.py>'" << endl;  
        return 1;
    }

#ifdef GUI
    guiMain(argc, argv);    
#else
    vector<string> args;
    for (int i=1; i<argc; i++) args.push_back(argv[i]);
    runScript(args);
#endif        		
	return 0;
}
