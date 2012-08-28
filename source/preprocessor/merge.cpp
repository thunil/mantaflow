/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor merge file gen
 *
 ******************************************************************************/

#include <map>
#include <iostream>
#include <cstdlib>
#include "prep.h"
using namespace std;

string generateMerge(const string& text) {
    // split lines 
    vector<string> lines, replLines;
    lines.push_back("");
    for (size_t i=0; i<text.size(); i++) {
        if (text[i] == '\n') {
            lines.push_back("");            
        } else
            *lines.rbegin() += text[i];
    }
    
    // collect templates
    map<string,vector<string> > templates;
    for(vector<string>::iterator it = lines.begin(); it != lines.end(); it++) {
        if (!it->compare(0,10,"@template ")) {
            size_t p = it->find(' ', 10);
            string tclass = it->substr(10, p-10);
            string expr = it->substr(p+1);
            templates[tclass].push_back(expr);
        }        
    }
    
    // instantiate templates
    for(vector<string>::iterator it = lines.begin(); it != lines.end(); it++) {
        if (it->empty()) continue;
        if (!it->compare(0,10,"@instance ")) {
            size_t p = it->find(' ', 10);
            string tclass = it->substr(10, p-10);
            string targ = it->substr(p+1);
            p = targ.find(' ');
            string aliasn = targ.substr(p+1);
            targ = targ.substr(0,p);
            vector<string>& exprs = templates[tclass];
            for(vector<string>::iterator ie = exprs.begin(); ie != exprs.end(); ie++) {
                string ex = *ie;
                replaceAll(ex, "$", targ);
                replaceAll(ex, "@", aliasn);
                replLines.push_back(ex);
            }
        } else if ((*it)[0] != '@') {
            replLines.push_back(*it);
        }
    }

    //sort in heaps
    string declares, calls, includes, converters;
    for(vector<string>::iterator it = replLines.begin(); it != replLines.end(); it++) {
        if (it->size() < 2) continue;
        if (!it->compare(0, 8, "#include")) {
            if (includes.find(*it) == string::npos)
                includes += (*it)+"\n";
        } else if (!it->compare(0, 10, "template<>")) {
            converters += (*it) + "\n";
        } else if (!it->compare(0, 9, "PbWrapper"))
            calls += "\t\t" + (*it) + "\n";
        else
            declares += (*it) + "\n";
    }
    
    // write header
    string newText = "\n\n\n\n\n// This file is generated using the MantaFlow preprocessor merge stage (prep merge). Do not edit.\n\n\n\n\n";
    newText += "#include \"pclass.h\"\n";
    newText += includes + "\n";
    newText += "namespace Manta {\n";
    newText += converters + declares + "\n";
    newText += "struct _pbRegister {\n";
    newText += "\t_pbRegister(int dummy) {\n";
    newText += calls;
    newText += "\t}\n};\n";
    newText += "static const Manta::_pbRegister _doRegister(12);\n";
    newText += "} //namespace\n";
    
    return newText;
}