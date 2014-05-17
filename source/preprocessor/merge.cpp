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
#include <sstream>
#include <cstdlib>
#include "prep.h"
using namespace std;

struct Chain {
    Chain() {};
    Chain(const string& tpl, const string& target, const string& targetTpl) : 
        tpl(tpl), target(target), targetTpl(targetTpl) {}
    string tpl, target, targetTpl;
};

struct Request {
    Request(const string& c, const string& t) : cls(c), tpl(t), base("") {}
    string cls, tpl, base;
};

struct RegFile {
    RegFile(const string& name) : filename(name),active(false) {}

    string filename;
    ostringstream out, header;
    vector<Request> req;
    bool active;
};

struct ClassInfo {
    map<string,bool> tplDone;
    vector<string> snippets;
};

static map<string, ClassInfo> classes;
static map<string, Chain> chains;
static vector<RegFile *> regFiles;

// compare template arguments
string mapArgs(const string& inst, const string& match, const string& target) {
    vector<string> inArg = split(inst,',');
    vector<string> maArg = split(match,',');
    vector<string> taArg = split(target,',');
    vector<string> destArg = taArg;
    
    for (int i=0; i<maArg.size(); i++) {
        for (int j=0; j<taArg.size(); j++) {
            if (maArg[i] == taArg[j]) {
                destArg[j] = inArg[i];
            }
        }
    }
    stringstream s;
    for (int i=0; i<destArg.size(); i++) {
        s << destArg[i];
        if (i != destArg.size()-1) s << ',';
    }
    return s.str();
}

void resolveChains(RegFile& file) {
    for (int i=0; i<file.req.size(); i++) {
        Request& req = file.req[i];
        map<string, Chain>::iterator it = chains.find(req.cls);
        if (it != chains.end()) {
            Chain& chain = it->second;
            string tpl = mapArgs(req.tpl, chain.tpl, chain.targetTpl);
            req.base = tpl;
            file.req.push_back(Request(chain.target, tpl));
        }
    }
}

void resolveRequests(RegFile& file) {
    // sort request by class
    map<string, vector<Request*> > sortedReqs;
    for (int i=0; i<file.req.size(); i++) {
        Request& req = file.req[i];
        ClassInfo& info = classes[req.cls];
        if (!info.tplDone[req.tpl]) {
            info.tplDone[req.tpl] = true;
            sortedReqs[req.cls].push_back(&req);
        }
    }

    // process requests
    int idx = 0;
    for(map<string,vector<Request*> >::iterator it = sortedReqs.begin(); it != sortedReqs.end(); ++it) {
        ClassInfo& info = classes[it->first];
        file.out << "#ifdef _C_" << it->first << '\n';
        for (int i=0; i<it->second.size(); i++) {
            Request& req = *(it->second[i]);
            for (int j=0; j<info.snippets.size(); j++) {
                stringstream idxStr;
                idxStr << idx;
                const string table[] = {"CT", req.tpl, "BT", req.base, "IDX", idxStr.str(), "@end"};
                file.out << replaceSet(info.snippets[j], table) << '\n';
                idx++;
            }
        }
        file.out << "#endif\n";
        file.active = true;
    }
}

// create data structure from regfiles
void parseLine(const string& line, RegFile& file) {
    if (line.empty()) return;

    vector<string> parts = split(line,'^');
    string cls = parts[0].substr(1);
    
    if (line[0] == '+')
        classes[cls].snippets.push_back(parts[1]);
    else if (line[0] == '>') 
        file.req.push_back(Request(cls,parts[1]));
    else if (line[0] == '@')
        chains[cls] = Chain(parts[1],parts[2],parts[3]);
    else if (line[0] == '#')
        file.header << line << '\n';
    else {
        file.out << line << '\n';
        file.active = true;
    }
}

void generateMerge(int num, char* files[]) {
    // parse files
    for (int i=0; i<num; i++) {
        regFiles.push_back(new RegFile(files[i]));

        string text = readFile(files[i]);
        replaceAll(text,"\r","");
        vector<string> lines = split(text,'\n');

        for (int j=0; j<lines.size(); j++)
            parseLine(lines[j], *regFiles.back());
    }

    // process and save files
    for (int i=0; i<num; i++) {
        resolveChains(*regFiles[i]);
        resolveRequests(*regFiles[i]);

        string text = "";
        if (regFiles[i]->active) {
            text  = "\n\n\n\n\n// DO NOT EDIT !\n";
            text += "// This file is generated using the MantaFlow preprocessor (prep link).";
            text += "\n\n\n\n\n";
            text += regFiles[i]->header.str();
            text += "namespace Manta {\n";
            text += regFiles[i]->out.str();
            text += "}";
        }
        string filename = regFiles[i]->filename + ".cpp";
        // only write if content is different
        if (!fileExists(filename) || readFile(filename) != text)
            writeFile(filename, text);
        delete regFiles[i];
    }
}
