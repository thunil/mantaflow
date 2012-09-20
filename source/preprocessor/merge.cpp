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
    string lookupClass;
    string lookupArgs;
    string chainClass;
    string chainArgs;
    string complete;
};

struct ChainedInstance {
    string tclass, targ, aliasn;
    string chain;
};

vector<string> splitString(string s, char sep) {
    string l="";
    vector<string> tokens;
    for (size_t i=0; i<s.size(); i++) {
        if (s[i]==sep) {
            tokens.push_back(l);
            l="";
        } else
            l+=s[i];
    }
    if (!l.empty()) tokens.push_back(l);
    return tokens;
}

void applyChain(vector<Chain>& chains, vector<ChainedInstance>& instances, const string& id, const string& aliasnDummy) {
    string classname = id.substr(0,id.find('@'));
    string templArgs = id.substr(id.find('@')+1);
    
    // match to chain
    for (size_t i=0; i<chains.size(); i++) {
        if (chains[i].lookupClass == classname) {
            // class found, create argument replacement table
            vector<string> argName = splitString(chains[i].lookupArgs,',');
            vector<string> argDef = splitString(templArgs,',');
            if (argName.size()!=argDef.size())
                errMsg(0,"internal merge error, chain replace inconsistency.");
            
            string arg=chains[i].chainArgs;
            for (size_t j=0; j<argName.size(); j++)
                replaceAll(arg,argName[j],argDef[j]);
            
            // create chained instance                
            // TODO: do this recursively 
            ChainedInstance inst;
            inst.tclass = chains[i].chainClass;
            inst.targ = arg;
            inst.aliasn = "_" + chains[i].chainClass + "_" + arg;
            inst.chain = chains[i].complete;
            instances.push_back(inst);
        }
    }
}
    
void applyInstance(map<string,vector<string> >& templates, vector<string>& replLines, map<string,bool>& didInstance, const string& tclass, const string& targ, const string& aliasn ) {
    // avoid duplicates   
    string id = tclass + "@" + targ;
    if (didInstance.find(id) == didInstance.end()) 
        didInstance[id] = true;
    else
        return;
            
    // apply instance
    if (templates.find(tclass) == templates.end()) 
        errMsg(0, "template instancing of '" + tclass + "<"+ targ + ">' failed: class '" + tclass + "' not found");
    
    // generate converters
    replLines.push_back(createConverters(tclass + "<" + targ + ">", "", " ", "\n"));

    // generate template register calls
    vector<string>& exprs = templates[tclass];
    vector<string> arglist = splitString(targ,',');                
    for(vector<string>::iterator ie = exprs.begin(); ie != exprs.end(); ie++) {
        string ex = *ie;
        replaceAll(ex, "$$", targ);
        replaceAll(ex, "@", aliasn);
        // indivual arg replacement
        for(size_t i=0; i<arglist.size(); i++) {
            stringstream s;
            s << "$" << i << "$";
            if (ex.find(s.str()) == ex.npos) break;
            replaceAll(ex, s.str(), arglist[i]);
        }
        replLines.push_back(ex);
    }
}

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
    
    // collect and instantiate templates (two-pass)
    map<string,vector<string> > templates;
    map<string, string> instances;
    vector<Chain> chains;
    for(vector<string>::iterator it = lines.begin(); it != lines.end(); it++) {
        if (it->empty()) continue;
        if (!it->compare(0,10,"@instance ")) {
            // extract args
            size_t p = it->find(' ', 10);
            string tclass = it->substr(10, p-10);
            string targ = it->substr(p+1);
            p = targ.find(' ');
            string aliasn = targ.substr(p+1);
            targ = targ.substr(0,p);
            
            // is instance already registered ?
            string id = tclass + "@" + targ;
            if (instances.find(id) == instances.end())
                instances.insert(pair<string,string>(id,aliasn));
            else {
                // update alias name ?
                if (aliasn[0] != '_' && instances[id][0] == '_')
                    instances[id] = aliasn;                
            }
        } else if (!it->compare(0,10,"@template ")) {
            size_t p = it->find(' ', 10);
            string tclass = it->substr(10, p-10);
            string expr = it->substr(p+1);
            templates[tclass].push_back(expr);
        } else if (!it->compare(0,7,"@chain ")) {
            // extract args
            size_t p0 = it->find(' ', 7);
            size_t p1 = it->find(' ', p0+1);
            size_t p2 = it->find(' ', p1+1);
            Chain c;
            c.complete = *it;
            c.lookupClass = it->substr(7, p0-7);
            c.lookupArgs = it->substr(p0+1, p1-p0-1);
            c.chainClass = it->substr(p1+1, p2-p1-1);
            c.chainArgs = it->substr(p2+1);
            chains.push_back(c);
        }        
    }
    
    // apply chains
    vector<ChainedInstance> chainedInstances;
    for (map<string,string>::iterator it = instances.begin(); it != instances.end(); ++it)
        applyChain(chains, chainedInstances, it->first, it->second);
    
    map<string, bool> didInstance;
    for(vector<string>::iterator it = lines.begin(); it != lines.end(); it++) {
        if (it->empty()) continue;
        if (!it->compare(0,7,"@chain ")) {
            for (vector<ChainedInstance>::iterator c=chainedInstances.begin(); c!=chainedInstances.end(); ++c) {
                if (c->chain == *it)
                    applyInstance(templates, replLines, didInstance, c->tclass, c->targ, c->aliasn);
            }
        } else if (!it->compare(0,10,"@instance ")) {
            // extract args
            size_t p = it->find(' ', 10);
            string tclass = it->substr(10, p-10);
            string targ = it->substr(p+1);
            targ = targ.substr(0,targ.find(' '));
            string aliasn = instances[tclass + "@" + targ];
            
            applyInstance(templates, replLines, didInstance, tclass, targ, aliasn);
        } else if (!it->compare(0,10,"@filename ")) {
            gFilename = it->substr(10);
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