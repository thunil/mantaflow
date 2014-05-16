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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <stack>
#include "prep.h"
using namespace std;

void errMsg(int line, const string& text) {
    cerr << gFilename << ":" << line << ": error: " << text << endl;
    exit(1);
}

void debMsgHelper(int line, const string& text) {
    cout << gFilename << ":" << line << ": info: " << text << endl;
}

bool fileExists(const string& name) {
    ifstream t(name.c_str());
    if (!t.good()) return false;
    t.close();
    return true;
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

Sink::Sink(const string& outfile):
    filename(outfile)
{
    isHeader = outfile.compare(outfile.size()-3, 3, ".py") == 0 ||
               outfile.compare(outfile.size()-2, 2, ".h") == 0;
}

void Sink::write() {
    writeFile(filename, inplace.str());
    if (isHeader && !gDocMode) {
        writeFile(filename + ".reg", link.str());
    }
}

vector<string> split(const string& text, char sep) {
    vector<string> bins;
    string cur;
    for (int i=0; i<text.size(); i++) {
        if (text[i]==sep) {
            bins.push_back(cur);
            cur = "";
        } else
            cur += text[i];
    }
    bins.push_back(cur);
    return bins;
}

void replaceAll(string& source, string const& find, string const& replace)
{
    for(string::size_type i = 0; (i = source.find(find, i)) != std::string::npos;)
    {
        source.replace(i, find.length(), replace);
        i += replace.length() - find.length() + 1;
    }
}

// Helpers for replaceSet
static string getBracketArg(const string& a, int &pos) {
    string ret="";
    pos++;
    for (;pos<a.size(); pos++) {
        if (a[pos]!='(' && a[pos]!=' ' && a[pos]!='$') break;
    }
    for (; pos<a.size(); pos++) {
        if (a[pos]==')') return ret;
        ret += a[pos];
    }
    return "";
}

inline bool compareKW(const string& a, int& pos, const string& kw) {
    if (a.compare(pos+1,kw.size(),kw) == 0) {
        pos += kw.size() ;
        return true;
    }
    return false;
}

string replaceSet(const string& templ, const string table[]) {
    vector<string> key, value;
    for (int i=0;table[i] != "@end";i+=2) {
        key.push_back(table[i]);
        value.push_back(table[i+1]);
    }
    stack<bool> conditionStack;
    conditionStack.push(true);
    stringstream s;
    for (int i=0; i<templ.size(); i++) {
        char c = templ[i];
        if (c=='@') {
            if (compareKW(templ,i,"IF")) {
                string cond = getBracketArg(templ,i);
                vector<string>::iterator it = find(key.begin(),key.end(),cond);
                bool res = (it != key.end()) && (!value[it-key.begin()].empty());
                conditionStack.push(res);
            } else if (compareKW(templ,i,"END")) {
                conditionStack.pop();
            } else if (compareKW(templ,i,"ELSE")) {
                conditionStack.top() = !conditionStack.top();
            }
            continue;
        }
        if (!conditionStack.top()) continue;
        if (c=='$') {
            int oldi=i;
            for (int j=0; j<key.size(); j++) {
                if (compareKW(templ,i,key[j])) {
                    s << value[j];
                    break;
                }
            }
            if (oldi != i) continue;
        }
        s << templ[i];
    }
    return s.str();
}