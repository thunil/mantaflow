/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor Parsing
 *
 ******************************************************************************/

#include <iostream>
#include <cstdlib>
#include "prep.h"

using namespace std;

// ignore whitspaces, return number of lines
int consumeWS(const vector<Token>& tokens, size_t& index) {
    int linebreaks = 0;
    while (tokens[index].type == TkWhitespace || tokens[index].type == TkComment) {
        if (tokens[index].text.find('\n') != string::npos)
            linebreaks++;
        index++;
        if (index >= tokens.size())
            errMsg(-1, "Preprocessor ran out of tokens. This shouldn't happen.");
    }
    return linebreaks;
}

// parse argument list tokens of the form (int a = 12, bool* b, Grid<float>& grid)
vector<Argument> parseArgs(const vector<Token>& tokens, size_t& index, bool expectType, int &lb) {    
    vector<Argument> args;
    const int firstStage = expectType ? 0 : 4;
    
    // ignore trailing whitespaces, remove first bracket
    lb += consumeWS(tokens, index);
    if (tokens[index].type != TkBracketL) return args;
    
    // tokenizer already established that open brackets are closed eventually, so no need to check that
    Argument cur = Argument();    
    int stage = firstStage; // 0: type, 1: template open / modifiers / name 2: template, 3: template close 4: name / modifiers, 5: assignOp 6+: defaults
    int number = 0, opened=1;
    do {
        index++;    
        TokenType t = tokens[index].type;
        string text = tokens[index].text;
        int line = tokens[index].line;
        
        if (t==TkBracketL) opened++;
        if (t==TkBracketR) opened--;
        //cout << index << ":" << stage << "," << t << "," << text << " " << opened << endl;
        
        if (t == TkWhitespace || t == TkComment) {
            if (text[0]=='\n') lb++;
            else if (text[0]!='\r') cur.complete += text;            
            continue;
        }
        if ((opened <1 && t == TkBracketR) || (opened<=1 && t == TkComma)) { 
            // argument complete, add to list
            if (stage == 5 || stage > 6) {
                args.push_back(cur);
                cur = Argument();
                number++;
                cur.number = number;
            }
            else if (stage != firstStage)
                errMsg(line, "incomplete argument! expect argument of the form '[type] name [= default]'");            
            
            stage = firstStage; 
            continue;
        }
        cur.complete += text;
        switch (stage) {
            case 0:
                // type name
                if (t != TkDescriptor)
                    errMsg(line, "incomplete argument! expect argument of the form '[type] name [= default]'");
                if (text == "const") {
                    cur.isConst = true;                    
                } 
                else {
                    cur.type += text;
                    stage++;
                }
                break;
            case 1:
                if (t == TkColon && tokens[index+1].type == TkColon) { // type had namespace
                    cur.type += "::";
                    cur.complete += ":";
                    index++;
                    stage=0;
                    break;
                } // else : fallthough to name
            case 4:
                // template open / name / modifiers
                if (t == TkTBracketL && stage == 1) {
                    stage++;
                } else if (t == TkRef) {
                    cur.isRef = true;
                    stage = 4;
                } else if (t == TkPointer) {
                    cur.isPointer = true;
                    stage = 4;
                } else if (t == TkDescriptor) {
                    cur.name = text;
                    stage = 5;
                } else 
                    errMsg(line, "incomplete argument! expect argument of the form '[type] name [= default]'");
                break;
            case 2:
                // template
                if (t != TkDescriptor)
                    errMsg(line, "incomplete argument type! expect type<Template> [*|&]");
                cur.templ = text;
                stage++;
                break;
            case 3:
                // template close
                if (t != TkTBracketR)
                    errMsg(line, "incomplete argument type! expect type<Template> [*|&]");
                stage++;
                break;
            case 5:
                // assign op
                if (t != TkAssign)
                    errMsg(line, "incomplete argument! expect argument of the form '[type] name [= default]'");
                stage++;
                break;
            default:
                // default values, stage 6+
                if (t == TkNone || t == TkDescriptor || t==TkColon || t==TkComma || t==TkBracketL || t==TkBracketR)
                    cur.value += text;
                else 
                    errMsg(line, "incomplete argument! expect argument of the form '[type] name [= default]'");
                stage++;
            break;
        }        
    } while(tokens[index].type != TkBracketR || opened>0);
    
    index++; // consume closing bracket
    
    // ignore WS
    lb += consumeWS(tokens, index);    
    return args;
}

// Parse syntax KEYWORD(opt1, opt2, ...) STATEMENTS [ {} or ; ]    
string parseBlock(const string& kw, const vector<Token>& tokens, int line) {    
    
    // parse keyword options
    Keyword key = checkKeyword(kw);
    size_t index = 0;
    int lb = 0;
    
    vector<Argument> options = parseArgs(tokens, index, false, lb);
    
    if (tokens[index].type != TkDescriptor)
        errMsg(line, "malformed preprocessor keyword block. Expected '" + kw + "(opt1, opt2, ...) STATEMENTS [{}|;]'");
    
    // keyword dependent processing. 
    // No need to check for index overflow, as long as test each type, since sequence will end in ; or {}
    if (key == KwPlugin) {
        string type = tokens[index++].text;        
            
        lb += consumeWS(tokens, index);
        if (tokens[index].type != TkDescriptor)
            errMsg(line, "Malformed PLUGIN, expected: PLUGIN [type|class] name(args...) { code }");
        string name = tokens[index++].text;        
        
        vector<Argument> args = parseArgs(tokens, index, true, lb);
        
        if (tokens[index].type != TkCodeBlock || index+1 != tokens.size())
            errMsg(line, "Malformed PLUGIN, expected: PLUGIN [class] name(args...) { code }");
        
        return processPlugin(lb, name, options, args, tokens[index].text, line, type);
    }
    else if (key == KwKernel) {
        string name = tokens[index++].text;
        string templ = "";
        
        // resolve template class or instantiation
        if (name == "template") {
            // template kernel
            if (tokens[index].type != TkTBracketL)
                errMsg(line, "syntax error. expected KERNEL(...) template<class T> [class] ...");            
            do {
                index++;
                lb += consumeWS(tokens, index); 
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "syntax error. expected KERNEL(...) template<class T> [class] ...");
                if (!templ.empty()) templ += ", ";
                templ += tokens[index++].text + " ";            
                lb += consumeWS(tokens, index);
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "syntax error. expected KERNEL(...) template<class T> [class] ...");
                templ += tokens[index++].text;
                lb += consumeWS(tokens, index);
            } while(tokens[index].type == TkComma);
            if (tokens[index++].type != TkTBracketR)
                errMsg(line, "syntax error. expected KERNEL(...) template<class T> [class] ...");
            lb += consumeWS(tokens, index);
            if (tokens[index].type != TkDescriptor)
                errMsg(line, "syntax error. expected KERNEL(...) template<class T> [class] ...");
            name = tokens[index++].text;
            lb += consumeWS(tokens, index);            
        }                
        
        bool isClass = name == "class" || name == "struct";        
        if (isClass) {
            lb += consumeWS(tokens, index);
            if (tokens[index].type != TkDescriptor)
                errMsg(line, "Malformed KERNEL, expected: KERNEL struct name(args...) { code }");
            
            name = tokens[index++].text;
        } 
        vector<Argument> args = parseArgs(tokens, index, true, lb);
        
        if (tokens[index].type != TkCodeBlock || index+1 != tokens.size())
            errMsg(line, "Malformed KERNEL, expected: KERNEL(opts...) name(args...) { code }");
        
        return processKernel(lb, name, options, templ, args, tokens[index].text, line, isClass);
    }
    else if (key == KwPython) {
        string templ = "";
        string type = tokens[index++].text; 
        lb += consumeWS(tokens, index);
        
        // resolve template class or instantiation
        if (type == "template") {
            if (tokens[index].type == TkTBracketL) {
                // template class
                index++;
                lb += consumeWS(tokens, index); 
                if (tokens[index].type != TkDescriptor || (tokens[index].text != "class" && tokens[index].text != "typename"))
                    errMsg(line, "syntax error. expected PYTHON template<class T> class ...");
                index++;
                lb += consumeWS(tokens, index);
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "syntax error. expected PYTHON template<class T> class ...");
                templ = tokens[index++].text;
                lb += consumeWS(tokens, index);
                if (tokens[index++].type != TkTBracketR)
                    errMsg(line, "syntax error. expected PYTHON template<class T> class ...");
                lb += consumeWS(tokens, index);
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "syntax error. expected PYTHON template<class T> class ...");
                type = tokens[index++].text;
                lb += consumeWS(tokens, index);
            }
            else if (tokens[index].type == TkDescriptor && tokens[index].text == "class") {
                // template instantiation
                index++;
                lb += consumeWS(tokens, index);
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "syntax error. expected PYTHON class name<T>;");
                string classname = tokens[index++].text;
                lb += consumeWS(tokens, index);
                if (tokens[index++].type != TkTBracketL)
                    errMsg(line, "syntax error. expected PYTHON class name<T>;");
                lb += consumeWS(tokens, index);
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "syntax error. expected PYTHON class name<T>;");
                templ = tokens[index++].text;
                lb += consumeWS(tokens, index);
                if (tokens[index++].type != TkTBracketR)
                    errMsg(line, "syntax error. expected PYTHON class name<T>;");
                lb += consumeWS(tokens, index);
                if (tokens[index].type != TkSemicolon || index+1 != tokens.size())
                    errMsg(line, "syntax error. expected PYTHON class name<T>;");
                processPythonInstantiation(lb, classname, templ, line);
                return "";
            }
            else 
                errMsg(line, "syntax error. expected PYTHON class name<T>; or PYTHON template<class T> class name...");                
        }
        
        // resolve extended function return type
        if (type == "const") {
            if (tokens[index].type != TkDescriptor)
                errMsg(line, "malformed preprocessor keyword block.");
            type += " " + tokens[index++].text;
        }
        while (tokens[index].type == TkPointer || tokens[index].type == TkRef) {
            type += tokens[index++].text;
            lb += consumeWS(tokens, index);
        }
        if (tokens[index].type == TkTBracketL) {
            index++;
            lb += consumeWS(tokens, index);
            if (tokens[index].type != TkDescriptor)
                errMsg(line, "PYTHON return type template malformed");
            type += "<"+tokens[index++].text+">";
            lb += consumeWS(tokens, index);
            if (tokens[index++].type != TkBracketR)
                errMsg(line, "PYTHON return type template malformed");            
            lb += consumeWS(tokens, index);            
        }
        while (tokens[index].type == TkPointer || tokens[index].type == TkRef) {
            type += tokens[index++].text;
            lb += consumeWS(tokens, index);
        }
                
        if (type == "class") {
            if (tokens[index].type != TkDescriptor)
                errMsg(line, "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");
            string name = tokens[index++].text;
            
            lb += consumeWS(tokens, index);
            if (tokens[index++].type != TkColon)
                errMsg(line, "PYTHON class '" + name + "' must publicly derive from PbClass (or a subclass)");
            
            lb += consumeWS(tokens, index);
            if (tokens[index].type != TkDescriptor || tokens[index++].text != "public")
                errMsg(line, "PYTHON class '" + name + "' must publicly derive from PbClass (or a subclass)");
            
            lb += consumeWS(tokens, index);
            if (tokens[index].type != TkDescriptor)
                errMsg(line, "PYTHON class '" + name + "' must publicly derive from PbClass (or a subclass)");
            string baseclass = tokens[index++].text;
            
            lb += consumeWS(tokens, index);
            if (tokens[index].type == TkTBracketL) {
                index++;
                lb += consumeWS(tokens, index);
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "malformed template argument, expected baseclass<T>");
                baseclass += "<" + tokens[index++].text + ">";
                
                lb += consumeWS(tokens, index);
                if (tokens[index++].type != TkTBracketR)
                    errMsg(line, "malformed template argument, expected baseclass<T>");            
            }
            
            lb += consumeWS(tokens, index);
            if (tokens[index].type != TkCodeBlock || index+1 != tokens.size())
                errMsg(line, "malformed preprocessor keyword block. Expected 'PYTHON class name : public baseclass { code }'");
            
            return processPythonClass(lb, name, options, templ, baseclass, tokens[index].text, line);
        } 
        else if (type == gParent){
            // constructor
            string name = type;
            type = "";
            
            vector<Argument> args = parseArgs(tokens, index, true, lb);
            
            // if constructor list, collect until codeblock
            string cb = "";
            if (tokens[index].type == TkColon) {
                while(tokens[index].type != TkSemicolon && tokens[index].type != TkCodeBlock) {
                    cb += tokens[index].text;
                    index++;
                    if (index == tokens.size())
                        errMsg(line, "malformed constructor.");
                }
            }
            
            if ( (tokens[index].type != TkCodeBlock && tokens[index].type != TkSemicolon) || index+1 != tokens.size())
                errMsg(line, "malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]");
            return processPythonFunction(lb, name, type, args, cb, tokens[index].text, line);
        } else {
            // member (function)
            if (tokens[index].type != TkDescriptor)
                errMsg(line, "malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]'");
            string name = tokens[index++].text;
            
            lb += consumeWS(tokens, index);
            if (tokens[index].type == TkSemicolon) {
                // member variable
                return processPythonVariable(lb, name, options, type, line);
            }
            vector<Argument> args = parseArgs(tokens, index, true, lb);
            
            if ( (tokens[index].type != TkCodeBlock && tokens[index].type != TkSemicolon) || index+1 != tokens.size())
                errMsg(line, "malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]");
            return processPythonFunction(lb, name, type, args, "", tokens[index].text, line);
        }        
    }
    else 
        errMsg(-1, "preprocessor error, unimplemented Keyword '"+kw+"'. This shouldn't happen.");
    
    return "";
}

// find a function inside a code-snippet and replace the text
string replaceFunctionHeader(const string& text, const string& funcname, const string& header, const string& finalizer, int line) {
    size_t cursor = 0;
    int dline = 0;
    vector<Token> tokens = tokenizeBlock("class", text, cursor, dline, true);
    string newText = "";
    bool found = false;
    
    for (size_t i=0; i < tokens.size(); i++) {
        if (tokens[i].type == TkDescriptor && tokens[i].text == "void") {
            int oldi = i;
            i++;
            consumeWS(tokens, i);
            if (tokens[i].type == TkDescriptor && tokens[i].text == funcname) {            
                i++;
                consumeWS(tokens, i);
                if (tokens[i++].type != TkBracketL) 
                    errMsg(line, "malformed member function. Expected void " + funcname + "(...) {...}");
                while (tokens[i++].type != TkBracketR) {};
                consumeWS(tokens, i);
                if (tokens[i].type != TkCodeBlock)
                    errMsg(line, "malformed member function. Expected void " + funcname + "(...) {...}");
                
                // ok, replace
                newText += header;
                newText += tokens[i].text.substr(1, tokens[i].text.size()-2); // remove brackets
                newText += finalizer;
                found = true;
                continue;
            } else {
                i=oldi;
            }
        }
        newText += tokens[i].text;
    }
    
    if (!found)
        errMsg(line, "preprocessing custom class: can't find member function void " + funcname + "(...)");
    
    return newText;
}
