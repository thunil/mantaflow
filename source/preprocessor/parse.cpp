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

string stripWS(const string& s) {
    string t;
    for (size_t i=0; i<s.length(); i++) {
        if (s[i] != ' ' && s[i] != '\r' && s[i] != '\t' && s[i] != '\n') t += s[i];
    }
    return t;
}

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

// parse a single argument
Argument parseSingleArg(const vector<Token>& tokens, size_t& index, bool expectType, bool expectName, bool isList, int& lb) {
    const int firstStage = expectType ? 0 : 4;
    
    lb += consumeWS(tokens, index);
    
    // tokenizer already established that open brackets are closed eventually, so no need to check that
    Argument cur = Argument();    
    int stage = firstStage; // 0: type, 1: template open / modifiers / name 2: template, 3: template close 4: name / modifiers, 5: assignOp 6+: defaults
    int openTBrackets=0, openBrackets=0;
    bool endPossible = false, endToken = false;
        
    do {
        TokenType t = tokens[index].type;
        string text = tokens[index].text;
        int line = tokens[index].line;
        
        if (t==TkBracketL) openBrackets++;
        if (t==TkBracketR) openBrackets--;
        if (t==TkTBracketL) openTBrackets++;
        if (t==TkTBracketR) openTBrackets--;
        
        //cout << index << ":" << stage << "," << t << "," << text << " : " << cur.complete << endl;
        endPossible = openBrackets == 0 && openTBrackets == 0;
        endToken = false;
        if (t == TkWhitespace || t == TkComment) {
            if (text[0]=='\n') lb++;
            else if (text[0]!='\r') cur.complete += text;
        }
        else
        {
            cur.complete += text;
            switch (stage) {
                case 0:
                    // type name
                    assert(t==TkDescriptor, "incomplete argument! expect argument of the form '[type] name [= default]'");
                    if (text == "const") {
                        cur.isConst = true;     
                        endPossible = false;
                    } 
                    else {
                        cur.type += text;
                        stage++;
                        if (expectName) endPossible = false;
                    }
                    break;
                case 1:
                    if (t == TkColon && tokens[index+1].type == TkColon) { // type had namespace
                        cur.type += "::";
                        cur.complete += ":";
                        index++;
                        stage=0;
                        endPossible = false;
                        break;
                    } // else : fallthough to name
                case 4:
                    // template open / name / modifiers
                    if (t == TkTBracketL && stage == 1) {
                        stage = 2;
                        endPossible = false;
                    } else if (t == TkRef) {
                        cur.isRef = true;
                        stage = 4;
                        if (expectName) endPossible = false;
                    } else if (t == TkPointer) {
                        cur.isPointer = true;
                        stage = 4;
                        if (expectName) endPossible = false;
                    } else if (t == TkDescriptor) {
                        assert(expectName, "Argument '"+ cur.complete +"': only type expected");
                        cur.name = text;
                        stage = 5;
                    } else 
                        errMsg(line, "incomplete argument '" + cur.complete + "'");
                    break;
                case 2:
                    // template
                    endPossible = false;
                    assert(t == TkDescriptor, "incomplete argument type '" + cur.complete + "'! expect type<Template> [*|&]");
                    cur.templ += text;
                    stage++;
                    break;
                case 3:
                    // template close
                    endPossible = false;
                    if (t == TkTBracketR)
                        stage = 4;
                    else if (t == TkComma) {
                        cur.templ += ",";
                        stage = 2;
                    }
                    else 
                        errMsg(line, "incomplete argument type! expect type<Template> [*|&]");
                    break;
                case 5:
                    // assign op
                    endPossible = false;
                    assert(t == TkAssign, "incomplete argument! expect argument of the form '[type] name [= default]'");
                    stage++;
                    break;
                case 6:
                    // default values, stage 6+
                    if (t == TkNone || t == TkDescriptor || t==TkColon || t==TkComma || t==TkBracketL || t==TkBracketR)
                        cur.value += text;
                    else 
                        errMsg(line, "incomplete argument '" + cur.complete + "'");
                    break;
                default:
                    errMsg(line, "internal parser error, invalid stage");
                    break;
            }
        }
        index++;
        if (index >= tokens.size())
            errMsg(-1, "Preprocessor ran out of tokens. This shouldn't happen.");
        endToken = (!isList && tokens[index].type == TkDescriptor) ||
                   (isList && (tokens[index].type == TkBracketR || tokens[index].type == TkComma || tokens[index].type == TkTBracketR));
    } while(!endPossible || !endToken);
    
    //cout << "Complete Type : " << cur.complete << endl;
    
    return cur;
}


// parse argument list tokens of the form (int a = 12, bool* b, Grid<float>& grid)
// or <class T, int A> if isTemplate is specified
ArgList parseArgs(const vector<Token>& tokens, size_t& index, bool expectType, int &lb, bool isTemplate) {    
    ArgList args;
    TokenType openBracket = isTemplate ? TkTBracketL : TkBracketL;
    TokenType closeBracket = isTemplate ? TkTBracketR : TkBracketR;
    
    // ignore trailing whitespaces, remove first bracket
    lb += consumeWS(tokens, index);
    if (tokens[index].type != openBracket) return args;
    index++;
    lb += consumeWS(tokens, index);
    
    for(;;) {
        if (tokens[index].type == closeBracket) 
            break;
        
        Argument cur = parseSingleArg(tokens, index, expectType, true, true, lb);
        args.push_back(cur);
        
        if (tokens[index].type == closeBracket) 
            break;
        else if (tokens[index].type != TkComma)
            errMsg(tokens[index].line, "invalid argument list !");
        index++;
        lb += consumeWS(tokens, index);
    }
    
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
    
    ArgList options = parseArgs(tokens, index, false, lb, false);
    
    assert(tokens[index].type == TkDescriptor, "malformed preprocessor keyword block. Expected '" + kw + "(opt1, opt2, ...) STATEMENTS [{}|;]'");
    
    // keyword dependent processing. 
    // No need to check for index overflow, as long as test each type, since sequence will end in ; or {}
    if (key == KwKernel) {
        string text = tokens[index].text;
        ArgList templArgs;
        
        // resolve template class or instantiation
        if (text == "template") {
            // templateted kernel
            lb += consumeWS(tokens, ++index);
            assert (tokens[index].type == TkTBracketL, "syntax error. expected KERNEL(...) template<class T> [class] ...");
            templArgs = parseArgs(tokens, index, true, lb, true);
            assert(tokens[index].type == TkDescriptor, "syntax error. expected KERNEL(...) template<class T> [class] ...");
            text = tokens[index].text;            
        }                
        
        // parse return type 
        Argument retType = parseSingleArg(tokens, index, true, false, false, lb);
        
        assert (tokens[index].type == TkDescriptor, "Malformed KERNEL, expected: KERNEL ret_type/struct name(args...) { code }");
        string name = tokens[index++].text;
        ArgList args = parseArgs(tokens, index, true, lb, false);
        
        if (tokens[index].type != TkCodeBlock || index+1 != tokens.size())
            errMsg(line, "Malformed KERNEL, expected: KERNEL(opts...) ret_type name(args...) { code }");
        
        bool isClass = retType.type == "struct" || retType.type == "class";
        if (isClass) {
            retType.complete = "void";
            retType.type = "void";
        }
        return processKernel(lb, name, options, retType, templArgs, args, tokens[index].text, line, isClass);
    }
    else if (key == KwPython) 
    {
        ArgList templArgs;
        string type = tokens[index].text; 
        
        // resolve template alias
        if (type == "alias") {            
            // template instantiation
            
            lb += consumeWS(tokens, ++index);
            assert (tokens[index].type == TkDescriptor, "syntax error. expected PYTHON alias type<T> name;");
            string classname = tokens[index++].text;
            
            templArgs = parseArgs(tokens, index, false, lb, true);
            
            assert (tokens[index].type == TkDescriptor, "syntax error. expected PYTHON alias type<T> name;");
            string aliasname = tokens[index++].text;
            lb += consumeWS(tokens, index);
            
            if (tokens[index].type != TkSemicolon || index+1 != tokens.size())
                errMsg(line, "syntax error. expected PYTHON alias type<T> name;");            
            return processPythonInstantiation(lb, classname, templArgs, aliasname, line);
        }
        // resolve template class
        else if (type == "template") {
            
            lb += consumeWS(tokens, ++index);        
            if (tokens[index].type == TkTBracketL) {
                // class template                
                templArgs = parseArgs(tokens, index, true, lb, true);
                
                if (tokens[index].type != TkDescriptor)
                    errMsg(line, "syntax error. expected PYTHON template<class T> class ...");
                type = tokens[index].text;                
            }
            else 
                errMsg(line, "syntax error. expected PYTHON class name<T>; or PYTHON template<class T> class name...");                
        }
                
        if (type == "class") {
            lb += consumeWS(tokens, ++index);                    
            assert(tokens[index].type == TkDescriptor, "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");
            string name = tokens[index++].text;
            
            lb += consumeWS(tokens, index);
            assert(tokens[index++].type == TkColon, "PYTHON class '" + name + "' must publicly derive from PbClass (or a subclass)");
            
            lb += consumeWS(tokens, index);
            assert(tokens[index].type == TkDescriptor && tokens[index++].text == "public", "PYTHON class '" + name + "' must publicly derive from PbClass (or a subclass)");
            
            lb += consumeWS(tokens, index);
            assert(tokens[index].type == TkDescriptor, "PYTHON class '" + name + "' must publicly derive from PbClass (or a subclass)");
            string baseclassName = tokens[index++].text;
            
            lb += consumeWS(tokens, index);
            ArgList baseclassTempl;
            if (tokens[index].type == TkTBracketL) {
                baseclassTempl = parseArgs(tokens, index, false, lb, true);
            }
            
            lb += consumeWS(tokens, index);
            if (tokens[index].type != TkCodeBlock || index+1 != tokens.size())
                errMsg(line, "malformed preprocessor keyword block. Expected 'PYTHON class name : public baseclass { code }'");
            
            return processPythonClass(lb, name, options, templArgs, baseclassName, baseclassTempl, tokens[index].text, line);
        } 
        else if (type == gParent){
            lb += consumeWS(tokens, ++index);                    
            // constructor
            string name = type;
            type = "";
            
            ArgList args = parseArgs(tokens, index, true, lb, false);
            
            // if constructor list, collect until codeblock
            string cb = "";
            if (tokens[index].type == TkColon) {
                while(tokens[index].type != TkSemicolon && tokens[index].type != TkCodeBlock) {
                    cb += tokens[index].text;
                    index++;
                    assert(index != tokens.size(), "malformed constructor.");
                }
            }
            
            if ( (tokens[index].type != TkCodeBlock && tokens[index].type != TkSemicolon) || index+1 != tokens.size())
                errMsg(line, "malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]");
            return processPythonFunction(lb, name, type, args, cb, tokens[index].text, line);
        } else {
            // parse return type 
            Argument retType = parseSingleArg(tokens, index, true, false, false, lb);
            type = stripWS(retType.complete);

            // function or member function
            assert(tokens[index].type == TkDescriptor, "malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]'");
            string name = tokens[index++].text;
            
            lb += consumeWS(tokens, index);
            if (tokens[index].type == TkSemicolon) {
                // member variable
                return processPythonVariable(lb, name, options, type, line);
            }
            ArgList args = parseArgs(tokens, index, true, lb, false);
            
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
                assert(tokens[i++].type == TkBracketL, "malformed member function. Expected void " + funcname + "(...) {...}");
                while (tokens[i++].type != TkBracketR) {};
                consumeWS(tokens, i);
                assert(tokens[i].type == TkCodeBlock, "malformed member function. Expected void " + funcname + "(...) {...}");
                
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
    
    assert(found, "preprocessing custom class: can't find member function void " + funcname + "(...)");
    
    return newText;
}
