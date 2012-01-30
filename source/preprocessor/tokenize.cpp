/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor Tokenizer
 *
 ******************************************************************************/

#include <string>
#include <iostream>
#include <vector>

#include "prep.h"

using namespace std;

//*************************************************************************************
// Helpers

struct BracketStack {
    bool empty() { return stack.empty(); }
    char top() { return empty() ? '\0' : *(stack.rbegin()); }
    void push(char c) { stack += c; }
    char pop() { if (empty()) return '\0'; char c = *(stack.rbegin()); stack.erase(stack.end()-1); return c; }
    
    string stack;
};

//*************************************************************************************
// Lexing functions

// parse complete file.
// Defer parts to be processed to processTextBlock, and directly copy the rest
string processText(const string& text, int baseline) {
    string newText = "";
    
    // no real lexing yet, only track define and comment blocks
    string word = "";
    int line=baseline;
    bool comment=false, slComment=false, define=false, extendLine=false;
    for (size_t i=0; i<text.size()-1; i++) {
        char c = text[i];
        
        // track lines
        bool isEOL = !extendLine && (c=='\n' || c=='\r');
        if (c=='\\') 
            extendLine = true;
        else if (c!='\n' && c!='\r') 
            extendLine = false;
        if (c=='\n') line++;
        
        // track comments and defines
        if (comment) {
            if (c=='*' && text[i+1]=='/') comment = false;            
            newText += c;
        } 
        else if (slComment) {
            if (isEOL) slComment = false;
            newText += c;
        } 
        else if (define) {
            if (isEOL) define = false;
            newText += c;
        }
        else if (c=='/' && text[i+1]=='*') {
            comment = true;
            newText += word;
            word = "";
            newText += c;
        }
        else if (c=='/' && text[i+1]=='/') {
            slComment = true;
            newText += word;
            word = "";
            newText += c;
        }
        else if (c=='#') {
            define = true;
            newText += word;
            word = "";
            newText += c;            
        } 
        else {
            // tokenize keywords
            if (isNameChar(c))
                word += c;
            else {
                if (word.empty() || checkKeyword(word) == KwNone) {
                    newText += word;                    
                    newText += c;                
                } else {
                    vector<Token> tokens = tokenizeBlock(word, text, i, line, false);
                    newText += parseBlock(word, tokens, line);
                }
                word = "";         
            }
        }
    }
    newText += word;
    
    return newText;    
}

// tokenize and parse until keyword section ends
vector<Token> tokenizeBlock(const string& kw, const string& text, size_t& i, int& line, bool complete) {
    vector<Token> tokens;
    tokens.push_back(Token(TkWhitespace, line));
    BracketStack brackets;
        
    // tokenize loop
    bool comment=false, slComment=false, define=false, extendLine=false;
    int codeblockLevel = 0;
    for (; i<text.size(); i++) {
        char c = text[i];
        string& curstr = tokens.rbegin()->text;
        bool lastchar = i == text.size()-1 ;
        
        // track lines
        bool isEOL = !extendLine && (c=='\n' || c=='\r');
        if (c=='\\') 
            extendLine = true;
        else if (c!='\n' && c!='\r') 
            extendLine = false;
        if (c=='\n') line++;
        
        // track comments and defines
        if (comment) {
            if (!lastchar && c=='*' && text[i+1]=='/') comment = false;            
            curstr += c;
        }
        else if (slComment) {
            if (isEOL) slComment = false;
            curstr += c;
        } 
        else if (define) {
            if (isEOL) define = false;
            curstr += c;
        }
        else if (!lastchar && c=='/' && text[i+1]=='*') {
            comment = true;
            if (codeblockLevel==0) 
                tokens.push_back(Token(TkComment, line, c));
            else
                curstr += c;
        }
        else if (!lastchar && c=='/' && text[i+1]=='/') {
            slComment = true;
            if (codeblockLevel==0) 
                tokens.push_back(Token(TkComment, line, c));
            else
                curstr += c;
        }
        else if (c=='#') {
            define = true;
            if (codeblockLevel==0) 
                tokens.push_back(Token(TkComment, line, c));
            else
                curstr += c;
        }   
        
        // track codeblock
        else if (codeblockLevel > 0) {
            curstr += c;
            if (c == '}') codeblockLevel--;
            if (c == '{') codeblockLevel++;
            
            if (codeblockLevel == 0 && !complete) {
                // block finished, return to parseText
                return tokens;
            }
        }
        else if (c=='{') {
            codeblockLevel++;
            if (!brackets.empty())
                errMsg(line, "codeblock {} is not allowed inside brackets.");                
            tokens.push_back(Token(TkCodeBlock, line, c));
        }
                
        // track brackets
        else if (c=='(') {
            tokens.push_back(Token(TkBracketL, line, c));
            brackets.push(c);
        }
        else if (c==')') {
            if (brackets.pop() != '(') 
                errMsg(line, "bracket mismatch, can't close ')'");
            tokens.push_back(Token(TkBracketR, line, c));
        }
        else if (c=='<') {
            tokens.push_back(Token(TkTBracketL, line, c));
            brackets.push(c);
        }
        else if (c=='>') {
            if (brackets.pop() != '<')
                errMsg(line, "bracket mismatch, can't close '>'");
            tokens.push_back(Token(TkTBracketR, line, c));
        }
        
        // track symbol tokens
        else if (c==' ' || c=='\t' || c=='\r' || c=='\n') {
            tokens.push_back(Token(TkWhitespace, line, c));
        }
        else if (c=='=') {
            tokens.push_back(Token(TkAssign, line, c));
        }
        else if (c==',') {
            tokens.push_back(Token(TkComma, line, c));
        }
        else if (c=='*') {
            tokens.push_back(Token(TkPointer, line, c));
        }
        else if (c=='&') {
            tokens.push_back(Token(TkRef, line, c));
        }
        else if (c==':') {
            tokens.push_back(Token(TkColon, line, c));
        }
        else if (c==';') {
            if (!brackets.empty())
                errMsg(line, "Semicolon is not allowed inside brackets.");                
            
            tokens.push_back(Token(TkSemicolon, line, c));
                        
            // block finished, return to parseText
            if (!complete) return tokens;
        }
        
        // track descriptors
        else if (isNameChar(c)) {
            if (tokens.rbegin()->type != TkDescriptor)
                tokens.push_back(Token(TkDescriptor, line, c));
            else
                curstr += c;
        }
        // misc stuff: don't bother with more complex tokens
        else {
            if (tokens.rbegin()->type != TkNone) {
                tokens.push_back(Token(TkNone, line, c));
            } else
                curstr += c;
        }        
    }
    
    // EOF before block end ?
    if (!complete)
        errMsg(line, "EOF before block for keyword '" + kw + "' was closed.");
    
    return tokens;
}
