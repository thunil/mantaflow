/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Preprocessor Declarations
 *
 ******************************************************************************/

#ifndef _PREP_H
#define _PREP_H

#include <string>
#include <vector>

// Tokens
enum Keyword { KwNone = 0, KwKernel, KwPython };

enum TokenType { TkNone = 0, TkComment, TkWhitespace, TkCodeBlock, TkDescriptor, TkComma, TkBracketL, TkBracketR, TkTBracketL, TkTBracketR, TkAssign, TkPointer, TkRef, TkColon, TkSemicolon };

struct Token {
    Token(TokenType t, int l) : type(t), text(""), line(l) {}
    Token(TokenType t, int l, char c) : type(t), text(""), line(l) { text += c; }
    TokenType type;
    std::string text;
    int line;
};

struct Argument {
    Argument() : name(""),value(""),type(""),templ(""), complete(""), isRef(false), isPointer(false), isConst(false), number(0) {}
    std::string name, value, type, templ, complete;
    bool isRef, isPointer, isConst;
    int number;
    std::string getType() const {
        std::string s="";
        if (isConst) s+= "const ";
        s += type;
        if (!templ.empty()) s+= "<" + templ + " >";
        if (isRef) s+= "&";
        if (isPointer) s+= "*";            
        return s;
    }
    std::string getTypeName() const { return getType() + " " + name; }
};
typedef std::vector<Argument> ArgList;

inline std::string listArgs(const ArgList& args) {
    std::string s = "";
    for (size_t i=0; i<args.size(); i++) {
        s += args[i].complete;
        if (i!=args.size()-1) s+=",";
    }
    return s;
}

inline bool isNameChar(char c) {
    return (c>='A' && c<='Z') || (c>='a' && c<='z') || (c>='0' && c<='9') || c=='_';
}

inline bool isIntegral(const std::string& t) {
    return t=="int" || t=="char" || t=="unsigned" || t=="bool" || t=="float" || t=="long" || t=="double" || t=="Real" || t=="Vec3" || t=="Vec3i" || t=="string" || t=="std::string";
}

inline Keyword checkKeyword(const std::string& word) {
    if (word == "KERNEL") return KwKernel;
    if (word == "PYTHON") return KwPython;
    return KwNone;
}

#define assert(x, msg) if(!(x)){errMsg(line,msg);}

// from main.cpp
enum MType { MTNone = 0, MTTBB, MTOpenMP};
extern std::string gFilename;
extern bool gDebugMode;
extern MType gMTType;
extern bool gDocMode;
extern std::string gRegText;
extern std::string gParent;
void errMsg(int line, const std::string& text);
void replaceAll(std::string& text, const std::string& pattern, const std::string& repl);

// functions from merge.cpp
std::string generateMerge(const std::string& text);

// functions from tokenize.cpp
std::string processText(const std::string& text, int baseline);
std::vector<Token> tokenizeBlock(const std::string& kw, const std::string& text, size_t& cursor, int& line, bool complete);

// functions from parse.cpp
std::string parseBlock(const std::string& kw, const std::vector<Token>& tokens, int line);
std::string replaceFunctionHeader(const std::string& text, const std::string& funcname, const std::string& header, const std::string& finalizer, int line);
std::string stripWS(const std::string& s);

// functions from codegen_XXX.cpp
std::string createConverters(const std::string& name, const std::string& tb, const std::string& nl, const std::string& nlr);
std::string processKernel(int lb, const std::string& name, const ArgList& opts, Argument retType, const ArgList& returnArg, const ArgList& templArgs, const ArgList& args, const std::string& code, int line); 
std::string processPythonFunction(int lb, const std::string& name, const ArgList& opts, const std::string& type, const ArgList& args, const std::string& initlist, bool isInline, bool isConst, bool isVirtual, const std::string& code, int line);
std::string processPythonVariable(int lb, const std::string& name, const ArgList& opts, const std::string& type, int line); 
std::string processPythonClass(int lb, const std::string& name, const ArgList& opts, const ArgList& templArgs, const std::string& baseclassName, const ArgList& baseclassTempl, const std::string& code, int line); 
std::string processPythonInstantiation(int lb, const std::string& name, const ArgList& templArgs, const std::string& aliasname, int line); 
std::string buildline(int lb);

#endif // _PREP_H
