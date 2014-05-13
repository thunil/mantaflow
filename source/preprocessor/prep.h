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
#include <sstream>

// Helpers
struct BracketStack {
    bool empty() { return stack.empty(); }
    char top() { return empty() ? '\0' : *(stack.rbegin()); }
    void push_back(char c) { stack += c; }
    char pop() { if (empty()) return '\0'; char c = *(stack.rbegin()); stack.erase(stack.end()-1); return c; }
    
    std::string stack;
};

// Tokens
enum Keyword { KwNone = 0, KwKernel, KwPython };

enum TokenType { TkNone = 0, TkComment, TkWhitespace, TkCodeBlock, TkDescriptor, TkComma, TkBracketL, 
                 TkBracketR, TkTBracketL, TkTBracketR, TkAssign, TkPointer, TkRef, TkDoubleColon, TkSemicolon, 
                 TkSimpleType, TkTypeQualifier, TkConst, TkEnd, TkManta, TkUnsupportedKW, TkClass,
                 TkInline, TkTemplate, TkStatic, TkVirtual, TkString, TkPublic, TkColon };

const std::string unsupportedKeywords[] = {"and", "and", "and_eq", "auto", "bitand", "bitor", "break", 
    "catch", "const_cast", "continue", "default", "delete", "do", "dynamic_cast", "else", "enum", 
    "explicit", "export", "extern", "for", "friend", "goto", "if", "mutable", "namespace", "new", 
    "not", "not_eq", "or", "or_eq", "operator", "private", "protected", "register", 
    "reinterpret_cast", "return", "sizeof", "static_cast", "switch", "this", "throw", "try", "typedef", 
    "union", "using", "volatile", "while", "xor", "xor_eq", "" };

struct Token {
    Token(TokenType t, int l) : type(t), text(""), line(l) {}
    Token(TokenType t, int l, char c) : type(t), text(""), line(l) { text += c; }
    TokenType type;
    std::string text;
    int line;
};

// tracks a set of tokens, and the current position in this list
struct TokenPointer {
    TokenPointer(const std::vector<Token>& t) : parent(0),queue(t),ptr(0),oldPtr(0),linebreaks(0) {}
    TokenPointer(TokenPointer& t) : parent(&t),queue(t.queue),ptr(t.ptr),oldPtr(t.ptr),linebreaks(0) {}
    ~TokenPointer();
    TokenPointer *parent;
    const std::vector<Token>& queue;
    int ptr, oldPtr;
    int linebreaks;
    std::string completeText, minimalText;

    inline void reset() { linebreaks = 0; completeText=minimalText=""; consumeWhitespace(); }
    
    inline TokenType curType() { return done() ? TkEnd : cur().type; }
    inline int curLine() { if (queue.empty()) return -1; return done() ? queue.back().line : cur().line; }
    inline const Token& cur() { return queue[ptr]; }
    inline bool done() { return ptr >= queue.size(); }
    inline bool next() { oldPtr = ptr; ptr++; consumeWhitespace(); return !done(); }
    inline void previous() { ptr = oldPtr; consumeWhitespace(); }
    inline bool isLast() { return ptr = queue.size()-1;}
    void consumeWhitespace();
};

struct Type {
    Type() : isConst(false), isRef(false), isPointer(false) {};
    
    std::string name;
    bool isConst, isRef, isPointer;
    std::vector<Type> templateTypes;
};

struct Argument {
    Argument() : type() {};

    Type type;
    std::string name;
};
typedef std::vector<Argument> ArgList;

struct Function {
    Function() : returnType(),isInline(false),isVirtual(false),isConst(false),noParentheses(false) {}

    std::string name;
    Type returnType;
    bool isInline, isVirtual, isConst, noParentheses;
    std::vector<Type> templateTypes;
    std::vector<Argument> arguments;
};

struct Class {
    Class() {};

    std::string name;
    Type baseClass;
    std::vector<Type> templateTypes;  
};

/*
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
*/
/*
inline std::string listArgs(const ArgList& args) {
    std::string s = "";
    for (size_t i=0; i<args.size(); i++) {
        s += args[i].complete;
        if (i!=args.size()-1) s+=",";
    }
    return s;
}
*/
inline bool isNameChar(char c) {
    return (c>='A' && c<='Z') || (c>='a' && c<='z') || (c>='0' && c<='9') || c=='_';
}

inline bool isIntegral(const std::string& t) {
    return t=="int" || t=="char" || t=="unsigned" || t=="bool" || t=="float" || t=="long" || t=="double" || t=="Real" || t=="Vec3" || t=="Vec3i" || t=="string" || t=="std::string";
}

#define debMsg(l, msg) if(gDebugMode){ std::ostringstream out; out << msg; debMsgHelper(l, out.str()); }

// from main.cpp
enum MType { MTNone = 0, MTTBB, MTOpenMP};
extern std::string gFilename;
extern bool gDebugMode;
extern MType gMTType;
extern bool gDocMode;
extern std::string gRegText;
extern std::string gParent;
void errMsg(int line, const std::string& text);
void debMsgHelper(int line, const std::string& text);
void replaceAll(std::string& text, const std::string& pattern, const std::string& repl);

// functions from merge.cpp
std::string generateMerge(const std::string& text);

// functions from tokenize.cpp
std::string processText(const std::string& text, int baseline);

// functions from parse.cpp
std::string parseBlock(const std::string& kw, std::vector<Token>& tokens);
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
