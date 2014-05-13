/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
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

struct Text {
    int line0;
    std::string minimal, original;
    void reset() { minimal = original = ""; line0=0; }
    void add(Text* a) { minimal += a->minimal; original += a->original; }
    std::string linebreaks() const; 
};

// tracks a set of tokens, and the current position in this list
struct TokenPointer {
    TokenPointer(const std::vector<Token>& t, Text *track) : parent(0),queue(t),ptr(0),txt(track) { reset(); }
    TokenPointer(TokenPointer& t, Text *track) : parent(&t),queue(t.queue),ptr(t.ptr),txt(track) { reset(); }
    ~TokenPointer();
    TokenPointer *parent;
    const std::vector<Token>& queue;
    int ptr;
    Text *txt;

    inline void reset() { txt->reset(); if(!done()) txt->line0 = cur().line; consumeWhitespace(); }
    inline TokenType curType() { return done() ? TkEnd : cur().type; }
    inline const Token& cur() { return queue[ptr]; }
    inline TokenType previewType() { return (ptr+1 >= queue.size()) ? TkEnd : queue[ptr+1].type; }
    inline bool done() { return ptr >= queue.size(); }
    inline bool isLast() { return ptr = queue.size()-1;}
    void next();
    void consumeWhitespace();
    void errorMsg(const std::string& msg);
};

template <class T>
struct List : Text {
    std::vector<T> _data;
    inline size_t size() const { return _data.size(); }
    inline T& operator[](int i) { return _data[i]; }
    inline const T& operator[](int i) const { return _data[i]; }
    inline bool empty() const { return _data.empty(); }
    inline T& back() { return _data.back(); }
    inline void push_back(const T& a) { _data.push_back(a); }
};

struct Type : Text {
    Type() : isConst(false), isRef(false), isPointer(false) {};
    
    std::string name;
    bool isConst, isRef, isPointer;
    List<Type> templateTypes;

    bool operator==(const Type& a) const;
    std::string build() const;
};

struct Argument : Text {
    Argument() : type(),index(-1) {};
    
    Type type;
    int index;
    std::string name, value;
    std::string completeText, minimalText;
};

struct Function : Text {
    Function() : returnType(),isInline(false),isVirtual(false),isConst(false),noParentheses(false) {}

    std::string name;
    Type returnType;
    bool isInline, isVirtual, isConst, noParentheses;
    List<Type> templateTypes;
    List<Argument> arguments;
};

struct Class : Text {
    Class() {};

    std::string name;
    Type baseClass;
    List<Type> templateTypes;  
};

struct Block : Text {
    Block() {};

    std::string initList;
    Type aliasType;
    std::string aliasName;
    Class cls;
    Function func;
    List<Argument> options;
    List<Argument> reduceArgs;
};

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
std::string parseBlock(const std::string& kw, const std::vector<Token>& tokens);

// functions from codegen_XXX.cpp
std::string createConverters(const std::string& name, const std::string& tb, const std::string& nl, const std::string& nlr);
std::string processKernel(const Block& block, const std::string& code);
std::string processPythonFunction(const Block& block, const std::string& code);
std::string processPythonVariable(const Block& block, const std::string& code);
std::string processPythonClass(const Block& block, const std::string& code);
std::string processPythonInstantiation(const Block& block, const std::string& code);

#endif // _PREP_H
