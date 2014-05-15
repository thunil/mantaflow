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
#include <fstream>

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
    std::string linebreaks() const; 
    virtual std::string dynamicClass() { return ""; }
    void prequel(const Text* a) { minimal = a->minimal + minimal; original = a->original + original; }
};

// tracks a set of tokens, and the current position in this list
struct TokenPointer {
    TokenPointer(const std::vector<Token>& t, Text *track) : parent(0),queue(t),ptr(0),txt(track) { reset(); }
    TokenPointer(TokenPointer& t, Text *track) : parent(&t),queue(t.queue),ptr(t.ptr),txt(track) { reset(); }
    TokenPointer *parent;
    const std::vector<Token>& queue;
    int ptr;
    Text *txt;

    inline void reset() { txt->reset(); consumeWhitespace(); if(!done()) txt->line0 = cur().line;  }
    inline TokenType curType() { return done() ? TkEnd : cur().type; }
    inline const Token& cur() { return queue[ptr]; }
    inline bool done() { return ptr >= queue.size(); }
    inline bool isLast() { return ptr = queue.size()-1;}
    void forward(const std::string& minimal, const std::string& original, int offset);
    void next();
    void consumeWhitespace();
    void errorMsg(const std::string& msg);
    std::string backtrace();
};

template <class T>
struct List : Text {
    std::vector<T> _data;
    std::string listText;
    inline size_t size() const { return _data.size(); }
    inline T& operator[](int i) { return _data[i]; }
    inline const T& operator[](int i) const { return _data[i]; }
    inline bool empty() const { return _data.empty(); }
    inline T& back() { return _data.back(); }
    inline void push_back(const T& a) { _data.push_back(a); }
    virtual std::string dynamicClass() { return "List"; }
};

struct Type : Text {
    Type() : isConst(false), isRef(false), isPointer(false) {};
    
    std::string name;
    bool isConst, isRef, isPointer;
    List<Type> templateTypes;

    inline bool isTemplated() const { return !templateTypes.empty(); }
    std::string tplString() const;
    bool operator==(const Type& a) const;
    std::string build() const;
    virtual std::string dynamicClass() { return "Type"; }
};

struct Argument : Text {
    Argument() : type(),index(-1) {};
    
    Type type;
    int index;
    std::string name, value;
    std::string completeText, minimalText;
    virtual std::string dynamicClass() { return "Argument"; }
};

struct Function : Text {
    Function() : returnType(),isInline(false),isVirtual(false),isConst(false),noParentheses(false) {}

    std::string name;
    Type returnType;
    bool isInline, isVirtual, isConst, noParentheses;
    List<Type> templateTypes;
    List<Argument> arguments;
    std::string callString() const;
    virtual std::string dynamicClass() { return "Function"; }
};

struct Class : Text {
    Class() {};

    std::string name;
    Type baseClass;
    List<Type> templateTypes;  
    std::string tplString() const;
    std::string fullName() const { return isTemplated() ? (name+"<"+tplString()+">") : name; }
    inline bool isTemplated() const { return !templateTypes.empty(); }
    virtual std::string dynamicClass() { return "Class"; }
};

struct Block : Text {
    Block() {};

    std::string initList;
    Class cls;
    const Class *parent;
    Function func;
    List<Argument> options;
    List<Argument> reduceArgs;
    virtual std::string dynamicClass() { return "Block"; }
};

// defined in util.cpp
struct Sink {
    Sink(const std::string& file);
    void write();

    std::ostringstream inplace;
    std::ostringstream link;
    bool isHeader;
private:
    std::string filename;
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

// functions from util.cpp
void generateMerge(int num, char* files[]);
void errMsg(int line, const std::string& text);
void debMsgHelper(int line, const std::string& text);
void replaceAll(std::string& text, const std::string& pattern, const std::string& repl);
std::string readFile(const std::string&);
void writeFile(const std::string& name, const std::string& text);
std::string replaceSet(const std::string& templ, const std::string table[]);
std::vector<std::string> split(const std::string& text, char sep);

// functions from merge.cpp
std::string generateMerge(const std::string& text);

// functions from tokenize.cpp
void processText(const std::string& text, int baseline, Sink& sink, const Class *parent);

// functions from parse.cpp
void parseBlock(const std::string& kw, const std::vector<Token>& tokens, const Class *parent, Sink& sink);

// functions from codegen_XXX.cpp
std::string createConverters(const std::string& name, const std::string& tb, const std::string& nl, const std::string& nlr);
std::string processKernel(const Block& block, const std::string& code);
void processPythonFunction(const Block& block, const std::string& code, Sink& sink);
void processPythonVariable(const Block& block, Sink& sink);
void processPythonClass(const Block& block, const std::string& code, Sink& sink);
void processPythonInstantiation(const Block& block, const Type& aliasType, const std::string& aliasName, Sink& sink);

#endif // _PREP_H
