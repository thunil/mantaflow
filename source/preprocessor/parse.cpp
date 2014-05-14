/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2014 Tobias Pfaff, Nils Thuerey 
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
#include <algorithm>
#include "prep.h"

using namespace std;

List<Type> parseTypeList(TokenPointer& parentPtr);

//*************************************************************
// Helpers

bool Type::operator==(const Type& a) const {
    if (a.templateTypes.size() != templateTypes.size()) return false;
    for (size_t i=0; i<templateTypes.size(); i++)
        if (templateTypes[i].name != a.templateTypes[i].name) return false;
    return a.name==name && a.isConst==isConst && a.isRef==isRef && a.isPointer==isPointer;
}

string Type::build() const {
    string s="";
    if (isConst) s+= "const ";
    s += name;
    if (!templateTypes.empty())
        s+= templateTypes.minimal;
    if (isRef) s+= "&";
    if (isPointer) s+= "*";            
    return s;
}

string Text::linebreaks() const {
    int num = count(original.begin(), original.end(), '\n') - 
        count(minimal.begin(),minimal.end(), '\n');
    string s="";
    for (int i=0; i< num; i++)
        s += '\n';
    return s;
}


//*************************************************************
// parsers

#define tkAssert(x,msg) {if(!(x)){tk.errorMsg(msg);}}
#define typeAssert(x) tkAssert(x," is not a valid type.")
#define argAssert(x) tkAssert(x," is not a valid argument.")

string parseRunaway(TokenPointer& parentPtr) {
    Text text;
    TokenPointer tk(parentPtr, &text);

    // don't validate, just track bracket level
    BracketStack stack;
    for(;!tk.done();tk.next()) {
        if (stack.empty() && (tk.curType() == TkComma || 
            tk.curType() == TkBracketR || tk.curType() == TkTBracketR)) {
                break;
        }
        if (tk.curType() == TkBracketL || tk.curType() == TkTBracketL) {
            stack.push_back(tk.cur().text[0]);
        } else if (tk.curType() == TkBracketR) {
            argAssert(stack.pop() == '(');
        } else if (tk.curType() == TkTBracketR) {
            argAssert(stack.pop() == '<');
        }
    }
    argAssert(stack.empty());
    return text.minimal;
}

void parsePointer(TokenPointer& tk, Type& cur) {
    if (tk.done()) return;
    if (tk.cur().type == TkPointer) {
        cur.isPointer = true;
        tk.next();
    } else if (tk.cur().type == TkRef) {
        cur.isRef = true;
        tk.next();
    }
}

Type parseType(TokenPointer& parentPtr) {
    Type cur = Type();
    TokenPointer tk(parentPtr, &cur);

    // constness
    if (tk.curType() == TkConst) {
        cur.isConst = true;
        tk.next();
    }
    typeAssert(!tk.done());
    
    // signed / unsigned
    if (tk.curType() == TkTypeQualifier) {
        cur.name = tk.cur().text;
        tk.next();
        if (tk.done())
            return cur;
        typeAssert(tk.curType() == TkSimpleType);
        cur.name += " " + tk.cur().text;
        tk.next();
        parsePointer(tk, cur);
        return cur;
    }

    // template argument
    if (tk.curType() == TkClass) {
        tk.next();
    }

    typeAssert(tk.curType() == TkDescriptor || tk.curType() == TkSimpleType);
    cur.name = tk.cur().text;
    tk.next();

    // namespace
    if (tk.curType() == TkDoubleColon) {
        cur.name += "::";
        tk.next();
        typeAssert(tk.curType() == TkDescriptor || tk.curType() == TkSimpleType);
        cur.name += tk.cur().text;
        tk.next();
    }

    // template
    if (tk.curType() == TkTBracketL) {
        cur.templateTypes = parseTypeList(tk);
    }

    parsePointer(tk, cur);
    return cur;
}

Argument parseArgument(TokenPointer& parentPtr, bool requireName, bool requireType) {
    Argument cur = Argument();
    TokenPointer tk(parentPtr, &cur);

    if (requireType)
        cur.type = parseType(tk);

    if (tk.curType() != TkDescriptor && !requireName) 
        return cur;

    argAssert(tk.curType() == TkDescriptor);
    cur.name = tk.cur().text;
    tk.next();

    // default value ?
    if (tk.curType() == TkAssign || tk.curType() == TkBracketL) {
        if (tk.curType() == TkAssign)
            tk.next();
        cur.value = parseRunaway(tk);
    }
    return cur;
}

List<Type> parseTypeList(TokenPointer& parentPtr) {
    List<Type> list;
    TokenPointer tk(parentPtr, &list);
    
    tkAssert(tk.curType() == TkTBracketL, "expect template opening bracket");
    tk.next();
    if (tk.curType() != TkTBracketR) {
        for(;;) {
            list.push_back(parseType(tk));
            if (tk.curType() == TkTBracketR) 
                break;
            tkAssert(tk.curType() == TkComma, "expect comma or closing bracket");
            tk.next();
        }
    }
    list.listText = list.minimal.substr(1);
    tkAssert(tk.curType() == TkTBracketR, "expect template closing bracket");
    tk.next();
    return list;
}

List<Argument> parseArgumentList(TokenPointer& parentPtr, bool requireName, bool requireType) {
    List<Argument> list;
    TokenPointer tk(parentPtr, &list);
    
    tkAssert(tk.curType() == TkBracketL, "expect opening bracket");
    tk.next();
    if (tk.curType() != TkBracketR) {
        for(int idx=0;;idx++) {
            list.push_back(parseArgument(tk, requireName, requireType));
            list.back().index = idx;
            if (tk.curType() == TkBracketR) 
                break;
            tkAssert(tk.curType() == TkComma, "expect comma or closing bracket");
            tk.next();
        }
    }
    list.listText = list.minimal.substr(1);
    tkAssert(tk.curType() == TkBracketR, "expect closing bracket");
    tk.next();
    return list;
}

Function parseFunction(TokenPointer& parentPtr, bool requireNames, bool requireType, bool requireArgs) {
    Function cur;
    TokenPointer tk(parentPtr, &cur);

    // templated
    if (tk.curType() == TkTemplate) {
        tk.next();
        cur.templateTypes = parseTypeList(tk);            
    }

    for (;tk.curType() == TkInline || tk.curType() == TkVirtual; tk.next()) {
        if (tk.curType() == TkInline) cur.isInline = true;
        if (tk.curType() == TkVirtual) cur.isVirtual = true;
    }

    if (requireType)
        cur.returnType = parseType(tk);
    tkAssert(tk.curType() == TkDescriptor, "malformed function/kernel");
    cur.name = tk.cur().text;
    tk.next();
    if (requireArgs || tk.curType() == TkBracketL)
        cur.arguments = parseArgumentList(tk, requireNames, true);
    else
        cur.noParentheses = true;

    if (tk.curType() == TkConst) {
        cur.isConst = true;
        tk.next();
    }

    return cur;
}

// Parse syntax KEYWORD(opt1, opt2, ...) STATEMENTS [ {} or ; ]    
string parseBlock(const string& kw, const vector<Token>& tokens) {
    Block block = Block();
    TokenPointer tk(tokens, &block);

    // parse keyword options
    if (tk.curType() == TkBracketL)
        block.options = parseArgumentList(tk, true, false);

    if (kw == "KERNEL") {
        List<Type> templTypes;

        // templated kernel
        if (tk.curType() == TkTemplate) {
            tk.next();
            templTypes = parseTypeList(tk);            
        }
        
        // return values
        while (tk.curType() == TkDescriptor && tk.cur().text == "returns") {
            tk.next();
            tkAssert(tk.curType() == TkBracketL, "expext opening bracket");
            tk.next();
            block.reduceArgs.push_back(parseArgument(tk, true, true));
            tkAssert(tk.curType() == TkBracketR, "expect closing bracket");
            tk.next();            
        }

        block.func = parseFunction(tk, true, true, true);
        if (!templTypes.empty())
            block.func.templateTypes = templTypes;

        tkAssert(tk.curType() == TkCodeBlock && tk.isLast(), 
            "Malformed KERNEL, expected KERNEL(opts...) ret_type name(args...) { code }");

        return processKernel(block, tk.cur().text);
    }
    else if (kw == "PYTHON")
    {
        // template instantiation / alias
        if (tk.curType() == TkDescriptor && tk.cur().text == "alias") {
            tk.next();
            Type aliasType = parseType(tk);
            tkAssert(tk.curType() == TkDescriptor, "malformed preprocessor block. Expected 'PYTHON alias cname pyname; '");
            string aliasName = tk.cur().text;
            tk.next();
            tkAssert(tk.curType() == TkSemicolon && tk.isLast(), "malformed preprocessor block. Expected 'PYTHON alias cname pyname;'");
            return processPythonInstantiation(block, aliasType, aliasName);
        }
        List<Type> templTypes;

        // resolve template class
        if (tk.curType() == TkTemplate) {
            tk.next();
            templTypes = parseTypeList(tk);            
        }

        if (tk.curType() == TkClass && tk.cur().text != "typename") {
            block.cls.templateTypes = templTypes;
            tk.next();
            tkAssert(tk.curType() == TkDescriptor, "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");
            block.cls.name = tk.cur().text;
            tk.next();
            tkAssert(tk.curType() == TkColon, "PYTHON class must publicly derive from PbClass (or a subclass)");
            tk.next();
            tkAssert(tk.curType() == TkPublic, "PYTHON class must publicly derive from PbClass (or a subclass)");
            tk.next();
            tkAssert(tk.curType() == TkDescriptor, "PYTHON class must publicly derive from PbClass (or a subclass)");
            block.cls.baseClass = parseType(tk);
            tkAssert(tk.curType() == TkCodeBlock && tk.isLast(), "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");

            return processPythonClass(block, tk.cur().text);
        }
        else
        {
            bool isConstructor = tk.curType() == TkDescriptor && gParent == tk.cur().text;
            block.func = parseFunction(tk, false, !isConstructor, false);

            if (isConstructor && tk.curType() == TkColon) {
                // read till end
                block.initList = " : ";
                while(!tk.done() && tk.curType() != TkSemicolon && tk.curType() != TkCodeBlock) {
                    block.initList += tk.cur().text;
                    tk.next();
                }
                tkAssert(!tk.done(), "Constructor initializer list not limited");
            }

            if (tk.curType() == TkSemicolon && block.func.noParentheses) {
                tkAssert(tk.curType() == TkSemicolon && tk.isLast(), 
                    "malformed preprocessor keyword block. Expected 'PYTHON type varname;'");
               return processPythonVariable(block);
            }
            tkAssert((tk.curType() == TkCodeBlock || tk.curType() == TkSemicolon) && tk.isLast(), 
                "malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]'");
            return processPythonFunction(block, tk.cur().text);
        }

    }
    return "";
}
