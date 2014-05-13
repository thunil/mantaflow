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

vector<Type> parseTypeList(TokenPointer&);

#define tkAssert(x,msg) if(!(x)){errMsg(tk.curLine(),tk.completeText+": "+msg);}
#define typeAssert(x) tkAssert(x," is not a valid type.")
#define argAssert(x) tkAssert(x," is not a valid argument.")

Type parseType(TokenPointer& parentPtr) {
    Type cur = Type();
    TokenPointer tk(parentPtr);

    // constness
    if (tk.curType() == TkConst) {
        cur.isConst = true;
        tk.next();
    }
    typeAssert(!tk.done());
    
    // signed / unsigned
    if (tk.curType() == TkTypeQualifier) {
        cur.name = tk.cur().text;
        if (!tk.next())
            return cur;
        typeAssert(tk.curType() == TkSimpleType);
        cur.name += " " + tk.cur().text;
        tk.next();
        parsePointer(tk, cur);
        return cur;
    }

    typeAssert(tk.curType() == TkDescriptor || tk.curType() == TkSimpleType || tk.curType() == TkClass);
    cur.name = tk.cur().text;
    tk.next();

    // namespace
    if (tk.curType() == TkDoubleColon) {
        cur.name += "::";
        tk.next();
        typeAssert(tk.curType() == TkDescriptor);
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
    TokenPointer tk(parentPtr);

    if (requireType)
        cur.type = parseType(tk);

    if (tk.curType() != TkDescriptor && !requireName) 
        return cur;
    argAssert(tk.curType() == TkDescriptor);
    cur.name = tk.cur().text;
    tk.next();

    // default value ?
    if (tk.curType() == TkAssign) {
        tk.next();
        argAssert(tk.curType() == TkDescriptor);
        tk.next();

        // don't validate, just track bracket level
        BracketStack stack;
        for(;!tk.done();tk.next()) {
            if (tk.curType() == TkBracketL || tk.curType() == TkTBracketL) {
                stack.push_back(tk.cur().text[0]);
            } 
            else if (stack.empty() && (tk.cur().type == TkComma || 
                                       tk.curType() == TkBracketR || tk.curType() == TkTBracketR)) {
                tk.previous();
                break;
            } else if (tk.curType() == TkBracketR) {
                argAssert(stack.pop() == '(');
            } else if (tk.curType() == TkTBracketR) {
                argAssert(stack.pop() == '<');
            }
        }
        argAssert(stack.empty());
    }
    return cur;
}

vector<Type> parseTypeList(TokenPointer& parentPtr) {
    TokenPointer tk(parentPtr);
    vector<Type> list;
    
    typeAssert(tk.curType() == TkTBracketL);
    tk.next();
    for(;;) {
        list.push_back(parseType(tk));
        if (tk.curType() == TkTBracketR) 
            break;
        typeAssert(tk.curType() == TkComma);
    }
    typeAssert(tk.curType() == TkTBracketR);
    tk.next();
    return list;
}

vector<Argument> parseArgumentList(TokenPointer& parentPtr, bool requireName, bool requireType) {
    TokenPointer tk(parentPtr);
    vector<Argument> list;
    
    typeAssert(tk.curType() == TkBracketL);
    tk.next();
    for(;;) {
        list.push_back(parseArgument(tk, requireName, requireType));
        if (tk.curType() == TkBracketR) 
            break;
        typeAssert(tk.curType() == TkComma);
    }
    typeAssert(tk.curType() == TkBracketR);
    tk.next();
    return list;
}

Function parseFunction(TokenPointer& parentPtr, bool requireNames, bool requireType, bool requireArgs) {
    TokenPointer tk(parentPtr);
    Function cur;

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

    if (tk.curType() == TkConst)
        cur.isConst = true;

    return cur;
}

// Parse syntax KEYWORD(opt1, opt2, ...) STATEMENTS [ {} or ; ]    
string parseBlock(const string& kw, const vector<Token>& tokens) {    
    TokenPointer tk(tokens);

    // parse keyword options
    ArgList options = parseArgumentList(tk, true, false);
    
    if (kw == "KERNEL") {
        ArgList returnArgs;
        std::vector<Type> templTypes;

        // templated kernel
        if (tk.curType() == TkTemplate) {
            tk.next();
            templTypes = parseTypeList(tk);            
        }
        
        // return values
        while (tk.curType() == TkDescriptor && tk.cur().text == "returns") {
            tk.next();
            typeAssert(tk.curType() == TkBracketL);
            tk.next();
            returnArgs.push_back(parseArgument(tk, true, true));
            typeAssert(tk.curType() == TkBracketR);
            tk.next();            
        }

        Function kernel = parseFunction(tk, true, true, true);
        if (!templTypes.empty())
            kernel.templateTypes = templTypes;

        tkAssert(tk.curType() == TkCodeBlock && tk.isLast(), 
            "Malformed KERNEL, expected KERNEL(opts...) ret_type name(args...) { code }");

        // todo procKernel
    }
    else if (kw == "PYTHON")
    {
        // template instantiation / alias
        if (tk.curType() == TkDescriptor && tk.cur().text == "alias") {

            /*

            lb += consumeWS(tokens, ++index);
            assert (tokens[index].type == TkDescriptor, "syntax error. expected PYTHON alias type<T> name;");
            string classname = tokens[index++].text;
            
            templArgs = parseArgs(tokens, index, false, lb, true);
            
            assert (tokens[index].type == TkDescriptor, "syntax error. expected PYTHON alias type<T> name;");
            string aliasname = tokens[index++].text;
            lb += consumeWS(tokens, index);
            
            if (tokens[index].type != TkSemicolon || index+1 != tokens.size())
                errMsg(line, "syntax error. expected PYTHON alias type<T> name;");            
            return processPythonInstantiation(lb, classname, templArgs, aliasname, line);*/
        }
        vector<Type> templTypes;

        // resolve template class
        if (tk.curType() == TkTemplate) {
            tk.next();
            templTypes = parseTypeList(tk);            
        }

        if (tk.curType() == TkClass && tk.cur().text != "typename") {
            Class cur;
            cur.templateTypes = templTypes;
            tk.next();
            tkAssert(tk.curType() == TkDescriptor, "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");
            cur.name = tk.cur().text;
            tk.next();
            tkAssert(tk.curType() == TkColon, "PYTHON class must publicly derive from PbClass (or a subclass)");
            tk.next();
            tkAssert(tk.curType() == TkPublic, "PYTHON class must publicly derive from PbClass (or a subclass)");
            tk.next();
            tkAssert(tk.curType() == TkDescriptor, "PYTHON class must publicly derive from PbClass (or a subclass)");
            cur.baseClass = parseType(tk);
            tk.next();
            tkAssert(tk.curType() == TkCodeBlock && tk.isLast(), "malformed preprocessor keyword block. Expected 'PYTHON class name : public X {}'");

            // todo proc class
        }
        else
        {
            bool isConstructor = tk.curType() == TkDescriptor && gParent == tk.cur().text;
            Function cur = parseFunction(tk, false, !isConstructor, false);

            if (tk.curType() == TkSemicolon && cur.noParentheses) {
                // proc python variable                
            }
            tkAssert((tk.curType() == TkCodeBlock || tk.curType() == TkSemicolon) && tk.isLast(), 
                "malformed preprocessor keyword block. Expected 'PYTHON type funcname(args) [{}|;]'");

            // proc member func
        }

    }
    return "";
}
