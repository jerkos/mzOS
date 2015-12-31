/*
Parser.cpp: parser for molecular formulas

Copyright (c) 2004, 2005 Perttu Haimi

All rights reserved.

Redistribution and use in source and binary forms,
with or without modification, are permitted provided
that the following conditions are met:

    * Redistributions of source code must retain the
      above copyright notice, this list of conditions
      and the following disclaimer.
    * Redistributions in binary form must reproduce
      the above copyright notice, this list of conditions
      and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors
      may be used to endorse or promote products derived
      from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cwctype> // isnumber() etc.
#include <cwchar> // wcstol() wide chars to numeric conversion

#include "parser.h"
#include <cassert>

using std::wstring;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Parser::Parser(ElemMap const & e) 
  : _elements(e) 
{};

// Will throw exception if frml in not correct.
void Parser::parse_formula(const wstring & frml)
{
  _tokenizer.set_input(frml);
  _tokenizer.get_token();
  Parser::ElementListNode * tmp = parse();
  
  if(_tokenizer._ttype != EOS) {
    delete tmp; 
    _tokenizer.error(L"End of input expected");
  }
  delete tmp;
}

void Parser::write_compound(const wstring & frml, FormMap & c)
{
  c.clear();

  _tokenizer.set_input(frml);
  _tokenizer.get_token();
  Parser::ElementListNode * tmp = parse();
  
  if(_tokenizer._ttype != EOS) {
    delete tmp;	
    _tokenizer.error(L"End of input expected");
  }

  tmp->fill_compound(1, c);
  cleanUpFormula(c);
  delete tmp;	
}

// Parse a chemical formula.
// Formula can contain parenthesis.
Parser::ElementListNode * Parser::parse()
{
  Parser::ElementListNode * inner;
  Parser::ElementListNode * outer = new Parser::ElementListNode();
  while(_tokenizer._ttype == LPAREN ||
	_tokenizer._ttype == ELEMENT) {
    if(_tokenizer._ttype == LPAREN) {
      try {
	_tokenizer.get_token();
      } catch (Parser::Error e) {
	delete outer;
	throw e;
      }
      inner = parse();
      if(_tokenizer._ttype != RPAREN) {
	delete outer;
	delete inner;
	_tokenizer.error(L"Expected right parenthesis");
      }
    } else {
      ElemMap::const_iterator i = _elements.find(_tokenizer._tval);
      if(i != _elements.end())
	inner = new ElementListNode((*i).second);
      else {
	delete outer;
	_tokenizer.error(L"\'" + _tokenizer._tval + L"\' is not an element");
      }
    }
    _tokenizer.get_token();
    if(_tokenizer._ttype == NUM) {
      inner->set_count(wcstol(_tokenizer._tval.c_str(), NULL, 0));
      _tokenizer.get_token();
    }
    outer->append(inner);
  }
  return outer;
}

// how to make this a member?
static const wstring ENDTOKEN(L"%END%");

void Parser::Tokenizer::error(wstring msg)
{
  msg += L":\n";
  _input.resize(_input.length() - ENDTOKEN.length());
  msg += _input + L"\n";
  msg += wstring(_last, L' ') + L"^\n";
  throw Parser::Error(msg);
}

void Parser::Tokenizer::set_input(const wstring & input)
{
  _i = 0;
  _input = input;
  _input += ENDTOKEN;
}

void Parser::Tokenizer::get_token()
{
  States s = START;
  _last = _i;
  for(;;)
    {
      wchar_t c = _input[_i];
      switch(s) {
      case START:
	if(iswspace(c)) {
	  ; // skip whitespace
	} else if(iswdigit(c) || c == L'-') {
	  _ttype = NUM;	
	  _tval = c;
	  s = NUMBER;
	} else if(iswupper(c)) {
	  _ttype = ELEMENT;
	  _tval = c;
	  s = ALPHA;
	} else if(c == L'[') {
	  _ttype = ELEMENT;
	  _tval = L"";
	  s = LITERAL;
	} else if(c == L'(') {
	  _ttype = LPAREN;
	  s = STOP;
	} else if(c == L')') {
	  _ttype = RPAREN;
	  s = STOP;
	} else if(c == L'%') {
	  _ttype = EOS;
	  s = END;
	} else {
	  error(L"Parse error");
	}
	break;
      case ALPHA:
	if(iswlower(c)) {
	  _tval += c;
	} else {
	  return; // we're done
	}
	break;
      case LITERAL:
	if(c == L'%') {
	  error(L"Missing left bracket");
	} else if(c != L']') {
	  _tval += c;
	} else {
	  ++_i; // skip the ']'
	  return; // we're done
	}	
	break;
      case NUMBER:
	if(iswdigit(c)) {
	  _tval += c;
	} else {
	  return; // we're done
	}
	break;
      case END:
	if(_input.find(ENDTOKEN) == _i - 1) {
	  s = STOP;
	} else {
	  error(L"Illegal character: \'%\'");
	}
	break;
      case STOP:
	return; // exit
      }
      ++_i;
    }
}


Parser::ElementListNode::ElementListNode(size_t elem_index)
  : _count(1)	 
{
  ElementNode * en = new ElementNode(elem_index);
  _NodeList.push_front(en);
}


Parser::ElementListNode::~ElementListNode()
{
  for(NodeList::iterator i = _NodeList.begin();
      i != _NodeList.end();
      ++i)
    {
      delete *i;
    }
}

void Parser::ElementListNode::fill_compound(ulong weight, FormMap & c) const
{
  ulong totalweight = _count * weight;
  
  for(NodeList::const_iterator i = _NodeList.begin();
      i != _NodeList.end();
      ++i)
    {
      (*i)->fill_compound(totalweight, c);
    }

}

size_t Parser::ElementListNode::length()
{
  return _NodeList.size();
}

