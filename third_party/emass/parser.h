// parser.h: interface for the Parser class.
//
// Copyright 2005 Perttu Haimi
//

#if !defined(PARSER_H)
#define PARSER_H

#include <list>
#include <string>

#include "formula.h"

typedef unsigned long ulong;

class Parser  
{
public:
  Parser(ElemMap const & e);
  
  void parse_formula(const std::wstring & frml);
  void write_compound(const std::wstring & frml, FormMap & c);
  
  // Exception class
  class Error
    {
    public:
      Error(std::wstring message)
	: _message(message) {};
      std::wstring _message;
    };
  
 private:
  enum TokenType {LPAREN, RPAREN, EOS, NUM, ELEMENT};
  
  // Helper class for handling the input std::wstring
  class Tokenizer
    {
    public:
      Tokenizer()
	: _ttype(EOS), _tval(L""), _i(0), _last(0) {};
      void get_token();
      void set_input(const std::wstring & input);
      void error(std::wstring msg);
      
      // Current tokentype and value
      TokenType _ttype;
      std::wstring _tval;
    private:
      enum States {START, ALPHA, LITERAL, NUMBER, END, STOP};
      unsigned long _i;
      unsigned long _last;
      std::wstring _input;
    };
  friend class Parser::Tokenizer;
  
  // Helper classes for forming the parse tree
  class Node
    {
    public:
      virtual ~Node() {};
      virtual void fill_compound(ulong weight, FormMap & c) const = 0;
    };
  
  typedef std::list<Node *> NodeList;
  
  class ElementNode : public Node
    {
    public:
      ElementNode() {};
      ElementNode(size_t elem_index)
	: _elem_index(elem_index) {};
      void fill_compound(ulong weight, FormMap & c) const 
	{ c[_elem_index] += weight; };
    private:
      size_t _elem_index;
    };
  
  class ElementListNode : public Node
    {
    public:
      size_t length();
      ElementListNode() 
	: _count(1) {};
      ElementListNode(size_t elem_index);
      ~ElementListNode();
      void fill_compound(ulong weight, FormMap & c) const;
      void append(Node * n) { _NodeList.push_front(n); };
      void set_count(ulong count) { _count = count; };
    private:
      NodeList _NodeList;
      ulong _count;
    };
  
  // Private member function
  Parser::ElementListNode * parse();	
  
  // Private data
  ElemMap const & _elements;
  Parser::Tokenizer _tokenizer;
};

#endif // PARSER_H
