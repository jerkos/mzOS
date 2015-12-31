#if !defined(FORMULA_H)
#define FORMULA_H

#include <map>
#include <string>

// map from element abbreviation to index in the elements table
typedef std::map<std::wstring, size_t> ElemMap;

// map from element index to the count of occurences in the formula
typedef std::map<size_t, long> FormMap;

void cleanUpFormula(FormMap & form);

void subtractFormula(FormMap & form, const FormMap & sub);

void addFormula(FormMap & form, const FormMap & add);

void negateFormula(FormMap & form);

bool isRealFormula(const FormMap & form);

#endif // FORMULA_H
