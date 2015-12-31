/*
formula.cpp: minimalistic representation of molecular
             formulas

Copyright (c) 2005 Perttu Haimi

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

#include "formula.h"
#include <vector>

/*
 * Removes elements having zero count.
 */
void cleanUpFormula(FormMap & form) {
  
  typedef std::vector<size_t> Erasable;

  Erasable e;
  for(FormMap::iterator i = form.begin(); i != form.end(); ++i)
    if(i->second == 0)
      e.push_back(i->first);

  for(Erasable::iterator i = e.begin(); i != e.end(); ++i)
      form.erase(*i);
}


/*
 * Substracts FormMap sub from FormMap form.
 */
void subtractFormula(FormMap & form, const FormMap & sub) {
  FormMap tmp = sub;
  negateFormula(tmp);
  addFormula(form, tmp);
}


/*
 * Adds each element count from FormMap add to FormMap form.
 */
void addFormula(FormMap & form, const FormMap & add) {
  for(FormMap::const_iterator i = add.begin(); i != add.end(); ++i)
    form[i->first] += i->second;
  cleanUpFormula(form);
}


/*
 * Negates the sign of the count of each element in the formula.
 */
void negateFormula(FormMap & form) {
  for(FormMap::iterator i = form.begin(); i != form.end(); ++i)
    i->second = -(i->second);
}


/*
 * Checks that all atom counts are positive.
 */
bool isRealFormula(const FormMap & form) {
  for(FormMap::const_iterator i = form.begin(); i != form.end(); ++i) {
    if(i->second < 0)
      return false;
  }
  return true;
}
