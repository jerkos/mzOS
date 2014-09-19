/*
emass.cpp: Calculation of accurate masses and 
           intensities of isotopic peaks

Based on an algorithm developed by Alan L. Rockwood.

Published in 
Rockwood, A.L. and Haimi, P.: "Efficent calculation of 
Accurate Masses of Isotopic Peaks",
Journal of The American Society for Mass Spectrometry
JASMS 03-2263, in press

Copyright (c) 2005 Perttu Haimi and Alan L. Rockwood

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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cwchar>
#include <cfloat>

#include "parser.h"
#include "getopt.h"

using std::wcin;
using std::wcout;
using std::wcerr;
using std::endl;
using std::wifstream;
using std::wistringstream;
using std::ios;

using std::string;
using std::wstring;

using std::getline;
using std::swap;

using std::abs;

struct peak
{
  double mass;
  double rel_area;
};

typedef std::vector<peak> Pattern;              // index: peak_number
typedef std::vector<Pattern> SuperAtomList;        // index: bit_number 
typedef std::vector<SuperAtomList> SuperAtomData;  // index: element_number

SuperAtomData sad;
ElemMap em;

const double ELECTRON_MASS = 0.00054858;
const double DUMMY_MASS = -10000000;

int init_data(string filename, SuperAtomData & sad, ElemMap & em)
{
  wifstream f(filename.c_str());
  
  if(f.fail())
    return 0;

  sad.clear();
  em.clear();

  wstring line;
  ulong elemindex = 0;
  int state = 0;
  while(getline(f, line)) {
    wistringstream ist(line);
    wstring element;
    switch(state) {
    case 0: // new element
      ist >> element;
      em[element] = elemindex;
      sad.push_back(SuperAtomList(1));
      sad.back().reserve(8); // reserve room for 8 superatoms
      elemindex++;
      state = 1;
      break;
    case 1: // isotope
      peak p;
      Pattern & idist = sad.back()[0];
      if(ist >> p.mass >> p.rel_area) {
	// fill the gaps in the patterns with zero abundancy peaks
	if(idist.size() > 0) {
	  double prevmass = idist.back().mass;
	  for(int i = 0; i < int(p.mass - prevmass - 0.5); i++) {
	    peak filler;
	    filler.mass = DUMMY_MASS;
	    filler.rel_area = 0;
	    idist.push_back(filler);
	  }
	}
	// insert the peak
	idist.push_back(p);                                                        
      } else  
	state = 0; // no more isotope data
      break;
    }
  }  
  f.close();
  return 1;
}

// Merge two patterns to one.
void convolute_basic(Pattern & h, const Pattern & g, const Pattern & f)
{
  h.clear();
  size_t g_n = g.size();
  size_t f_n = f.size();
  if(g_n == 0 || f_n == 0)
     return;
  for(size_t k = 0; k < g_n + f_n - 1; k++) {
    double sumweight = 0, summass = 0;
    size_t start = k < (f_n - 1) ? 0 : k - f_n + 1; // max(0, k-f_n+1)
    size_t end = k < (g_n - 1) ? k : g_n - 1;       // min(g_n - 1, k)
    for(size_t i = start; i <= end; i++) {
      double weight = g[i].rel_area * f[k - i].rel_area;
      double mass = g[i].mass + f[k - i].mass;
      sumweight += weight;
      summass += weight * mass;
    }
    peak p;
    if(sumweight == 0)
      p.mass = DUMMY_MASS;
    else
      p.mass = summass / sumweight;
    p.rel_area = sumweight;
    h.push_back(p);
  }
}

// Prune the small peaks from both sides but
// leave them within the pattern.
void prune(Pattern & f, double limit)
{
  // prune the front
  Pattern::iterator i = f.begin();
  while(i != f.end()) {
    if((*i).rel_area > limit)
      break;
    i++;
  }
  f.erase(f.begin(), i);

  // prune the end
  while(1) {
    if(f.size() == 0)
      break;
    if(f.back().rel_area > limit)
      break;
    f.pop_back();
  } 
}

//                script: RMS (weighted) (ATGC)1000
void print_pattern(Pattern & result, int digits)
{
  // find the maximum
  double max_area = 0;
  double sum_area = 0;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(max_area < (*i).rel_area)
      max_area = (*i).rel_area;
    sum_area += (*i).rel_area;
  }
  if(max_area == 0)
    return; // empty pattern

  wcout.setf(ios::fixed);
  wcout.precision(digits);
  double print_limit = pow(10.0, -digits) / 2;
  //wcout.precision(30);
  //double print_limit = 0.000001 / 2;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    double mass = (*i).mass;
    double rel_area = (*i).rel_area;
    double val_perc = rel_area / max_area * 100;
    //double val_norm = rel_area / sum_area;
    if(mass != DUMMY_MASS && val_perc >= print_limit)
      wcout << mass << L" " << val_perc << endl;
    //wcout << mass << L" " << val_perc << L" " << val_norm << endl;
  }
}

void calculate(Pattern & tmp, Pattern & result, FormMap & fm, 
	       double limit, long charge)
{
  for(FormMap::iterator i = fm.begin(); i != fm.end(); i++) {
    size_t atom_index = (*i).first;
    SuperAtomList sal = sad[atom_index];
    ulong n = (*i).second;
    ulong j = 0;
    while(n > 0) {
      size_t sz = sal.size();
      if(j == sz) { 
	sal.resize(sz + 1); // Make new superatom from previous
                            // largest superatom. We are trying to
                            // avoid copying on assignment here.
	convolute_basic(sal[j], sal[j - 1], sal[j - 1]);
	prune(sal[j], limit);
      }
      if(n & 1) { // digit is 1, convolute result
	convolute_basic(tmp , result, sal[j]);
	prune(tmp, limit);
	swap(tmp, result);  // Hopefully the swap implementation
                            // will not copy all elements.
      }
      n >>= 1; 
      j++;
    }
  }
  
  // take charge into account
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(charge > 0)
      (*i).mass = (*i).mass / abs(charge) - ELECTRON_MASS;
    else if (charge < 0)
      (*i).mass = (*i).mass / abs(charge) + ELECTRON_MASS;
  }
}

int main(int argc, char * argv[])
{
  char * isotopefn = "ISOTOPE.DAT";
  double limit = 0;
  int digits = 6;
  int opt;
  while((opt = getopt(argc, argv, "i:l:d:")) != -1) {
    if(opt == 'i')
      isotopefn = optarg;
    else if(opt == 'l')
      limit = atof(optarg);
    else if(opt == 'd') {
      digits = atoi(optarg);
      if(digits < 1)
	digits = 1;
      if(digits > 30)
	digits = 30;
    } else
      return 2;
  }

  if(!init_data(isotopefn, sad, em)) {
    string tmp(isotopefn);
    wcerr << "Could not open isotope file: " 
	  << wstring(tmp.begin(), tmp.end()) << endl;
    return 1;
  }

  Parser p(em);

  FormMap fm;

  Pattern result(1);
  Pattern tmp;

  wstring line;
  while(getline(wcin, line)) {

    size_t sep_pos = line.find(L',');
    long charge = 0;
    if(sep_pos < line.size() - 1) {
      wstring charge_str = wstring(line, sep_pos + 1, line.size());
      charge = wcstol(charge_str.c_str(), NULL, 0);
    }

    wstring formula(line, 0, sep_pos);
    wcout << L"formula: " << formula 
	  << L" charge : " << charge 
	  << L" limit: "
	  << std::setiosflags(ios::scientific) 
	  << limit 
	  << std::resetiosflags(ios::scientific)
	  << endl;

    // read a formula
    fm.clear();
    try {
      p.write_compound(formula, fm);
    } catch (Parser::Error e) {
      wcerr << e._message;
      continue; // skip to next compound
    }

    // initialize the result
    result.resize(1);
    result.front().mass = 0.0;
    result.front().rel_area = 1.0;

    calculate(tmp, result, fm, limit, charge);
    print_pattern(result, digits);
    
  }

  return 0;
}
