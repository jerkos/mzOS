#if !defined(EMASS_H)
#define EMASS_H

#include <vector>
#include <string>

#include "formula.h"

struct peak
{
  double mass;
  double rel_area;
};

typedef std::vector<peak> Pattern;              // index: peak_number
typedef std::vector<Pattern> SuperAtomList;        // index: bit_number 
typedef std::vector<SuperAtomList> SuperAtomData;  // index: element_number

typedef Pattern::iterator pit;
typedef Pattern::const_iterator cpit;

void init_data(std::string filename, SuperAtomData & sad, ElemMap & em);
void convolute_basic(pit h, cpit g, size_t ng, 
		     cpit f, size_t nf);
void convolute_karat(pit h, cpit g, size_t ng, 
		     cpit f, size_t nf);
void merge(pit res, cpit a, size_t na, 
	   cpit b, size_t nb, pit sum, size_t nsum);
void calc_sum(cpit a, size_t shift, pit res, size_t na);
void convolute_lopsided(pit h,
			cpit g, size_t ng,
			cpit f, size_t nf);
void prune(Pattern & f, double limit);
void normalize(Pattern & f);
void calculate(Pattern & tmp, Pattern & result, FormMap & fm, 
	       double limit, long charge, SuperAtomData sad);
void print_pattern(Pattern & result, int digits);

#endif // EMASS_H
