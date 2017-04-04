/*========================================================================================================
                                                   random.h
==========================================================================================================

 Routines for the generation of pseudo-random numbers
 
 C++-code accompanying:
 
        authors and G. Sander van Doorn
        title
 
 Written by:
        G. Sander van Doorn
        Groningen Institute for Evolutionary Life Sciences (Gelifes)
        University of Groningen
        the Netherlands
 
 Program version
        xx/xx/20xx	: ...
 
 =====================================================================================================*/

#ifndef random_h
#define random_h

#include <random>
#include "utils.h"

namespace rnd
{ 
    void set_seed();
	void set_seed(unsigned int);
	extern std::mt19937_64 rng;
	int integer(int);
	size_t integer(size_t);
	bool bernoulli(double = 0.5);
	int binomial(int, double = 0.5);
	size_t binomial(size_t, double = 0.5);
	size_t poisson(double = 1.0);
	double uniform();
	double normal(double = 0.0, double = 1.0);
	double exponential(double = 1.0);
}

#endif //#ifndef random_h