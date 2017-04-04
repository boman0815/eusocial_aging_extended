/*========================================================================================================
                                                random.cpp
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

========================================================================================================*/

#include <chrono>
#include <sstream>
#include <iomanip>
#include "random.h"

namespace rnd {
    std::mt19937_64 rng;
	
    void set_seed(unsigned seed)
    {
        std::ostringstream oss;
        oss << "random seed set to " << seed;
        echo(oss.str());
        rng.seed(seed);
    }
    
	void set_seed()
    {
        unsigned seed = static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        set_seed(seed);
    }

	int integer(int n) { return std::uniform_int_distribution<int>(0, --n)(rng); }
	size_t integer(size_t n) { return std::uniform_int_distribution<size_t>(0u, --n)(rng); }
	bool bernoulli(double p) { return std::bernoulli_distribution(p)(rng); }
	int binomial(int n, const double p) { return std::binomial_distribution<int>(n, p)(rng); }
	size_t binomial(size_t n, const double p) { return std::binomial_distribution<size_t>(n, p)(rng); }
	size_t poisson(double lambda) { return std::poisson_distribution<size_t>(lambda)(rng); }
	double uniform() { return std::uniform_real_distribution<double>(0.0, 1.0)(rng); }
	double normal(double mu, double sigma) { return std::normal_distribution<double>(mu, sigma)(rng); }
	double exponential(double lambda) { return std::exponential_distribution<double>(lambda)(rng); }
}