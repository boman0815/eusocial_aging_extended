/*========================================================================================================
                                                   corand.h
==========================================================================================================

 Routines for the generation of correlated pseudo-random numbers
 
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

#ifndef corand_h
#define corand_h

#include "random.h"
#include "linalg.h"

namespace rnd
{
    class MultiNormal
    {
    public:
        MultiNormal(const Matrix&, const Vector&);
        MultiNormal(const Matrix&);
        Vector operator()() const;
        friend std::ostream& operator<< (std::ostream&, const MultiNormal&);
    private:
        const Vector::size_type n;
        const Vector mu;
        Vector sigma;
        Matrix U;
    };
}
#endif /* #ifndef corand_h */
