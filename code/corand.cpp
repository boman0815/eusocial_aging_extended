/*========================================================================================================
                                                   corand.cpp
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

#include "corand.h"

namespace rnd {
    MultiNormal::MultiNormal(const Matrix &Sigma, const Vector &mean) : n(mean.size()), mu(mean)
    {
        verify(n == Sigma.sizec() && n == Sigma.sizer());
        std::vector<Vector> tmp;
        sigma = Sigma.eigenvaluesSymPart(tmp);
        for(double &x : sigma) {
            if(x < 0.0) {
                x = 0.0;
                warning(CURRENT_FUNCTION, "variance/covariance matrix is not positive definite");
            }
            x = sqrt(x);
        }
        U = Matrix(tmp);
    }
    MultiNormal::MultiNormal(const Matrix &Sigma) : n(Sigma.sizec()), mu(Vector(Sigma.sizec(), 0.0))
    {
        verify(n == Sigma.sizer());
        std::vector<Vector> tmp;
        sigma = Sigma.eigenvaluesSymPart(tmp);
        for(double &x : sigma) {
            if(x < 0.0) x = 0.0;
            x = sqrt(x);
        }
        U = Matrix(tmp);
    }
    Vector MultiNormal::operator()() const
    {
        Vector x(n);
        for(Vector::size_type i = 0u; i < n; ++i) x[i] = rnd::normal(0.0, sigma[i]);
        return mu + U * x;
    }
    std::ostream& operator<< (std::ostream &os, const MultiNormal &obj)
    {
        os << obj.sigma << '\n';
        os << '\n' << obj.U;
        return os;
    }
}

