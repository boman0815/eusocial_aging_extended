/*========================================================================================================
 linalg.h
 ==========================================================================================================
 
 Vector and matrix algebra
 
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

#include "linalg.h"
#include "utils.h"

/*=======================================================================================================
 implementation of the class Vector
 ========================================================================================================*/

Vector& Vector::operator=(const Vector &other)
{
    if(this != &other) {
        size_type n = other.size();
        resize(n);
        double *p = data();
        double const *q = other.data();
        while(n) {
            *p = *q;
            ++p; ++q;
            --n;
        }
    }
    return *this;
}

Vector& Vector::operator+= (const Vector &other)
{
    size_type n = size();
    verify(n == other.size());
    double *p = data();
    double const *q = other.data();
    while(n) {
        *p += *q;
        ++p; ++q;
        --n;
    }
    return *this;
}

Vector& Vector::operator-= (const Vector &other)
{
    size_type n = size();
    verify(n == other.size());
    double *p = data();
    double const *q = other.data();
    while(n) {
        *p -= *q;
        ++p; ++q;
        --n;
    }
    return *this;
}

double Vector::operator*(const Vector &other) const
{
    size_type n = size();
    verify(n == other.size());
    double const *p = data(), *q = other.data();
    double innerProduct = 0.0;
    while(n) {
        innerProduct += (*p) * (*q);
        ++p; ++q;
        --n;
    }
    return innerProduct;
}

Vector& Vector::operator*=(const double &factor)
{
    for(std::vector<double>::iterator it = begin(); it != end(); ++it) *it *= factor;
    return *this;
}

double Vector::total() const
{
    double sum = 0.0;
    for(std::vector<double>::const_iterator it = begin(); it != end(); ++it) sum += *it;
    return sum;
}

std::ostream& operator<< (std::ostream &os, const Vector &vec)
{
    bool fill = false;
    for(const double &x : vec) {
        fill ? os << os.fill() << x : os << x;
        fill = true;
    }
    return os;
}

/*=======================================================================================================
 implementation of the class Matrix
 ========================================================================================================*/

Matrix::Matrix(const std::vector<std::vector<double> > &mat) :
szr(mat.size()), szc(mat.front().size()), isLUdecomp(false)
{
    resize(szr * szc);
    double *p = data();
    for(const std::vector<double> &row : mat)
        for(const double &x : row) {
            *p = x;
            ++p;
        }
}

Matrix::Matrix(const std::vector<Vector> &mat) :
szc(mat.size()), szr(mat.front().size()), isLUdecomp(false)
{
    resize(szr * szc);
    double *p = data();
    for(size_type i = 0u; i < szr; ++i)
        for(size_type j = 0u; j < szc; ++j, ++p) *p = mat[j][i];
}

Matrix::Matrix(const Vector &diag) :
szr(diag.size()), szc(diag.size()), std::vector<double>(diag.size() * diag.size(), 0.0), isLUdecomp(false)
{
    setDiag(diag);
}

Matrix::~Matrix()
{
    if(isLUdecomp) {
        delete [] indx;
        for(size_type i = 0u; i < szr; ++i) delete [] LU[i];
        delete [] LU;
    }
}

Matrix& Matrix::operator=(const Matrix &other)
{
    if(this != &other) {
        szc = other.szc;
        szr = other.szr;
        size_type n = szc * szr;
        resize(n);
        double *p = data();
        double const *q = other.data();
        while(n) {
            *p = *q;
            ++p; ++q;
            --n;
        }
    }
    return *this;
}

Matrix& Matrix::operator+= (const Matrix &other)
{
    verify(szc == other.szc && szr == other.szr);
    size_type n = size();
    double *p = data();
    double const *q = other.data();
    while(n) {
        *p += *q;
        ++p; ++q;
        --n;
    }
    return *this;
}

Matrix& Matrix::operator-= (const Matrix &other)
{
    verify(szc == other.szc && szr == other.szr);
    size_type n = size();
    double *p = data();
    double const *q = other.data();
    while(n) {
        *p -= *q;
        ++p; ++q;
        --n;
    }
    return *this;
}

Matrix& Matrix::operator*=(const double &factor)
{
    for(std::vector<double>::iterator it = begin(); it != end(); ++it) *it *= factor;
    return *this;
}

Matrix Matrix::transpose() const
{
    Matrix result(szc, szr);
    double const *p = data();
    double *q = result.data();
    for(size_type i = 0u; i < szr; ++i)
        for(size_type j = 0u; j < szc; ++j, ++p) {
            size_type k = j * szr + i;
            *(q + k) = *p;
        }
    return result;
}

Vector Matrix::operator*(const Vector &x) const
{
    verify(szc == x.size());
    Vector result(szr, 0.0);
    double const *p = data();
    double *r = result.data();
    for(size_type i = 0u; i < szr; ++i, ++r) {
        double const *q = x.data();
        for(size_type j = 0u; j < szc; ++j, ++p, ++q) *r += (*p) * (*q);
    }
    return result;
}

Vector operator*(const Vector &x, const Matrix &M)
{
    verify(M.szr == x.size());
    const Matrix::size_type n = M.szc;
    Vector result(n, 0.0);
    double *r = result.data();
    for(Matrix::size_type j = 0u; j < n; ++j, ++r) {
        double const *q = x.data();
        double const *p = M.data() + j;
        for(Matrix::size_type i = 0u; i < M.szr; ++i, ++q, p += n) *r += (*q) * (*p);
    }
    return result;
}

Matrix Matrix::operator*(const Matrix &other) const
{
    verify(szc == other.szr);
    Matrix result(szr, other.szc, 0.0);
    double *r = result.data();
    for(size_type i = 0u; i < szr; ++i)
        for(size_type j = 0u; j < other.szc; ++j, ++r)
        {
            double const *q = other.data() + j, *p = data() + i * szc;
            for(size_type k = 0u; k < szc; ++k, ++p, q += other.szc) *r += (*p) * (*q);
        }
    return result;
}

void Matrix::setRow(const size_type &i, const Vector &vec)
{
    verify(szc == vec.size());
    double *p = data() + i * szc;
    const double *q = vec.data();
    for(size_type j = 0u; j < szc; ++j, ++p, ++q) *p = *q;
}

void Matrix::setCol(const size_type &j, const Vector &vec)
{
    verify(szr == vec.size());
    double *p = data() + j;
    const double *q = vec.data();
    for(size_type i = 0u; i < szr; ++i, p += szc, ++q) *p = *q;
}

void Matrix::setDiag(const Vector &vec)
{
    verify(szr == szc && szc == vec.size());
    double *p = data();
    const double *q = vec.data();
    for(size_type i = 0u; i < szr; ++i, p += szc + 1, ++q) *p = *q;
}

Vector Matrix::getRow(const size_type &i) const
{
    Vector vec(szc);
    const double *p = data() + i * szc;
    double *q = vec.data();
    for(size_type j = 0u; j < szc; ++j, ++p, ++q) *q = *p;
    return vec;
}

Vector Matrix::getCol(const size_type &j) const
{
    Vector vec(szr);
    const double *p = data() + j;
    double *q = vec.data();
    for(size_type i = 0u; i < szr; ++i, p += szc, ++q) *q = *p;
    return vec;
}

std::vector<Vector> Matrix::getRows() const
{
    std::vector<Vector> vecs(szr, Vector(szc));
    const double *p = data();
    for(size_type i = 0u; i < szr; ++i)
        for(size_type j = 0u; j < szc; ++j, ++p) vecs[i][j] = *p;
    return vecs;
}

std::vector<Vector> Matrix::getCols() const
{
    std::vector<Vector> vecs(szc, Vector(szr));
    const double *p = data();
    for(size_type i = 0u; i < szr; ++i)
        for(size_type j = 0u; j < szc; ++j, ++p) vecs[j][i] = *p;
    return vecs;
}

Vector Matrix::getDiag() const
{
    verify(szr == szc);
    Vector vec(szr);
    const double *p = data();
    double *q = vec.data();
    for(size_type i = 0u; i < szr; ++i, p += szc + 1, ++q) *q = *p;
    return vec;
}


std::ostream& operator<< (std::ostream &os, const Matrix &mat)
{
    Matrix::size_type k = 0u;
    for(const double &x : mat) {
        if(k % mat.szc == 0u) os << x;
        else os << os.fill() << x;
        ++k;
        if(k % mat.szc == 0u && k % mat.size() != 0u) os << '\n';
    }
    return os;
}

Vector Matrix::eigenvaluesSymPart(std::vector<Vector> &V) const
//Computes all eigenvalues and eigenvectors of the symmetric part A of the matrix. The function
//returns a vector containing the eigenvalues of A. V is a matrix whose columns contain, on output,
//the normalized eigenvectors of A.
{
    Matrix A(this->symPart());
    
    //create a nxn indentity matrix
    V = std::vector<Vector>(szc, Vector(szr, 0.0));
    for (size_type ip = 0u; ip < szc; ++ip) V[ip][ip] = 1.0;
    
    //initialize b and d to the diagonal of A; the vector z will accumulate terms of the form tapq as in equation (11.1.14) in Numerical Recipes.
    Vector b(A.getDiag()), d(b), z(szc, 0.0);
    
    int nrot = 0;
    for (int i = 1; i <= 50; ++i)
    {
        double sm = 0.0;
        for (size_type ip = 0u; ip < szc - 1u; ++ip)
        { //Sum off-diagonal elements.
            for (size_type iq = ip + 1u; iq < szc; ++iq)
                sm += fabs(A(ip,iq));
        }
        if ((float)sm == 0.0) //The normal return, which relies on quadratic convergence to machine underflow.
        {
            //normalize eigenvectors
            for (size_type ip = 0u; ip < szc; ++ip) {
                if(V[ip][0] < 0.0) V[ip] *= -1.0;
                V[ip].normalize();
            }
            //sort eigenvalues
            for (size_type ip = 0u; ip < szc - 1u; ++ip)
                for (size_type iq = ip + 1u; iq < szc; ++iq)
                    if(d[iq] > d[ip]) {
                        double tmp = d[ip];
                        d[ip] = d[iq];
                        d[iq] = tmp;
                        Vector tmpvec = V[ip];
                        V[ip] = V[iq];
                        V[iq] = tmpvec;
                    }
            return d;
        }
        double tresh = i < 4 ? 0.2 * sm / (szc * szc) : 0.0;
        for (size_type ip = 0u; ip < szc - 1u; ++ip)
        {
            for (size_type iq = ip + 1u; iq < szc; ++iq)
            {
                double g = 100.0 * fabs(A(ip, iq));
                //After four sweeps, skip the rotation if the off-diagonal element is small.
                if (i > 4 && (float)(fabs(d[ip]) + g) == (float)fabs(d[ip])
                    && (float)(fabs(d[iq]) + g) == (float)fabs(d[iq]))
                    A(ip, iq) = 0.0;
                else
                    if (fabs(A(ip,iq)) > tresh)
                    {
                        double h = d[iq] - d[ip], t;
                        if ((float)(fabs(h) + g) == (float)fabs(h)) t = (A(ip, iq)) / h;
                        else
                        {
                            double theta = 0.5 * h / (A(ip,iq)); //Equation (11.1.10) in Numerical Recipes.
                            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                            if (theta < 0.0) t = -t;
                        }
                        double c = 1.0 / sqrt(1 + t * t);
                        double s= t * c;
                        double tau = s /(1.0 + c);
                        h = t * A(ip, iq);
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        A(ip, iq) = 0.0;
                        for (size_type j = 0u; j < ip; ++j)
                        { //Case of rotations 1 ≤ j < p.
                            g=A(j,ip);
                            h=A(j,iq);
                            A(j,ip) = g - s * (h + g * tau);
                            A(j,iq) = h + s * (g - h * tau);
                        }
                        for (size_type j = ip + 1u; j < iq; ++j)
                        { //Case of rotations p < j < q.
                            g = A(ip,j);
                            h = A(j,iq);
                            A(ip,j) = g - s * (h + g * tau);
                            A(j,iq) = h + s * (g - h * tau);
                        }
                        for (size_type j = iq + 1u; j < szc; ++j)
                        { //Case of rotations q < j ≤ n.
                            g = A(ip,j);
                            h = A(iq,j);
                            A(ip,j) = g - s * (h + g * tau);
                            A(iq,j) = h + s * (g - h * tau);
                        }
                        for (size_type j = 0u; j < szc; ++j)
                        {
                            //indices are reversed in the following statements
                            //because V is a vector of column Vectors
                            g = V[ip][j]; //V(j,ip);
                            h = V[iq][j]; //V(j,iq);
                            V[ip][j] = g - s * (h + g * tau); //V(j,ip) = g - s * (h + g * tau);
                            V[iq][j] = h + s * (g - h * tau); //V(j,iq) = h + s * (g - h * tau);
                        }
                        ++nrot;
                    }
                
            }
        }
        for(size_type ip = 0u; ip < szc; ++ip)
        {
            b[ip] += z[ip];
            d[ip] = b[ip]; //Update d with the sum of tapq,
            z[ip] = 0.0; //and reinitialize z.
        }
    }
    error("Too many iterations", CURRENT_FUNCTION);
    return d;
}


void Matrix::LUdecomposition() const
//  This routine computes the LU decomposition of a rowwise permutation of a matrix. indx[1..n] is a state vector that records the row
//  permutation effected by the partial pivoting; d is stored as ±1 depending on whether the number of row interchanges was even or odd,
//  respectively. This routine is used in combination with linearSolve to solve linear equations or invert a matrix.
{
    const double TINY = 1.0e-25;
    verify(szc == szr);
    Vector vv(szc);                         //vv stores the implicit scaling of each row.
    determinant = 1.0;                      //No row interchanges yet.
    isLUdecomp = true;
    indx = new size_type[szc];
    LU = new double* [szc];
    double const *p  = data();
    for(size_type i = 0u; i < szc; ++i) {   //Loop over rows to get the implicit scaling information.
        LU[i] = new double[szc];
        double temp, big = 0.0;
        for(size_type j = 0u; j < szc; ++j, ++p) {
            LU[i][j] = *p;
            if((temp = fabs(LU[i][j])) > big) big = temp;
        }
        if(big == 0.0) error("Singular matrix", CURRENT_FUNCTION); //No nonzero largest element.
        vv[i] = 1.0 / big;                  //Save the scaling
    }
    
    for(size_type j = 0u; j < szc; ++j) {   //This is the loop over columns of Crout’s method.
        for(size_type i = 0u; i < j; ++i) { //This is equation (2.3.12) in Numerical Recipes except for i = j.
            double sum = LU[i][j];
            for(size_type k = 0u; k < i; ++k) sum -= LU[i][k] * LU[k][j];
            LU[i][j] = sum;
        }
        double big = 0.0, dum;              //Initialize for the search of the largest pivot element
        size_type imax;
        for(size_type i = j; i < szc; ++i) {//This is i >= j of equations (2.3.12) and (2.3.13) in Numerical Recipes
            double sum = LU[i][j];
            for (size_type k = 0u; k < j; ++k) sum -= LU[i][k] * LU[k][j];
            LU[i][j] = sum;
            if((dum = vv[i] * fabs(sum)) >= big) {  //Is the figure of merit for the pivot better than the best so far?
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {                            //Do we need to interchange rows?
            for(size_type k = 0u; k < szc; ++k) {   //Yes, do so...
                dum = LU[imax][k];
                LU[imax][k] = LU[j][k];
                LU[j][k] = dum;
            }
            determinant *= -1.0;                    //...and change the parity of d.
            vv[imax]=vv[j];                         //Also interchange the scale factor.
        }
        indx[j] = imax;
        if (LU[j][j] == 0.0) LU[j][j] = TINY;
        //If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
        //For some applications on singular matrices, it is desirable to substitute TINY for zero.
        if (j != szc - 1) {                         //Now, finally, divide by the pivot element.
            double dum = 1.0/(LU[j][j]);
            for (size_type i = j + 1u; i < szc; ++i) LU[i][j] *= dum;
        }
    }                                               //Go back for the next column in the reduction.
}

Vector Matrix::linearSolve(Vector b) const
//  Solves the set of n linear equations A·x = b, using the LU decomposition of the matrix A.
//  This routine takes into account the possibility that b will begin with many zero elements,
//  so it is efficient for use in matrix inversion.
{
    verify(szc == b.size() && szc == szr);
    if(!isLUdecomp) LUdecomposition();
    bool nonZero = false;
    for (size_type i = 0u, ii; i < szc; ++i) {  //When ii is set to a positive value, it will become the
        size_type ip = indx[i];                 //index of the first nonvanishing element of b. We now do
        double sum = b[ip];                     //the forward substitution, equation (2.3.6). The only new
        b[ip] = b[i];                           //wrinkle is to unscramble the permutation as we go.
        if (nonZero)
            for (size_type j = ii; j < i;j++) sum -= LU[i][j] * b[j];
        else if (sum) {                         //A nonzero element was encountered, so from now on we will
            ii = i;                             //have to do the sums in the loop above.
            nonZero = true;
        }
        b[i] = sum;
    }
    for(size_type i = szc; ;) { //Now we do the backsubstitution, equation (2.3.7).
        --i;
        double sum = b[i];
        for(size_type j = i + 1u; j < szc; ++j) sum -= LU[i][j] * b[j];
        b[i] = sum / LU[i][i];  //Store a component of the solution vector x.
        if(i == 0u) break;
    }                           
    return b;
}

Matrix Matrix::inverse() const
{
    Matrix result(szc, szr);
    for(size_type j = 0u; j < szr; ++j) {
        Vector col(szc, 0.0);
        col[j] = 1.0;
        result.setCol(j, linearSolve(col));
    }
    return result;
}

double Matrix::det() const
{
    if(!isLUdecomp) {
        LUdecomposition();
        for(size_type i = 0u; i < szc; ++i) determinant *= LU[i][i];
    }
    return determinant;
}
