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

#ifndef linalg_h
#define linalg_h
#include <vector>
#include <iostream>
#include <cmath>

class Matrix;      // all kinds of matrix operations 

class Vector: public std::vector<double>
{
public:
    Vector(const size_type &size = 1u, const double &val = 0.0) : std::vector<double>(size, val) {}
    Vector(const std::vector<double> &vec) : std::vector<double>(vec) {}
    friend std::ostream& operator<< (std::ostream&, const Vector&);
    Vector& operator=(const Vector&);
    Vector& operator+=(const Vector&);
    Vector& operator-=(const Vector&);
    double  operator*(const Vector&) const;
    Vector& operator*=(const double&);
    
    Vector operator+ (const Vector &other) const
    {
        Vector result(*this);
        return result += other;
    }
    Vector operator- (const Vector &other) const
    {
        Vector result(*this);
        return result -= other;
    }
    Vector operator*(const double &factor) const
    {
        Vector result(*this);
        return result *= factor;
    }
    friend Vector operator*(const double &factor, const Vector &x) {return x * factor;}
    Vector& operator/=(const double &factor) {return (*this) *= (1.0 / factor);}
    Vector operator/(const double &factor) const
    {
        Vector result(*this);
        return result *= (1.0 / factor);
    }
    double norm() const {return sqrt(this->operator*(*this));}
    double total() const;
    void normalize() {this->operator/=(this->norm());}
    friend double distance(const Vector &x, const Vector &y)
    {
        Vector tmp(x - y);
        return tmp.norm();
    }
private:
};

class Matrix: public std::vector<double>
{
public:
    Matrix(const size_type &sizer = 1u, const size_type &sizec = 1u, const double &val = 0.0) :
    szr(sizer), szc(sizec), std::vector<double>(sizer * sizec, val), isLUdecomp(false) {}
    Matrix(const std::vector<std::vector<double> >&);
    Matrix(const std::vector<Vector>&);
    Matrix(const Vector&);
    Matrix(const Matrix &mat) : szr(mat.szr), szc(mat.szc), std::vector<double>(mat), isLUdecomp(false) {}
    ~Matrix();
    friend std::ostream& operator<< (std::ostream&, const Matrix&);
    Matrix& operator=(const Matrix&);
    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);
    Matrix& operator*=(const double&);
    Vector operator*(const Vector&) const;
    friend Vector operator*(const Vector&, const Matrix&);
    Matrix operator*(const Matrix&) const;
    Matrix operator+(const Matrix &other) const
    {
        Matrix result(*this);
        return result += other;
    }
    Matrix operator-(const Matrix &other) const
    {
        Matrix result(*this);
        return result -= other;
    }
    Matrix operator*(const double &factor) const
    {
        Matrix result(*this);
        return result *= factor;
    }
    friend Matrix operator*(const double &factor, const Matrix &x) {return x * factor;}
    Matrix& operator/=(const double &factor) {return (*this) *= (1.0 / factor);}
    Matrix operator/(const double &factor) const
    {
        Matrix result(*this);
        return result *= (1.0 / factor);
    }
    Matrix transpose() const;
    Matrix symPart() const {return 0.5 * (*this + this->transpose());}
    Matrix aSymPart() const {return 0.5 * (*this - this->transpose());}
    Vector eigenvaluesSymPart(std::vector<Vector>&) const;
    Vector linearSolve(Vector) const;
    Matrix inverse() const;
    double trace() const {return getDiag().total();}
    double det() const;
    double& operator()(const size_type &i, const size_type &j) {return *(data() + i * szc + j);}
    double operator()(const size_type &i, const size_type &j) const {return *(data() + i * szc + j);}
    void setRow(const size_type&, const Vector&);
    void setCol(const size_type&, const Vector&);
    void setDiag(const Vector&);
    Vector getRow(const size_type&) const;
    Vector getCol(const size_type&) const;
    std::vector<Vector> getRows() const;
    std::vector<Vector> getCols() const;
    Vector getDiag() const;
    size_type sizec() const {return szc;}
    size_type sizer() const {return szr;}
private:
    void LUdecomposition() const;
    size_type szr, szc;
    mutable size_type *indx;
    mutable bool isLUdecomp;
    mutable double determinant, **LU;
};


//End of #ifdef linalg_h
#endif

