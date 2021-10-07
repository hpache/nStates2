#ifndef prototypeVector_H
#define prototypeVector_H


#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;

struct Vector{

    int n;
    double Omega_01;
    SparseVector<dcomplex> b;

    Vector(int dimensions, double Omega_p);
    SparseVector<dcomplex> getVector();
};

#endif