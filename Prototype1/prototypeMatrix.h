#ifndef prototypeMatrix_H
#define prototypeMatrix_H

#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>

#define _MATH_DEFINES
#include <math.h>

using namespace Eigen;

struct superMatrix{

    int n;
    double delta_p;
    double delta_c;
    double k_p;
    double k_c;
    double v;
    double Omega_01;
    double Omega_12;

    superMatrix(int dimensions, double dP, double dC, double kP, double kC, double V, double Omega_p, double Omega_c);
    SparseMatrix<dcomplex> M(double delta_c, double v, double Gamma_10, double Gamma_21);
};

#endif