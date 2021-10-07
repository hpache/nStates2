#ifndef bVector_H
#define bVector_H
#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>
#define _Math_Defines
#include <math.h>
using namespace Eigen;



struct finalVector{
   int n;
   SparseVector<dcomplex> b;
   double Omega_01;



   finalVector(int dimensions, double Omega_01);
   SparseVector<dcomplex> getVector();
};
 
#endif