#ifndef superMatrix_H
#define superMatrix_H
#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>
#define _Math_Defines
#include <math.h>
using namespace Eigen;



struct finalMatrix{


   int n;
   double delta_p;
   double delta_c;
   double k_p;
   double k_c;
   double v;
   double Omega_01;
   double Omega_12;


   finalMatrix(int n, double delta_p, double delta_c, double k_p, double k_c, double v, double Omega_01, double Omega_12);
   SparseMatrix<dcomplex> M(double delta_p, double delta_c, double k_p, double k_c, double v, double Omega_01, double Omega_12);
};
 
#endif