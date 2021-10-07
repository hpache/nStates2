#include "superMatrix.h"
#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>
#define _Math_Defines
#include <math.h>
using namespace Eigen;



finalMatrix::finalMatrix(int n, double delta_p, double delta_c, double k_p, double k_c, double v, double Omega_01, double Omega_12){
   n = n;
   delta_p = delta_p;
   delta_c = delta_c;
   k_p = k_p;
   k_c = k_c;
   v = v;
   Omega_01 = Omega_01;
   Omega_12 = Omega_12;
   SparseMatrix<dcomplex> A(n,n);
}

SparseMatrix<dcomplex> finalMatrix::M(double delta_p, double delta_c, double k_p, double k_c, double v, double Omega_01, double Omega_12){
   A.coeffRef(0,0) = dcomplex(-1.0*Gamma_10,0);
   A.coeffRef(0,1) = dcomplex(Gamma_21,0);
   A.coeffRef(0,2) = dcomplex(0,-1.0*Omega_01);
   A.coeffRef(0,4) = dcomplex(0,1.0*Omega_01);
   A.coeffRef(0,5) = dcomplex(0,1.0*Omega_12);
   A.coeffRef(0,7) = dcomplex(0,-1.0*Omega_12);
   A.coeffRef(1,1) = dcomplex(-1.0*Gamma_21,0);
   A.coeffRef(1,5) = dcomplex(0,-1.0*Omega_12);
   A.coeffRef(1,7) = dcomplex(0,1.0*Omega_12);
   A.coeffRef(2,0) = dcomplex(0,-2.0*Omega_01);
   A.coeffRef(2,1) = dcomplex(0,-1.0*Omega_01);
   A.coeffRef(2,2) = dcomplex(-0.5*Gamma_10,1.0*(delta_p+k_p*v));
   A.coeffRef(2,3) = dcomplex(0,1.0*Omega_12);
   A.coeffRef(3,2) = dcomplex(0,1.0*Omega_12);
   A.coeffRef(3,3) = dcomplex(-0.5*Gamma_21,1.0*(delta_c+delta_p-k_c*v+k_p*v));
   A.coeffRef(3,5) = dcomplex(0,-1.0*Omega_01);
   A.coeffRef(4,0) = dcomplex(0,2.0*Omega_01);
   A.coeffRef(4,1) = dcomplex(0,1.0*Omega_01);
   A.coeffRef(4,4) = dcomplex(-0.5*Gamma_10,-1.0*(delta_p+k_p*v));
   A.coeffRef(4,6) = dcomplex(0,-1.0*Omega_12);
   A.coeffRef(5,0) = dcomplex(0,1.0*Omega_12);
   A.coeffRef(5,1) = dcomplex(0,-1.0*Omega_12);
   A.coeffRef(5,3) = dcomplex(0,-1.0*Omega_01);
   A.coeffRef(5,5) = dcomplex(-0.5*Gamma_10-0.5*Gamma_21,1.0*(delta_c+delta_p-k_c*v+k_p*v-(delta_p+k_p*v)));
   A.coeffRef(6,4) = dcomplex(0,-1.0*Omega_12);
   A.coeffRef(6,6) = dcomplex(-0.5*Gamma_21,-1.0*(delta_c+delta_p-k_c*v+k_p*v));
   A.coeffRef(6,7) = dcomplex(0,1.0*Omega_01);
   A.coeffRef(7,0) = dcomplex(0,-1.0*Omega_12);
   A.coeffRef(7,1) = dcomplex(0,1.0*Omega_12);
   A.coeffRef(7,6) = dcomplex(0,1.0*Omega_01);
   A.coeffRef(7,7) = dcomplex(-0.5*Gamma_10-0.5*Gamma_21,1.0*(delta_p+k_p*v-(delta_c+delta_p-k_c*v+k_p*v)));
   return A;
}