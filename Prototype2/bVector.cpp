#include "bVector.h"
#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>
#define _Math_Defines
#include <math.h>
using namespace Eigen;



finalVector::finalVector(int dimensions, double Omega_01){
   n = dimensions;
   Omega_01 = Omega_01;
   SparseVector<dcomplex> b(n);
}

SparseVector<dcomplex> finalVector::getVector(){
   b.coeffRef(2) = dcomplex(0,-1.0*Omega_01);
   b.coeffRef(4) = dcomplex(0,1.0*Omega_01);
   return b;
}