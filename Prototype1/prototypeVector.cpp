/**
 * Henry Pacheco Cachon
 * Created: 22 September 2021
 * This file is a prototype vector file. This will work as a basic format for vector struct 
 * creation from the python program
 * **/

#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

#include "prototypeVector.h"

using namespace Eigen;


Vector::Vector(int dimensions, double Omega_01){
    n = dimensions;
    Omega_01 = Omega_01;
    SparseVector<dcomplex> b(n);
}

SparseVector<dcomplex> Vector::getVector(){

    b.coeffRef(2) = dcomplex(0, -1.0 * Omega_01);
    b.coeffRef(4) = dcomplex(0, 1.0 * Omega_01);

    return b;
}

