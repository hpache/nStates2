#include <Eigen/Sparse>
#include <iostream>
#include <complex>
#include <fstream>

#define _MATH_DEFINES
#include <math.h>

#include "prototypeMatrix.h"
#include "prototypeVector.h"

using namespace Eigen;

int main(){

    // Creating testVector object 
    Vector testVector(8, 0.5);
    // Testing if field is correct
    std::cout << "Vector Dimension: " << testVector.n << std::endl;

    // Creating testMatrix object
    superMatrix testMatrix(8, 0, 0, 1.0, 1.2, 0, 0.5, 0.5);
    // Testing if field is correct
    std::cout << "Matrix Dimension: " << testMatrix.n << std::endl;

}