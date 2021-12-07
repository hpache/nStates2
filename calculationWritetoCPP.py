# Henry Pacheco Cachon
# Created 12 October 2021
# File writes the calculation part of the cpp program
import sys
import subprocess
from sympy import *

# First run writetocpp
cmd = "writetocpp.py"
subprocess.call(["python3", cmd])

# Now we Initialize some important constants

# Initializing symbols related to probe and coupling lasers
dp0, dc0, kP, kC, v = symbols("delta_p delta_c k_p k_c v")

# Initializing symbols related to rabi frequencies
rabi01, rabi12 = symbols("Omega_01 Omega_12")


# Getting allowed transitions from arguments
transitions = [trans.split(",") for trans in sys.argv[1:]]

# Initializing constants
constants = [dp0, dc0, kP, kC, v, rabi01, rabi12] \
             + [Symbol("Gamma_{}{}".format(pair[0],pair[1])) for pair in transitions]

# Creating header lines
headerLines = ["#include <Eigen/Sparse>\n", "#include <iostream>\n", "#include <complex>\n", 
               "#include <fstream>\n", "#define _Math_Defines\n", 
               "#include <math.h>\n", "using namespace Eigen;\n"]

headerLines = ["#include \"superMatrix.h\"\n", "#include \"bVector.h\"\n"] + headerLines


# Lines containing all our calculation constants
calcConstants = ["\n\n\n\n"]
calcConstants = calcConstants + ["double kB = 1.381e-23;\n", "double AMU = 1.66e-27;\n", "double mAtom = 39;\n",
"double temperature = 323;\n", "double cellLength = 7.0; // This is in cm\n", 
"double u = sqrt(2 * kB * temperature / (mAtom * AMU)) * 100;\n", "double logP = 9.967 - 4646 / temperature;\n",
"double pressure = pow(10, logP);\n", "double density = (pressure) / (kB * temperature) * 1e-06;\n\n\n"]

# Constants directly related to our matrix and vectors and that don't depend on other constants
newGlobalVar = "double {} = {};\n"
askFormat = "{}: "

for constant in constants:

    # k_p and k_c are calculated from lamba
    if constant == Symbol("k_p") or constant == Symbol("k_c"):
        continue

    # v is the varying variable, don't need to give it a value
    if constant == Symbol("v"):
        calcConstants.append("double {};\n".format(constant))
        continue

    prompt = askFormat.format(constant)
    valInput = input(prompt)
    calcConstants.append(newGlobalVar.format(constant, valInput))

# Inputting wavelengths
lambdaPInput = input("lambda_p: ")
lambdaCInput = input("lambda_c: ")
lambdasFormat = "double {}: {} * 1e-07;\n"
calcConstants.append(lambdasFormat.format("lambdaP",lambdaPInput))
calcConstants.append(lambdasFormat.format("lambdaC",lambdaCInput))

# Inputting wavenumbers
kFormat = "double k_{} : 1 / {} * 1e-06;\n"
calcConstants.append(kFormat.format("p", "lambdaP"))
calcConstants.append(kFormat.format("c", "lambdaC"))

# Inputting weights
alphaFormat = "alpha{}_0 = (3 / (2 * M_PI)) * pow(lambda{}, 2) * density;\n"
calcConstants.append(alphaFormat.format("P","P"))
calcConstants.append(alphaFormat.format("C","C"))

# Inputting dimension
dimensionIn = input("number of levels: ")
calcConstants.append("double n = ({} * {}) - 1".format(dimensionIn, dimensionIn))

# Inputting detunings and ranges
initialDetuning = input("initialDetuning: ")
finalDetuning = input("finalDetuning: ")
numDetuning = input("numberDetunings: ")
initialVel = input("initialVelocity: ")
finalVel = input("finalVelocity: ")
numVel = input("numberVelocity: ")
calcConstants.append(newGlobalVar.format("initialDetuning", initialDetuning))
calcConstants.append(newGlobalVar.format("finalDetuning", finalDetuning))
calcConstants.append(newGlobalVar.format("numberDetunigns", numDetuning))
calcConstants.append(newGlobalVar.format("initialVelocity", initialVel))
calcConstants.append(newGlobalVar.format("finalVelocity", finalVel))
calcConstants.append(newGlobalVar.format("numberVelocity", numVel))

# Appending filename
filename = input("file name: ")
calcConstants.append("std::string filename = \"{}\";\n".format(filename))


# Writing eigen lines
eigenLines = ["\n\n\n"] + ["SparseLU<SparseMatrix<dcomplex>, COLAMDOrdering<int> > solver;\n",
"SparseMatrix<dcomplex> B;\n", "SparseVector<dcomplex> x(15);\n"]


# Write main method
main = ["\n\n\n", "int main() {\n"]
# Initialize calculation params
main = main + ["   double dv;\n", "   double dc;\n", "   double vAtom;\n", "   double coupling;\n",
               "   double weight;\n"]
# Initialize file writing
main = main + ["   std::ofstream myFile;\n", "   myFile.open(fileName);\n", 
               "   myFile << \"Coupling Detunings,TransmissionP_D \\n\";\n"]

# Start calculation body 1
main = main + ["   SparseVector<dcomplex> alphaP_D(numberDetunings);\n", 
               "   dv = (finalVelocity - initialVelocity) / (numberVelocity - 1);\n",
               "   dc = (finalDetuning - initialDetuning) / (numberDetunings - 1);\n"]
main = main + ["   finalVector b(n, Omega_01);\n"]

# Start calculation body 2
main = main + ["   for (int i = 0; i < numberDetunings; i++){\n", "      coupling = initialDetuning + i * dc;\n",
               "      for (int k = 0; k < numberVelocity; k++){\n", 
               "         v = initialVelocity + k * dv;\n", 
               "         weight = 1 / (sqrt(M_PI) * u) * exp(-1.0 * (pow(vAtom, 2) / pow(u, 2))) * dv;\n",
               "         finalMatrix BMatrix(int n, double delta_p, double delta_c, double k_p, double k_c, double v, double Omega_01, double Omega_12, double Gamma_21, double Gamma_10);\n",
               "         B = BMatrix.M(double delta_p, double delta_c, double k_p, double k_c, double v, double Omega_01, double Omega_12, double Gamma_21, double Gamma_10);\n",
               "         solver.analyzePattern(B);\n",
               "         solver.factorize(B);\n",
               "         x = solver.solve(b.getVector());\n",
               "         alphaP_D.coeffRef(i) += weight * alphaP_0 * (Gamma_10 / Omega_01) * ((dcomplex(0, 1) * x.coeffRef(5)).real());\n",
               "      }\n",
               "      myFile << coupling << "," << exp(-cellLength * alphaP_D.coeffRef(i)).real() << \"\\n\";\n",
               "   }\n"]

# Start ending
main = main + ["   std::cout << \"Finished\" << std::endl;\n"]
main = main + ["}\n"]

fullLines = headerLines + calcConstants + main

mainCPP = open("main.cpp", "w")
mainCPP.writelines(fullLines)
mainCPP.close()