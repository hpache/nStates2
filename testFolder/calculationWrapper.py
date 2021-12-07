# This class is meant to write the calculation cpp file

from sympy import *
import codeWrapper as cw

class calculationWrapper:

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


    ''' Constructor method '''
    def __init__(self, constants: list, dimension: int) -> None:
        self.constants = constants
        self.matrixParameters = None
        self.vectorParameters = None
        self.dimension = dimension
        self.location = None
    

    ''' Method sets the matrix function string '''
    def setMatrixParameters(self, constants: list) -> None:
        self.matrixParameters = ""

        for i in range(len(constants)):

            if i == len(constants) - 1:
                self.matrixParameters += constants[i]
            else:
                self.matrixParameters += constants[i] + ", "
        
    

    ''' Method sets the vector function string '''
    def setVectorParameters(self, constants: list) -> None:
        self.vectorParameters= ""

        for i in range(len(constants)):

            if i == len(constants) - 1:
                self.vectorParameters += constants[i]
            else:
                self.vectorParameters += constants[i] + ", "
    

    ''' Method sets the location for final solution '''
    def setLocation(self, index) -> None:
        self.location = index


    ''' This method appends to the calculation constants by asking for them '''
    def askConstants(self) -> None:

        for constant in self.constants:

            # k_p and k_c are calculated from lamba
            if constant == Symbol("k_p") or constant == Symbol("k_c"):
                continue

            # v is the varying variable, don't need to give it a value
            if constant == Symbol("v"):
                self.calcConstants.append("double {};\n".format(constant))
                continue

            prompt = self.askFormat.format(constant)
            valInput = input(prompt)
            self.calcConstants.append(self.newGlobalVar.format(constant, valInput))


    ''' This method sets the experiment parameters '''
    def setParameters(self) -> None:

        # Inputting wavelengths
        lambdaPInput = input("lambda_p: ")
        lambdaCInput = input("lambda_c: ")
        lambdasFormat = "double {} = {} * 1e-07;\n"
        self.calcConstants.append(lambdasFormat.format("lambdaP",lambdaPInput))
        self.calcConstants.append(lambdasFormat.format("lambdaC",lambdaCInput))

        # Inputting wavenumbers
        kFormat = "double k_{} = 1 / {} * 1e-06;\n"
        self.calcConstants.append(kFormat.format("p", "lambdaP"))
        self.calcConstants.append(kFormat.format("c", "lambdaC"))

        # Inputting weights
        alphaFormat = "double alpha{}_0 = (3 / (2 * M_PI)) * pow(lambda{}, 2) * density;\n"
        self.calcConstants.append(alphaFormat.format("P","P"))
        self.calcConstants.append(alphaFormat.format("C","C"))

        # Inputting dimension
        dimensionIn = input("number of levels: ")
        self.calcConstants.append("double n = ({} * {}) - 1; \n".format(dimensionIn, dimensionIn))

        # Inputting detunings and ranges
        initialDetuning = input("initialDetuning: ")
        finalDetuning = input("finalDetuning: ")
        numDetuning = input("numberDetunings: ")
        initialVel = input("initialVelocity: ")
        finalVel = input("finalVelocity: ")
        numVel = input("numberVelocity: ")
        self.calcConstants.append(self.newGlobalVar.format("initialDetuning", initialDetuning))
        self.calcConstants.append(self.newGlobalVar.format("finalDetuning", finalDetuning))
        self.calcConstants.append(self.newGlobalVar.format("numberDetunings", numDetuning))
        self.calcConstants.append(self.newGlobalVar.format("initialVelocity", initialVel))
        self.calcConstants.append(self.newGlobalVar.format("finalVelocity", finalVel))
        self.calcConstants.append(self.newGlobalVar.format("numberVelocity", numVel))
    

    ''' Writing the body of the cpp file '''
    def writeBody(self) -> None:

        # Appending filename
        filename = input("file name: ")
        self.calcConstants.append("std::string fileName = \"{}\";\n".format(filename))


        # Writing eigen lines
        eigenLines = ["\n\n\n"] + ["SparseLU<SparseMatrix<dcomplex>, COLAMDOrdering<int> > solver;\n",
        "SparseMatrix<dcomplex> B;\n", "SparseVector<dcomplex> x({});\n".format(self.dimension), 
        "SparseVector<dcomplex> b({});\n".format(self.dimension)]
        self.calcConstants = self.calcConstants + eigenLines


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
        main = main + ["   finalVector bVector(n, {});\n".format(self.vectorParameters), 
                       "   SparseVector<dcomplex> b = bVector.getVector(n, {});\n".format(self.vectorParameters)]

        # Start calculation body 2
        main = main + ["   for (int i = 0; i < numberDetunings; i++){\n", "      delta_c = initialDetuning + i * dc;\n",
                    "      for (int k = 0; k < numberVelocity; k++){\n", 
                    "         v = initialVelocity + k * dv;\n", 
                    "         weight = 1 / (sqrt(M_PI) * u) * exp(-1.0 * (pow(vAtom, 2) / pow(u, 2))) * dv;\n",
                    "         finalMatrix BMatrix(n, {});\n".format(self.matrixParameters),
                    "         SparseMatrix<dcomplex> B = BMatrix.M(n, {});\n".format(self.matrixParameters),
                    "         solver.analyzePattern(B);\n",
                    "         solver.factorize(B);\n",
                    "         x = solver.solve(b);\n",
                    "         alphaP_D.coeffRef(i) += weight * alphaP_0 * (Gamma_10 / Omega_01) * ((dcomplex(0, 1) * x.coeffRef({})).real());\n".format(self.location),
                    "      }\n",
                    "      myFile << delta_c << \",\" << exp(-cellLength * alphaP_D.coeffRef(i)).real() << \"\\n\";\n",
                    "   }\n"]

        # Start ending
        main = main + ["   std::cout << \"Finished\" << std::endl;\n"]
        main = main + ["}\n"]

        fullLines = self.headerLines + self.calcConstants + main

        mainCPP = open("main.cpp", "w")
        mainCPP.writelines(fullLines)
        mainCPP.close()


if __name__ == "__main__":

    test = cw.codeWrapper([1,1,1],[(1,0),(2,1)])
    test.addToConstantsList(Symbol("delta_p"))
    test.addToConstantsList(Symbol("delta_c"))
    test.addToConstantsList(Symbol("k_p"))
    test.addToConstantsList(Symbol("k_c"))
    test.addToConstantsList(Symbol("v"))
    test.createHamiltonian([0,Symbol("delta_p") + Symbol("k_p") * Symbol("v"), Symbol("delta_c") -  Symbol("k_c") * Symbol("v")])
    test.doPhysics()
    test.writeMatrixH()
    test.writeVectorH()
    test.writeVectorCPP()
    test.writeMatrixCPP()
    

    test1 = calculationWrapper(test.getConstants(), test.getDimension())
    test1.setMatrixParameters(test.getMatrixConstants())
    test1.setVectorParameters(test.getVectorConstants())
    test1.setLocation(4)
    test1.askConstants()
    test1.setParameters()
    test1.writeBody()