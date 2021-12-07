# Henry Pacheco Cachon
# Created 26 October 2021

import nStates as ns
from sympy import *
from sympy.printing.c import C11CodePrinter
from re import search
import numpy as np

# Class represents a code wrapper object
class codeWrapper:

    headerlines = ["#include <Eigen/Sparse>\n", "#include <iostream>\n", "#include <complex>\n", 
                   "#include <fstream>\n", "#define _Math_Defines\n", 
                   "#include <math.h>\n", "using namespace Eigen;\n"]
    
    newGlobalVar = "double {}"

    closingLines = ["};\n \n", "#endif"]

    matrixConstructor = "\n\n   finalMatrix(int n, "
    matrixFunction = "   SparseMatrix<dcomplex> M(int n, "

    vectorConstruct = "\n\n\n   finalVector(int dimensions, "
    vectorFunction = "   SparseVector<dcomplex> getVector(int n, "

    ''' Constructor method, input is the number of states, the hamiltonian, allowed transitions, constants, and condition'''
    def __init__(self, states: list, transitions: list) -> None:



        # Finding the total number of states
        self.numStates = sum(states)
        # Initializing list containing all the transitions
        self.transitions = transitions
        # Finding the dimension for our final arrays
        self.dimension = ( self.numStates ** 2 ) - 1

        # Initializing empty list containing our constants
        self.constants = [Symbol("Gamma_{}{}".format(pair[0],pair[1])) for pair in transitions]

        # Initializing list holding our states
        self.states = states

        # Initializing a variable to hold the hamiltonian
        self.Hamiltonian = None

        # Initializing variable holding our rho condition
        rho00 = 1
        for i in range(1, self.numStates):
            rho00 -= Symbol('rho_{}{}'.format(i,i))
        self.condition = {Symbol('rho_00') : rho00}

    ''' Method returns the list of constants '''
    def getMatrixConstants(self) -> list:
        
        output = []

        for constant in self.constants:

            output.append("{}".format(constant))
        
        return output


    ''' Method returns the list of constants unique to the vector '''
    def getVectorConstants(self) -> list:
        return self.vectorConstants

    
    ''' Method returns the dimension of the system '''
    def getDimension(self) -> int:
        return self.dimension
    
    
    ''' Auxiliary method, takes in a symbolic matrix and converts it to an eigen formatted matrix '''
    def _parseArray(self, symbolicArray, vector = False) -> list:

        # Preparing C code printer 
        printer = C11CodePrinter()

        if vector == False:
            A = MatrixSymbol('A',self.dimension, self.dimension)
        else:
            A = MatrixSymbol('b', self.dimension,1)

        # Converting the symbolic matrix to a C array similar to the matrix symbol A
        # For some reason this is a 1-d array as opposed to a 2-d array 
        test = printer.doprint(symbolicArray,assign_to=A)

        # printer.doprint() outputs a long string with each matrix element separated by a \n
        # The Lines variable is a list which separates each matrix element 
        Lines = test.split('\n')

        # Start row counter
        row = 0

        # Creating formatted strings for matrix elements
        # These formats are for the Eigen library for C++
        matrixAssign = ".coeffRef({},{})"
        vectorAssign = ".coeffRef({})"
        complexFormat = "dcomplex({},{});"

        # List containing the newly formatted matrix elements. 
        newLines = []
        for string in Lines:

            #Getting Indexes i.e. test[1], the index will be between the brackets
            l = string.index('[')
            r = string.index(']')
            number = int(string[l+1:r])

            x = number % self.dimension 
            y = row

            #Get values i.e. test[1] = 0.5*Omega_23; The value is between the = and ;
            eq = string.index('=')
            stop = string.index(';')
            value = string[eq+2:stop]
            
            # Ignoring all the zero values in order to make a sparse matrix
            if value != '0':

                # Finding the complex values within the matrix 
                if search('I',value):

                    # Looking for values like test[1] = 0.5*I*Gamma_10
                    if value.rfind(' ') == -1:
                        pro = value.index('*')
                        newComplex = complexFormat.format(0,value.replace('I*',''))
                    
                    # Looking for values like test[1] = 0.5*Gamma10 + I(0.5*Omega_01 + Omega_22)
                    else:
                        bigBoi = value.split(' ')
                        a = ''
                        b = ''
                        check = 0

                        for i in range(len(bigBoi)):
                            val = bigBoi[i]
                            if val.find("I") != -1:
                                check = i
                        
                        bigBoi[check:] = [''.join(bigBoi[check:])]

                        for i in range(len(bigBoi)):
                            val = bigBoi[i]
                            if val.find("I") != -1:
                                sign = bigBoi[i-1]
                                if sign == "+":
                                    b += val
                                else:
                                    b = sign + val
                            else:
                                a += val
                        
                        a = a[:-1]
                        b = b.replace("I*",'')
                        newComplex = complexFormat.format(a,b)

                # Looking for real numbers   
                else:
                    newComplex = complexFormat.format(value,0)

                # Creating final format of matrix element in Eigen syntax
                if vector == False:
                    newString = "   " + string[:l] + matrixAssign.format(y,x) + string[r+1:r+4] + newComplex +'\n'
                else:
                    newString = "   " + string[:l] + vectorAssign.format(x) + string[r+1:r+4] + newComplex +'\n'
                newLines.append(newString)

            # Moving on to the next row
            if x == self.dimension - 1:
                row +=1

        return newLines
        
    
    ''' Method gets all unique variables in the b vector '''
    def _distinctVariables(self,b) -> set:

        uniqueVars = set()
        stringConstants = ["{}".format(constant) for constant in self.constants]

        for i in range(len(b)):
            for j in range(len(stringConstants)):

                if search(stringConstants[j], b[i]):
                    uniqueVars.add(stringConstants[j])

        return uniqueVars


    ''' Aux method, creates omega matrix '''
    def _createOmegaMatrix(self) -> list:

        numberGround = self.states[0]
        numberInter = self.states[1]
        numberExcited = self.states[2]

        newConstants = set()

        omega = [[0 for i in range(self.numStates)] for j in range(self.numStates)]

        for i in range(numberGround):
            for j in range(numberInter):

                ground_inter = Symbol("Omega_{}{}".format(i,j + numberGround))
                omega[i][j + numberGround] = ground_inter
                omega[j + numberInter][i] = ground_inter
                newConstants.add(ground_inter)

        for i in range(numberInter):
            for j in range(numberExcited):
                inter_exc = Symbol("Omega_{}{}".format(i + numberGround, j + numberGround + numberInter))
                omega[i + numberGround][j + numberGround + numberInter] = inter_exc
                omega[j + numberGround + numberInter][i + numberGround] = inter_exc
                newConstants.add(inter_exc)
        
        for constant in newConstants:
            self.constants.append(constant)

        return omega


    ''' Aux method, creates free-atom hamiltonian '''
    def _createFreeAtomMatrix(self, energyLevels: list) -> list:

        freeAtom  = [[0 for i in range(self.numStates)] for j in range(self.numStates)]

        for i in range(self.numStates):
            for j in range(self.numStates):

                if i == j:
                    
                    freeAtom[i][j] = energyLevels[i]
        
        return freeAtom


    ''' Method allows user to add constants '''
    def addToConstantsList(self, newItem: Symbol) -> None:
        self.constants.append(newItem)
    

    ''' Method returns the hamiltonian '''
    def getHamiltonian(self) -> list:
        return self.Hamiltonian
    
    
    ''' Method returns the list of constants '''
    def getConstants(self) -> list:
        return self.constants

    ''' Method creates the hamiltonian for our system '''
    def createHamiltonian(self, energyLevels: list) -> None:

        omegaMatrix = self._createOmegaMatrix()
        freeAtomMatrix = self._createFreeAtomMatrix(energyLevels)
        hamiltonian = [[i for i in range(self.numStates)] for j in range(self.numStates)]

        for i in range(self.numStates):
            for j in range(self.numStates):
                hamiltonian[i][j] = omegaMatrix[i][j] + freeAtomMatrix[i][j]
        
        self.Hamiltonian = np.array(hamiltonian)

    
    ''' Do physics '''
    def doPhysics(self) -> None:

        model = ns.States(numberStates=self.numStates, hamiltonian=self.Hamiltonian, staticDictionary=self.condition)

        model.lindbladMatrix(self.transitions)
        model.densityMatrix()
        model.superMatrix()
        model.makeEquations()

        self.symbolicMatrix = model.changeParam(self.condition, symbolic=True)
        self.b = model.getVector()


    ''' Method writes the h file for the matrix struct '''
    def writeMatrixH(self):

        # Creating includes and define lines as well as namespace lines
        hFileFirstLine = ["#ifndef superMatrix_H\n", "#define superMatrix_H\n"] + self.headerlines
                        
        # Starting the struct creation
        structName = ["\n\n\nstruct finalMatrix{\n\n\n"]

        # Writing body of the struct
        # Start with global variables which is all of the constants
        globalVars = ["   int n;\n"] + ["   " + self.newGlobalVar.format(i) + ";\n" for i in self.constants]


        for i in range(len(self.constants)):

            if i == len(self.constants) - 1:
                self.matrixConstructor += self.newGlobalVar.format(self.constants[i]) + ");\n"
                self.matrixFunction += self.newGlobalVar.format(self.constants[i]) + ");\n"
            else:
                self.matrixConstructor += self.newGlobalVar.format(self.constants[i]) + ", "
                self.matrixFunction += self.newGlobalVar.format(self.constants[i]) + ", "

        # We add everything together to one list and write it to an h file
        fullHLines = hFileFirstLine + structName + globalVars + [self.matrixConstructor,self.matrixFunction] + self.closingLines

        matrixHfile = open("superMatrix.h", "w")
        matrixHfile.writelines(fullHLines)
        matrixHfile.close()


    ''' Method writes the h file for the vector struct '''
    def writeVectorH(self):

        # Get vector in Eigen format 
        self.vectorEigen = self._parseArray(self.b,vector = True)
        # Get distinct constants from vector
        self.vectorConstants = list(self._distinctVariables(self.vectorEigen))


        # Creating header lines
        vectorHFirstLines = ["#ifndef bVector_H\n", "#define bVector_H\n"] + self.headerlines

        # Create struct init
        vectorStruct = ["\n\n\nstruct finalVector{\n"]

        # Creating fields
        fields = ["   int n;\n",] + ["   " + self.newGlobalVar.format(c) + ";\n" for c in self.vectorConstants]


        # Adding to constructor method
        for i in range(len(self.vectorConstants)):

            if i == len(self.vectorConstants) - 1:
                self.vectorConstruct += self.newGlobalVar.format(self.vectorConstants[i]) + ");\n"
                self.vectorFunction += self.newGlobalVar.format(self.vectorConstants[i]) + ");\n"
            else:
                self.vectorConstruct += self.newGlobalVar.format(self.vectorConstants[i]) + ", "
                self.vectorFunction += self.newGlobalVar.format(self.vectorConstants[i]) + ", "
                

        # Adding all code lines together
        vectorFullLines = vectorHFirstLines + vectorStruct + fields + [self.vectorConstruct, self.vectorFunction] + self.closingLines

        # Writing h file
        vectorHFile = open("bVector.h","w")
        vectorHFile.writelines(vectorFullLines)
        vectorHFile.close()
    

    ''' Method writes the matrix cpp file '''
    def writeVectorCPP(self):

        # Write full constructor method
        vectorFullConstructor = ["\n\n\nfinalVector::" + self.vectorConstruct[6:-2] + "{\n", "   n = dimensions;\n"]
        # Write full getVector method
        vectorFullFunction = ["SparseVector<dcomplex> finalVector::" + self.vectorFunction[26:-2] + "{\n", "   SparseVector<dcomplex> b(n);\n"] + self.vectorEigen + ["   return b;\n","}"]
        

        # Add all fields and initialize them
        for i in range(len(self.vectorConstants)):
            vectorFullConstructor.append("   " + self.vectorConstants[i] + " = " + self.vectorConstants[i] + ";\n")
        vectorFullConstructor.append("}\n\n")

        # Add all the lines together
        vectorCPPFullLines = ["#include \"bVector.h\"\n"] + self.headerlines + vectorFullConstructor + vectorFullFunction

        # Write to cpp file
        vectorCPP = open("bVector.cpp","w")
        vectorCPP.writelines(vectorCPPFullLines)
        vectorCPP.close()
    

    ''' Method writes the vector cpp file '''
    def writeMatrixCPP(self):

        # Getting matrix in the eigenformat
        eigenMatrix = self._parseArray(self.symbolicMatrix, vector = False)

        # Writing full constructor method
        matrixFullConstructor = ["\n\n\nfinalMatrix::" + self.matrixConstructor[5:-2] + "{\n", "   n = n;\n"]
        # Writing full function method
        matrixFullFunction = ["SparseMatrix<dcomplex> finalMatrix::" + self.matrixFunction[26:-2] + "{\n"] + ["   SparseMatrix<dcomplex> A(n,n);\n"]
        matrixFullFunction = matrixFullFunction + eigenMatrix + ["   A.makeCompressed();\n","   return A;\n", "}"]

        for i in range(len(self.constants)):
            matrixConstantString = "{}".format(self.constants[i])
            matrixFullConstructor.append("   " + matrixConstantString + " = " + matrixConstantString + ";\n")
        matrixFullConstructor.append("}\n\n")

        matrixCPPFullLines = ["#include \"superMatrix.h\"\n"] + self.headerlines + matrixFullConstructor + matrixFullFunction

        matrixCPP = open("superMatrix.cpp", "w")
        matrixCPP.writelines(matrixCPPFullLines)
        matrixCPP.close()
        
        





if __name__ == "__main__":
    
    test = codeWrapper([1,1,1],[(1,0),(2,1)])
    test.addToConstantsList(Symbol("delta_p"))
    test.addToConstantsList(Symbol("delta_c"))
    test.createHamiltonian([0,Symbol("delta_p"),Symbol("delta_c")])
    test.doPhysics()
    test.writeMatrixH()
    test.writeVectorH()
    test.writeVectorCPP()
    test.writeMatrixCPP()
    print(test.getMatrixConstants())
    print(test.getVectorConstants())

    


    

