# Henry Pacheco Cachon
# Created: 27 September 2021
# The purpose of this file is to write a cpp matrix and vector file and get them ready to 
# be used by a main.cpp file

from turtle import clear
from nStates import *
from sympy import *
import numpy as np
from sympy.printing.c import C11CodePrinter
from re import search


''' Defining parser function '''
def parseArray(dimension, symbolicArray, vector = False):

    # Preparing C code printer 
    printer = C11CodePrinter()

    if vector == False:
        A = MatrixSymbol('A',dimension,dimension)
    else:
        A = MatrixSymbol('b',dimension,1)

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

        x = number % dimension 
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
        if x == dimension - 1:
            row +=1

    return newLines

''' Defining function which looks for all the distinct variables in a vector'''
def distinctVariables(vector, constants):

    uniqueVars = set()
    stringConstants = ["{}".format(constant) for constant in constants]

    for i in range(len(vector)):
        for j in range(len(stringConstants)):

            if search(stringConstants[j], vector[i]):
                uniqueVars.add(stringConstants[j])

    return uniqueVars


''' Begin by initializing the nstates matrix '''

# Initializing symbols related to probe and coupling lasers
dp0, dc0, kP, kC, v = symbols("delta_p delta_c k_p k_c v")

# Initializing symbols related to rabi frequencies
rabi01, rabi12 = symbols("Omega_01 Omega_12")

# Initializing the hamiltonian of the system
H = np.array([
    [0, rabi01, 0],
    [rabi01, dp0 + (kP * v), rabi12 ],
    [0, rabi12, (dp0 + (kP * v)) + (dc0 - (kC *v))]
    ])

# Initializing allowed transitions
transitions = [(2,1),(1,0)]

# Initializing substituiton of rho00
static = {Symbol('rho_00') : 1 - Symbol('rho_11') - Symbol('rho_22')}

# Initializing constants
constants = [dp0, dc0, kP, kC, v, rabi01, rabi12] + [Symbol("Gamma_{}{}".format(pair[0],pair[1])) for pair in transitions]

# Initializing three level model
threeLevel = States(numberStates = 3, hamiltonian = H, staticDictionary = static)


# Doing physics
threeLevel.lindbladMatrix(transitions)
threeLevel.densityMatrix()
threeLevel.superMatrix()
threeLevel.makeEquations()

#Getting symbolic version of the 8 x 8 matrix
symbolicMatrix = threeLevel.changeParam(static, symbolic = true)

# Also have to get the vector b
b = threeLevel.getVector()

# Initializing the dimensions of the symbolic matrix
dimension = (3 ** 2) - 1

''' Getting common lines ready '''
headerLines = ["#include <Eigen/Sparse>\n", "#include <iostream>\n", "#include <complex>\n", 
               "#include <fstream>\n", "#define _Math_Defines\n", 
               "#include <math.h>\n", "using namespace Eigen;\n"]

''' Constructing the h file for the matrix struct '''

# Creating includes and define lines as well as namespace lines
hFileFirstLine = ["#ifndef superMatrix_H\n", "#define superMatrix_H\n"] + headerLines
                
# Starting the struct creation
structName = ["\n\n\nstruct finalMatrix{\n\n\n"]

# Writing body of the struct
# Start with global variables which is all of the constants
newGlobalVar = "double {}"
globalVars = ["   int n;\n"] + ["   " + newGlobalVar.format(i) + ";\n" for i in constants]

# Next we write the constructor method and matrix function
constructor = "\n\n   finalMatrix(int n, "
function = "   SparseMatrix<dcomplex> M("


for i in range(len(constants)):

    if i == len(constants) - 1:
        constructor += newGlobalVar.format(constants[i]) + ");\n"
        function += newGlobalVar.format(constants[i]) + ");\n"
    else:
        constructor += newGlobalVar.format(constants[i]) + ", "
        function += newGlobalVar.format(constants[i]) + ", "

# Finally we write the closing bracket and finish the h file
closingLines = ["};\n \n", "#endif"]

# We add everything together to one list and write it to an h file
fullHLines = hFileFirstLine + structName + globalVars + [constructor,function] + closingLines

matrixHfile = open("superMatrix.h", "w")
matrixHfile.writelines(fullHLines)
matrixHfile.close()


''' Create H file for vector'''

# Get vector in Eigen format 
vectorEigen = parseArray(dimension,b,vector = True)
# Get distinct constants from vector
vectorConstants = list(distinctVariables(vectorEigen,constants))


# Creating header lines
vectorHFirstLines = ["#ifndef bVector_H\n", "#define bVector_H\n"] + headerLines

# Create struct init
vectorStruct = ["\n\n\nstruct finalVector{\n"]

# Creating fields
fields = ["   int n;\n", "   SparseVector<dcomplex> b;\n"] + ["   " + newGlobalVar.format(c) + ";\n" for c in vectorConstants]

# Create constructor and function
vectorConstruct = "\n\n\n   finalVector(int dimensions, "
vectorFunction = "   SparseVector<dcomplex> getVector();\n"

# Adding to constructor method
for i in range(len(vectorConstants)):

    if i == len(vectorConstants) - 1:
        vectorConstruct += newGlobalVar.format(vectorConstants[i]) + ");\n"
    else:
        vectorConstruct += newGlobalVar.format(vectorConstants[i]) + ", "

# Adding all code lines together
vectorFullLines = vectorHFirstLines + vectorStruct + fields + [vectorConstruct, vectorFunction] + closingLines

# Writing h file
vectorHFile = open("bVector.h","w")
vectorHFile.writelines(vectorFullLines)
vectorHFile.close()

''' Write vector.cpp file '''

# All header files are set, just need to write out the actual methods and constructors

# Write full constructor method
vectorFullConstructor = ["\n\n\nfinalVector::" + vectorConstruct[6:-2] + "{\n", "   n = dimensions;\n"]
# Write full getVector method
vectorFullFunction = ["SparseVector<dcomplex> finalVector::" + vectorFunction[26:-2] + "{\n"] + vectorEigen + ["   return b;\n","}"]

# Add all fields and initialize them
for i in range(len(vectorConstants)):
    vectorFullConstructor.append("   " + vectorConstants[i] + " = " + vectorConstants[i] + ";\n")
vectorFullConstructor.append("   SparseVector<dcomplex> b(n);\n")
vectorFullConstructor.append("}\n\n")

# Add all the lines together
vectorCPPFullLines = ["#include \"bVector.h\"\n"] + headerLines + vectorFullConstructor + vectorFullFunction

# Write to cpp file
vectorCPP = open("bVector.cpp","w")
vectorCPP.writelines(vectorCPPFullLines)
vectorCPP.close()


''' Write matrix.cpp file '''

# Again, header files are set, just need methods and constructors

# Getting matrix in the eigenformat
eigenMatrix = parseArray(dimension, symbolicMatrix, vector = False)

# Writing full constructor method
matrixFullConstructor = ["\n\n\nfinalMatrix::" + constructor[5:-2] + "{\n", "   n = n;\n"]
# Writing full function method
matrixFullFunction = ["SparseMatrix<dcomplex> finalMatrix::" + function[26:-2] + "{\n"] + eigenMatrix + ["   A.makeCompressed();\n","   return A;\n", "}"]

for i in range(len(constants)):
    matrixConstantString = "{}".format(constants[i])
    matrixFullConstructor.append("   " + matrixConstantString + " = " + matrixConstantString + ";\n")
matrixFullConstructor.append("   SparseMatrix<dcomplex> A(n,n);\n")
matrixFullConstructor.append("}\n\n")

matrixCPPFullLines = ["#include \"superMatrix.h\"\n"] + headerLines + matrixFullConstructor + matrixFullFunction

matrixCPP = open("superMatrix.cpp", "w")
matrixCPP.writelines(matrixCPPFullLines)
matrixCPP.close()