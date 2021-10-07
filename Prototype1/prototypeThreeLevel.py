# Henry Pacheco Cachon
# Created 20 September 2021
# This is a prototype for a python function that can output a symbolic matrix and be
# read by a java program.

from nStates import *
from sympy import *
import numpy as np
from sympy.printing.c import C11CodePrinter
from re import search


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
constants = [dp0, dc0, kP, kC, v, rabi01, rabi12]

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



# Preparing C code printer 
printer = C11CodePrinter()

# Creating a 35 x 35 matrix symbol 
A = MatrixSymbol('A',dimension,dimension)

# Converting the symbolic matrix to a C array similar to the matrix symbol A
# For some reason this is a 1-d array as opposed to a 2-d array 
test = printer.doprint(symbolicMatrix,assign_to=A)

# printer.doprint() outputs a long string with each matrix element separated by a \n
# The Lines variable is a list which separates each matrix element 
Lines = test.split('\n')

# Start row counter
row = 0

# Creating formatted strings for matrix elements
# These formats are for the Eigen library for C++
matrixAssign = ".coeffRef({},{})"
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
        newString = string[:l] + matrixAssign.format(y,x) + string[r+1:r+4] + newComplex +'\n'
        newLines.append(newString)

    # Moving on to the next row
    if x == dimension - 1:
        row +=1


# Creating all the global variables within the constants list
newVar = "double {};\n"
globalVars = [newVar.format(i) for i in constants] + ['\n','\n']
function = ["SparseMatrix<dcomplex> M(delta_c, v){\n","\n"]

# Putting the matrix and global variables into one large string
allLines = globalVars + function + ["   "+ string for string in newLines] + ['\n','   A.makeCompressed();\n','\n','   return A;\n \n','}']

'''# Writing the C++ code to a text document, can write to .cpp, but it won't compile anything
testFile = open("threeLevelMatrix.txt",'w')
testFile.writelines(allLines)
testFile.close()'''

pprint(threeLevel.getVector())