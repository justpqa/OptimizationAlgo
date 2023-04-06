# This is a class for storing matrix for LP in standard form
# import possible library
# This will only deal with condition that is already inequality
from collections import Counter
from copy import deepcopy
from typing import List
import numpy as np

# currently I need to use numpy as many matrix multiplication need to be done for the interior method
class InteriorMat:

    def __init__(self, matA: List[List[int]], vecb: List[int], vecc: List[int], vecx: List[str], signvec: List[str],isMax: bool = True, z: int = 0) -> None:
        self.A = np.asmatrix(matA)
        self.b = np.transpose(np.asmatrix([vecb]))
        self.c = np.asmatrix([vecc])
        self.z = z
        self.bounded = True
        self.feasible = True
        self.x = vecx
        self.res = np.array([0] * len(self.x))
        self.s = signvec
        self.isMax = 1 if isMax else 0
        self.numSurplus = 0
        self.numExcess = 0
    
    def convertStandard(self) -> None:
        if self.isMax:
            # all of them will be treated to >=
            for i in range(len(self.s)):
                if self.s[i] == ">=":
                    self.A[i,] = -1 * self.A[i,]
                    self.b[i,] = -1 * self.b[i,]
                    self.s[i] = "<="
                self.IneqtoEq(i)
        else:
            # all of them will be treated to <=
            for i in range(len(self.s)):
                if self.s[i] == "<=":
                    self.A[i,] = -1 * self.A[i,]
                    self.b[i,] = -1 * self.b[i,]
                    self.s[i] = ">="
                self.IneqtoEq(i)
    
    def convertDual(self):
        newA = np.transpose(self.A)
        newb = np.tranpose(self.c)
        newc = np.tranpose(self.b)
        numVar = len(newA[0])
        newy = ["y" + str(i) for i in range(numVar)] 
        isMax = False if self.isMax else True
        newsign = ["<=" for i in range(len(newA))] if isMax else [">=" for i in range(len(newA))]
        return InteriorMat(newA, newb, newc, newy, newsign, isMax)
    
    def interiorPoint(self) -> None:
        # currently have to use numpy matrix on this
        # maybe try to do this without numpy => not do this as this one is too long to do
        # if A is m*n (m constraints, n variables) => X, Z is n*n, W, Y is m*m
        dualMat = self.convertDual()
        A = self.A
        At = dualMat.A
        x = np.asmatrix([[0.0001 for i in range(len(A))]])
        y = np.asmatrix([[0.0001 for i in range(len(At))]])
        n = len(x)
        m = len(y)
        b = self.b
        c = dualMat.b
        w = np.ones(len(A))
        z = np.ones(len(At))
        W = np.identity(len(A))
        Z = np.identity(len(At))
        X = np.diag(x)
        Y = np.diag(y)
        duality_gap = (np.matmul(np.transpose(c), x)[0][0]  - np.matmul(np.transpose(b), y)[0][0]) / (1 + np.abs(np.matmul(np.transpose(b), y)[0][0]))
        u = 1.1
        while duality_gap >= 10**(-6):
            # now we solve for delta_x, delta_y, delta_w, delta_z (the variable in the matrix that we will solve will work this way)
            mat = np.zeros((2*n+2*m, 2*n+2*m))
            # first we fill values for delta_x
            mat[:m, :n] = A
            mat[m+n:m+2*n, :n] = Z
            mat[m:m+n, n:n+m] = At
            mat[m+2*n:, n:n+m] = W
            mat[:m, n+m:n+2*m] = np.identity(m)
            mat[m+2*n:, n+m:n+2*m] = Y
            mat[m:m+n, n+2*m:] = -1*np.identity(m)
            mat[m+n:m+2*n, n+2*m:] = X
            em = np.asmatrix([[1 for i in range(m)]])
            en = np.asmatrix([[1 for i in range(n)]])
            rhs = np.zeros[(2*n+2*m, 1)]
            temp = np.add(b, -np.matmul(A, np.transpose(np.array([x]))))
            temp = np.add(temp, -w)
            rhs[:m,] = temp
            temp = np.add(c, -np.matmul(At, np.transpose(np.array([y]))))
            temp = np.add(temp, z)
            rhs[m:m+n, ] = temp
            ZX = np.matmul(Z, X)
            temp = np.add(u*en, -np.matmul(ZX, en))
            rhs[m+n:m+2*n] = temp
            WY = np.matmul(W, Y)
            temp = np.add(u*em, -np.matmul(WY, em))
            rhs[m+2*n:] = temp
            res = np.matmul(rhs, np.linalg.inv(mat))
            # get res and test again with new c and x
            delta_x = res[:m,]
            delta_y = res[m:m+n,]
            x = np.add(x, delta_x)
            y = np.add(y, delta_y)
            duality_gap = (np.matmul(np.transpose(c), x)[0][0]  - np.matmul(np.transpose(b), y)[0][0]) / (1 + np.abs(np.matmul(np.transpose(b), y)[0][0]))
            u /= 10
            # dunno if this gonna work though as it is 2am already!!!
        
        if duality_gap != 0:
            self.feasible = False
        self.res = x[0]
        self.z = np.matmul(np.transpose(c), x)[0][0]
        dualMat.z = np.matmul(np.transpose(b), y)[0][0]
    
    def getRes(self) -> None:
        if not self.feasible:
                print("The LP is not feasible.")
        else:
            print("Optimized value: " + str(self.z))
            for i in range(len(self.x)):
                print("Value of " + str(self.x[i]) + " is: " + str(round(self.res[i],2)))