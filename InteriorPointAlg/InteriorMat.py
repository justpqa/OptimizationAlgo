# This is a class for storing matrix for LP in standard form
# import possible library
# This will only deal with condition that is already inequality
from collections import Counter
from copy import deepcopy
from typing import List
import numpy as np

# comment in 3/27/2023: need to check when will we have unbounded objective function
class InteriorMat:

    def __init__(self, matA: List[List[int]], vecb: List[int], vecc: List[int], vecx: List[str], signvec: List[str],isMax: bool = True, z: int = 0) -> None:
        self.A = matA
        self.b = vecb
        self.c = vecc
        self.z = z
        self.bounded = True
        self.feasible = True
        self.x = vecx
        self.s = signvec
        self.isMax = 1 if isMax else 0
        self.numSurplus = 0
        self.numExcess = 0
    
    def productNumVec(self, num: int, vec: List[int]) -> List[int]:
        res = deepcopy(vec)
        for i in range(len(vec)):
            res[i] *= num
        return res
    
    # modify from ineq to eq at condition inx (0-indexed array)
    def IneqtoEq(self, inx) -> None:
        # only add surplus and excess
        if self.s[inx] == "<=":
            self.numSurplus += 1
            varName = "s" + str(self.numSurplus)
            self.x.append(varName)
            self.c.append(0)
            for i in range(len(self.A)):
                if i == inx:
                    self.A[i].append(1)
                else:
                    self.A[i].append(0)
            self.s[inx] = "="
        else:
            self.numExcess += 1
            varName = "e" + str(self.numExcess)
            self.x.append(varName)
            self.c.append(0)
            for i in range(len(self.A)):
                if i == inx:
                    self.A[i].append(-1)
                else:
                    self.A[i].append(0)
            self.s[inx] = "="
    
    def convertStandard(self) -> None:
        if self.isMax:
            # all of them will be treated to >=
            for i in range(len(self.s)):
                if self.s[i] == ">=":
                    self.A[i] = self.productNumVec(-1, self.A[i])
                    self.b[i] *= -1
                    self.s[i] = "<="
                self.IneqtoEq(i)
        else:
            # all of them will be treated to <=
            for i in range(len(self.s)):
                if self.s[i] == "<=":
                    self.A[i] = self.productNumVec(-1, self.A[i])
                    self.b[i] *= -1
                    self.s[i] = ">="
                self.IneqtoEq(i)
    
    def transpose(self, mat):
        # iterate over 2d list mat to the length of an item
        res = []
        for i in range(len(mat[0])):
            # print(i)
            row =[]
            for item in mat:
                # appending to new list with values and index positions
                # i contains index position and item contains values
                row.append(item[i])
            res.append(row)
        return res
    
    def convertDual(self):
        newA = self.transpose(self.A)
        newb = self.c
        newc = self.b
        numVar = len(newA[0])
        newy = ["y" + str(i) for i in range(numVar)] 
        isMax = False if self.isMax else True
        newsign = ["<=" for i in range(len(newA))] if isMax else [">=" for i in range(len(newA))]
        return InteriorMat(newA, newb, newc, newy, newsign, isMax)
    
    def interiorPoint(self):
        # currently have to use numpy matrix on this
        # maybe try to do this without numpy
        # if A is m*n (m constraints, n variables) => X, Z is n*n, W, Y is m*m
        dualMat = self.convertDual()
        A = np.asmatrix(self.A)
        At = np.asmatrix(dualMat.A)
        x = np.transpose(np.asmatrix([[0.001 for i in range(len(A[0]))]]))
        y = np.transpose(np.asmatrix([[0.001 for i in range(len(At[0]))]]))
        n = len(x)
        m = len(y)
        b = np.transpose(np.asmatrix([self.b]))
        c = np.transpose(np.asmatrix([self.c]))
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
            # dunno if this gonna work though as it is 2am already!!!

