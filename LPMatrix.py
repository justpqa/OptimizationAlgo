# This is a class for storing matrix for LP in standard form
# import possible library
# This will only deal with condition that is already inequality
from collections import Counter
from copy import deepcopy

# comment in 3/27/2023: need to check when will we have unbounded objective function
class LPMatrix:
    def __init__(self, matA, vecb, vecc, vecx, isMax = True, z = 0):
        self.A = matA
        self.b = vecb
        self.c = [0] * len(vecc)
        for i in range(len(vecc)):
            self.c[i] = -vecc[i]
        self.z = z
        self.bounded = True
        self.x = vecx
        self.isMax = 1 if isMax else 0
    
    # for the tableau already in equality spot
    # find the next spot for pivoting
    def findCol(self):
        minNeg = 0
        maxPos = 0
        inx = -1
        if self.isMax == 1:
            for i in range(len(self.c)):
                if self.c[i] < minNeg:
                    minNeg = self.c[i]
                    inx = i
        else:
            for i in range(len(self.c)):
                if self.c[i] > maxPos:
                    maxPos = self.c[i]
                    inx = i
        return inx # if inx = -1 => we will know that our function has reached optimized value
    
    def findRow(self, col):
        minRatio = 99999
        inx = -1
        for i in range(len(self.A)):
            if self.A[i][col] != 0:
                ratio = self.b[i] / self.A[i][col]
                if ratio < minRatio and ratio > 0:
                    inx = i
                    minRatio = ratio
                elif ratio < minRatio and ratio == 0 and self.A[i][col] > 0:
                    inx = i
                    minRatio = ratio
                else:
                    continue
        return inx 
    
    def productNumVec(self, num, vec):
        res = deepcopy(vec)
        for i in range(len(vec)):
            res[i] *= num
        return res
    
    def minusVec(self, vec1, vec2):
        res = [0] * len(vec1)
        for i in range(len(vec1)):
            res[i] = vec1[i] - vec2[i]
        return res
    
    def simplexPivot(self, row, col):
        self.b[row] *= 1/(self.A[row][col])
        self.A[row] = self.productNumVec(1/(self.A[row][col]), self.A[row])
        for i in range(len(self.A)):
            if i != row:
                self.b[i] -= self.A[i][col] * self.b[row]
                minus = self.productNumVec(self.A[i][col], self.A[row])
                self.A[i] = self.minusVec(self.A[i], minus)
        self.z -= self.c[col] * self.b[row]
        minus = self.productNumVec(self.c[col], self.A[row])
        self.c = self.minusVec(self.c, minus)
    
    def simplexAlg(self):
        while self.findCol() != -1:
            c = self.findCol()
            r = self.findRow(c)
            if r == -1:
                self.bounded = False
                break
            else:
                self.simplexPivot(r, c)
    
    def getRes(self):
        if self.bounded:
            print("Optimized value: " + str(self.z))
            for i in range(len(self.x)):
                arr = [self.A[j][i] for j in range(len(self.A))]
                d = Counter(arr)
                if d[1] == 1 and d[0] == len(arr) - 1:
                    inx = arr.index(1)
                    val = self.b[inx]
                    print("Value of " + str(self.x[i]) + " is: " + str(round(val,2)))
                else:
                    print("Value of " + str(self.x[i]) + " is: 0") 
        else:
            print("The objective function is unbounded.")  