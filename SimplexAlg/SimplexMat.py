# This is a class for storing matrix for LP in standard form
# import possible library
# This will only deal with condition that is already inequality
from collections import Counter
from copy import deepcopy
from typing import List

# comment in 3/27/2023: need to check when will we have unbounded objective function
class SimplexMat:
    def __init__(self, matA: List[List[int]], vecb: List[int], vecc: List[int], vecx: List[str], signvec: List[str], isMax: bool = True, z: int = 0) -> None:
        self.A = matA
        self.b = vecb
        self.c = [0] * len(vecc)
        for i in range(len(vecc)):
            self.c[i] = -vecc[i]
        self.z = z
        self.bounded = True
        self.feasible = True
        self.x = vecx
        self.isMax = 1 if isMax else 0
        self.s = signvec
        self.numSurplus = 0
        self.numExcess = 0
        self.numArtificial = 0
    
    # for the tableau already in equality spot
    # find the next spot for pivoting
    def findCol(self) -> int:
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
    
    def findRow(self, col: int) -> int:
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
    
    def productNumVec(self, num: int, vec: List[int]) -> List[int]:
        res = deepcopy(vec)
        for i in range(len(vec)):
            res[i] *= num
        return res
    
    def minusVec(self, vec1: List[int], vec2: List[int]) -> List[int]:
        res = [0] * len(vec1)
        for i in range(len(vec1)):
            res[i] = vec1[i] - vec2[i]
        return res
    
    def addVec(self, vec1: List[int], vec2: List[int]) -> List[int]:
        res = [0] * len(vec1)
        for i in range(len(vec1)):
            res[i] = vec1[i] + vec2[i]
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
            self.numArtificial += 1
            varName = "a" + str(self.numArtificial)
            self.x.append(varName)
            if self.isMax:
                self.c.append(10000000)
            else:
                self.c.append(-10000000)
            for i in range(len(self.A)):
                if i == inx:
                    self.A[i].append(-1)
                    self.A[i].append(1)
                else:
                    self.A[i].append(0)
                    self.A[i].append(0)
            temp = self.productNumVec(-10000000, self.A[inx]) if self.isMax else self.productNumVec(10000000, self.A[inx])
            self.c = self.addVec(self.c, temp)
            self.z += -10000000 * self.b[inx] if self.isMax else 10000000 * self.b[inx]
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

    def simplexPivot(self, row: int, col: int) -> None:
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
    
    def checkFeasible(self) -> bool:
        for i in range(len(self.x)):
            if "a" in self.x[i]:
                arr = [self.A[j][i] for j in range(len(self.A))]
                d = Counter(arr)
                if d[1] == 1 and d[0] == len(arr) - 1:
                    inx = arr.index(1)
                    if self.b[inx] == 0:
                        continue
                    else:
                        return False
        return True
    
    def simplexAlg(self) -> None:
        self.convertStandard()
        while self.findCol() != -1:
            c = self.findCol()
            r = self.findRow(c)
            if r == -1:
                self.bounded = False
                break
            else:
                self.simplexPivot(r, c)
        self.feasible = self.checkFeasible()
    

    def getRes(self) -> None:
        if not self.bounded:
            print("The objective function is unbounded.")
        else:
            if not self.feasible:
                print("The LP is not feasible.")
            else:
                print("Optimized value: " + str(self.z))
                for i in range(len(self.x)):
                    if "x" in self.x[i]:
                        arr = [self.A[j][i] for j in range(len(self.A))]
                        d = Counter(arr)
                        if d[1] == 1 and d[0] == len(arr) - 1:
                            inx = arr.index(1)
                            val = self.b[inx]
                            print("Value of " + str(self.x[i]) + " is: " + str(round(val,2)))
                        else:
                            print("Value of " + str(self.x[i]) + " is: 0") 
