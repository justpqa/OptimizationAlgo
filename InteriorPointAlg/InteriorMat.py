# This is a class for storing matrix for LP in standard form
# import possible library
# This will only deal with condition that is already inequality
from collections import Counter
import MatMethod as mm
import math
from typing import List
import numpy as np

# currently I need to use numpy as many matrix multiplication need to be done for the interior method
class InteriorMat:

    def __init__(self, matA: List[List[int]], vecb: List[int], vecc: List[int], vecx: List[str], signvec: List[str],isMax: bool = True, z: int = 0) -> None:
        self.numVar = len(matA[0]) # this is n
        self.numCon = len(matA) # this is m
        self.A = matA
        self.b = [[b] for b in vecb]
        self.c = [[c] for c in vecc] # this is currently a column vector
        self.z = z
        self.x = vecx
        self.res = [0] * len(self.x)
        self.s = signvec
        self.bounded = True
        self.feasible = True
        self.isMax = 1 if isMax else 0
        self.numSurplus = 0
        self.numExcess = 0
    
    def standardize(self) -> None:
        if self.isMax:
            for i in range(len(self.s)):
                if self.s[i] == ">=":
                    self.A[i] = mm.productNumVec(-1, self.A[i])
                    self.b[i] = mm.productNumVec(-1, self.b[i])
                    self.s[i] = "<="
        else:
            for i in range(len(self.s)):
                if self.s[i] == "<=":
                    self.A[i] = mm.productNumVec(-1, self.A[i])
                    self.b[i] = mm.productNumVec(-1, self.b[i])
                    self.s[i] = ">="
    
    # try to implement Karmakar method instead
    def Karmarkarize(self):
        l = self.numCon + self.numVar + 1
        # Karmakarize an LP problem to the form: find x such that cx <= 0, x >= 0, Ax = 0, 1x = 1
        # initial solution in new problem is u = (1 1 .... 1) and then convert to v = (1/n 1/n .... 1/n)
        A = self.A # m*n
        At = mm.transpose(self.A) # n*m
        b = self.b # m*1
        bt = mm.transpose(self.b) # 1*m
        c = self.c #n*1
        ct = mm.transpose(self.c) # 1*n
        Im = mm.identity(self.numCon) # m*m
        In = mm.identity(self.numVar) # n*n
        ones_m = [[1] for i in range(self.numCon)] # 1*m
        ones_n = [[1] for i in range(self.numVar)] # 1*n

        # make c_hat for our Karmarkar method
        c_hat = [0 for i in range(2*l)]
        c_hat[-2] = 1

        # make A_hat for our Karmakar method (need to make 3 big row)
        big_row_1 = ct
        big_row_1 = mm.horz_concat(big_row_1, mm.productNumMat(-1, bt))
        temp = mm.addMat(mm.productNumMat(-1, mm.matmul(ct, ones_n)), mm.matmul(bt, ones_m))
        big_row_1 = mm.horz_concat(big_row_1, temp)
        big_row_1 = mm.horz_concat(big_row_1, [[0]])

        big_row_2 = A
        big_row_2 = mm.horz_concat(big_row_2, [[0 for j in range(self.numCon)] for i in range(self.numCon)])
        big_row_2 = mm.horz_concat(big_row_2, Im)
        big_row_2 = mm.horz_concat(big_row_2, [[0 for j in range(self.numVar)] for i in range(self.numCon)])
        temp = mm.minusMat(b, mm.productNumMat(-1, mm.matmul(A, ones_n)))
        big_row_2 = mm.horz_concat(big_row_2, temp)
        big_row_2 = mm.horz_concat(big_row_2, mm.productNumMat(-1, b))

        big_row_3 = [[0 for j in range(self.numVar)] for i in range(self.numVar)]
        big_row_3 = mm.horz_concat(big_row_3, At)
        big_row_3 = mm.horz_concat(big_row_3, [[0 for j in range(self.numCon)] for i in range(self.numVar)])
        big_row_3 = mm.horz_concat(big_row_3, mm.productNumMat(-1, In))
        temp = mm.minusMat(c, mm.productNumMat(-1, mm.matmul(At, ones_m)))
        big_row_3 = mm.horz_concat(big_row_3, temp)
        big_row_3 = mm.horz_concat(big_row_3, mm.productNumMat(-1, c))

        A_hat = big_row_1
        A_hat = mm.vert_concat(A_hat, big_row_2)
        A_hat = mm.vert_concat(A_hat, big_row_3)

        # make b_hat
        b_hat = [0 for i in range(l)]

        # make initial u
        u = ["u" + str(i) for i in range(1, 2*l+1)]
        res = [1 / (l) for i in range(2*l)]


        sign_u = ["=" for i in range(l)]

        final = InteriorMat(A_hat, b_hat, c_hat, u, sign_u, isMax = False)
        final.res = [1 / l for i in range(2*l)]

        return final 
    
    def solveK(self): # solve after Karmarkarize the problem
        q = 100 # q = ?
        alpha = 1/4
        r = 1 / math.sqrt(self.numVar * (self.numVar - 1))
        threshold = 1 / 2**(-q)
        x_new = self.res
        ct = mm.transpose(self.c)
        At = mm.transpose(self.A)
        while mm.matmul(ct, x_new)[0][0] / mm.matmul(ct, self.x)[0][0] >= threshold:
            D_new = mm.diag(self.res) 
            mid_big_term = D_new
            # need inverse function to continue
            # mid_big_term = DAt(AD^2At)^{-1}AD 
            # left_term = 1/n * 1_n * (1_n)T
            # p = (I - mid_big_term - left_term)Dct
            # norm_p = p / ||p||
            # z_new = 1 - alpha*r*(norm_p)
            # x_new = (1t*z_new)^{-1}Dz_new
                
        return 0
