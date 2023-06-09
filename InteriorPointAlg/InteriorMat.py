# This is a class for storing matrix for LP in standard form
# import possible library
# This will only deal with condition that is already inequality
from collections import Counter
import MatMethod as mm
import math
from typing import List

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
        self.fandb = True
        self.isMax = 1 if isMax else 0
        self.M = None
    
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
        M = 2**l
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
        c_hat = [0 for i in range(2*l-1)]
        c_hat[-1] = 1

        # make A_hat for our Karmakar method (need to make 3 big row)
        big_row_1 = ct
        big_row_1 = mm.horz_concat(big_row_1, mm.productNumMat(-1, bt))
        big_row_1 = mm.horz_concat(big_row_1, [[0 for i in range(self.numCon)]])
        big_row_1 = mm.horz_concat(big_row_1, [[0 for i in range(self.numVar)]])
        temp = mm.addMat(mm.productNumMat(-1, mm.matmul(ct, ones_n)), mm.matmul(bt, ones_m))
        big_row_1 = mm.horz_concat(big_row_1, temp)   
        big_row_1 = mm.horz_concat(big_row_1, [[0, 0]])

        big_row_2 = A
        big_row_2 = mm.horz_concat(big_row_2, [[0 for j in range(self.numCon)] for i in range(self.numCon)])
        big_row_2 = mm.horz_concat(big_row_2, Im)
        big_row_2 = mm.horz_concat(big_row_2, [[0 for j in range(self.numVar)] for i in range(self.numCon)])
        temp = mm.minusMat(b, mm.matmul(A, ones_n))
        big_row_2 = mm.horz_concat(big_row_2, temp)
        big_row_2 = mm.horz_concat(big_row_2, [[0] for i in range(self.numCon)]) # d1
        big_row_2 = mm.horz_concat(big_row_2, mm.productNumMat(-1, b)) # d2

        big_row_3 = [[0 for j in range(self.numVar)] for i in range(self.numVar)]
        big_row_3 = mm.horz_concat(big_row_3, At)
        big_row_3 = mm.horz_concat(big_row_3, [[0 for j in range(self.numCon)] for i in range(self.numVar)])
        big_row_3 = mm.horz_concat(big_row_3, mm.productNumMat(-1, In))
        temp = mm.addMat(mm.minusMat(c,  mm.matmul(At, ones_m)), ones_n)
        big_row_3 = mm.horz_concat(big_row_3, temp)
        big_row_3 = mm.horz_concat(big_row_3, [[0] for i in range(self.numVar)]) # d1
        big_row_3 = mm.horz_concat(big_row_3, mm.productNumMat(-1, c)) # d2
        
        big_row_4 = [[1 for i in range(2*l)]]
        big_row_4[-1] = -M
        
        big_row_5 = [[1 for i in range(2*l)]]

        A_hat = big_row_1
        A_hat = mm.vert_concat(A_hat, big_row_2)
        A_hat = mm.vert_concat(A_hat, big_row_3)
        A_hat = mm.vert_concat(A_hat, big_row_4)
        A_hat = mm.vert_concat(A_hat, big_row_5)

        # make b_hat
        b_hat = [0 for i in range(l)]

        # make initial u
        u = ["u" + str(i) for i in range(1, 2*l + 1)]
        res = [1/(2*l) for i in range(2*l)]


        sign_u = ["=" for i in range(l + 2)]

        final = InteriorMat(A_hat, b_hat, c_hat, u, sign_u, isMax = False)
        final.res = res
        
        self.M = M

        return final 
    
    def solveK(self) -> None: # solve after Karmarkarize the problem
        alpha = 1/4
        r = math.sqrt(self.numVar /(self.numVar-1))
        threshold = 10**(-100)
        x_new = [[r] for r in self.res]
        ct = mm.transpose(self.c)
        At = mm.transpose(self.A)
        ones_n = [[1] for i in range(len(self.res))]
        # print(mm.matmul(self.A, ones_n)) does not satisfy Ax = 0 for starting x
        while mm.matmul(ct, x_new)[0][0] >= threshold:
            D_new = mm.diag(x_new) 
            # calculate p
            # calculate middle part first
            middle_part = mm.matmul(mm.matmul(mm.matmul(self.A, D_new), D_new), At)
            middle_part = mm.inverse(middle_part)
            middle_part = mm.matmul(mm.matmul(D_new, At), middle_part)
            middle_part = mm.matmul(mm.matmul(middle_part, self.A), D_new)
            # calculate left_part
            left_part = mm.matmul(ones_n, mm.transpose(ones_n))
            left_part = mm.productNumMat(1/self.numVar, left_part)
            # calculate right_part
            right_part = mm.identity(self.numVar)
            # find p from 3 parts
            p = mm.minusMat(right_part, middle_part)
            p = mm.minusMat(p, left_part)
            p = mm.matmul(mm.matmul(p, D_new), self.c)
            p = mm.normalize(p)
            z_new = mm.minusMat(ones_n, mm.productNumMat(alpha*r, p))
            temp = mm.matmul(mm.matmul(mm.transpose(ones_n), D_new), z_new)
            x_new = mm.matmul(mm.productNumMat(1/temp[0][0], D_new), z_new)
        # need another step to convert to original variable
        if mm.matmul(ct, x_new)[0][0] == 0:
            self.fandb = True
            self.z = mm.matmul(ct, x_new)[0][0]
            self.res = mm.productNumVec([x_new[i][0] for i in range(len(self.res))] / self.M)
        else:
            self.fandb = False
    
    def getRes(self):
        if not self.fandb:
            print("The problem is either infeasible or unbounded.")
        else:
            for i in range(len(self.res)):
                print(self.x[i] + ": " + str(self.res[i]))
