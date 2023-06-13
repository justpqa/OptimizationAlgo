# This is a class for storing matrix for LP in standard form
# import possible library
# This will only deal with condition that is already inequality
from collections import Counter
import MatMethod as mm
import numpy as np # will try to remove this later
import math
from typing import List

# currently I need to use numpy as many matrix multiplication need to be done for the interior method
class InteriorMat:

    def __init__(self, matA: List[List[int]], vecb: List[int], vecc: List[int], vecx: List[str], signvec: List[str],isMax: bool = True, z: int = 0) -> None:
        self.numVar = len(matA[0]) # this is n
        self.numCon = len(matA) # this is m
        self.A = np.array(matA) #matA
        self.b = np.array([vecb]).T #[[b] for b in vecb]
        self.c = np.array([vecb]).T # this is currently a column vector
        self.z = z
        self.x = vecx
        self.res = np.array([0] * len(self.x))
        self.s = signvec
        self.fandb = True
        self.isMax = 1 if isMax else 0
        self.M = None
    
    def standardize(self) -> None:
        if self.isMax:
            for i in range(len(self.s)):
                if self.s[i] == ">=":
                    self.A[i, :] = -1 * self.A[i, :]
                    self.b[i, :] = -1 * self.b[i, :]
                    self.s[i] = "<="
        else:
            for i in range(len(self.s)):
                if self.s[i] == "<=":
                    self.A[i, :] = -1 * self.A[i, :]
                    self.b[i, :] = -1 * self.b[i, :]
                    self.s[i] = ">="
    
    # try to implement Karmakar method instead
    def Karmarkarize(self):
        l = self.numCon + self.numVar + 1
        M = 2**l
        # Karmakarize an LP problem to the form: find x such that cx <= 0, x >= 0, Ax = 0, 1x = 1
        # initial solution in new problem is u = (1 1 .... 1) and then convert to v = (1/n 1/n .... 1/n)
        A = self.A # m*n
        At = A.T
        b = self.b # m*1
        bt = b.T
        c = self.c #n*1
        ct = c.T # 1*n
        Im = np.identity(self.numCon) # m*m
        In = np.identity(self.numVar) # n*n
        ones_m = np.array([[1] for i in range(self.numCon)]) # 1*m
        ones_n = np.array([[1] for i in range(self.numVar)]) # 1*n

        # make c_hat for our Karmarkar method
        c_hat = np.array([[0 for i in range(2*l-1)]]).T
        c_hat[-1, 0] = 1

        # make A_hat for our Karmakar method (need to make 3 big row)
        big_row_1 = ct
        big_row_1 = np.hstack(big_row_1, -1*bt)
        big_row_1 = np.hstack(big_row_1, np.array([[0 for i in range(self.numCon)]]))
        big_row_1 = np.hstack(big_row_1, np.array([[0 for i in range(self.numVar)]]))
        temp = (-1) * np.matmul(ct, ones_n) + np.matmul(bt, ones_m)
        big_row_1 = np.hstack(big_row_1, temp)   
        big_row_1 = np.hstack(big_row_1, np.array([[0, 0]]))

        big_row_2 = A
        big_row_2 = np.hstack(big_row_2, np.array([[0 for j in range(self.numCon)] for i in range(self.numCon)]))
        big_row_2 = np.hstack(big_row_2, Im)
        big_row_2 = np.hstack(big_row_2, np.array([[0 for j in range(self.numVar)] for i in range(self.numCon)]))
        temp = b - np.matmul(A, ones_n)
        big_row_2 = np.hstack(big_row_2, temp)
        big_row_2 = np.hstack(big_row_2, np.array([[0] for i in range(self.numCon)])) # d1
        big_row_2 = np.hstack(big_row_2, (-1) * b) # d2

        big_row_3 = np.array([[0 for j in range(self.numVar)] for i in range(self.numVar)])
        big_row_3 = np.hstack(big_row_3, At)
        big_row_3 = np.hstack(big_row_3, np.array([[0 for j in range(self.numCon)] for i in range(self.numVar)]))
        big_row_3 = np.hstack(big_row_3, (-1) * In)
        temp = c - mm.matmul(At, ones_m) + ones_n
        big_row_3 = np.hstack(big_row_3, temp)
        big_row_3 = np.hstack(big_row_3, np.array([[0] for i in range(self.numVar)])) # d1
        big_row_3 = np.hstack(big_row_3, (-1)*c) # d2
        
        big_row_4 = np.array([[1 for i in range(2*l + 1)]])
        big_row_4[0, -1] = -M
        
        big_row_5 = np.array([[1 for i in range(2*l + 1)]])

        A_hat = big_row_1
        A_hat = np.vstack(A_hat, big_row_2)
        A_hat = np.vstack(A_hat, big_row_3)
        A_hat = np.vstack(A_hat, big_row_4)
        A_hat = np.vstack(A_hat, big_row_5)

        # make b_hat
        b_hat = np.array([[0 for i in range(l)]]).T

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
        x_new = np.array([[r] for r in self.res])
        ct = self.c.T
        At = self.A.T
        ones_n = np.array([[1] for i in range(len(self.res))])
        # print(mm.matmul(self.A, ones_n)) does not satisfy Ax = 0 for starting x
        while np.matmul(ct, x_new)[0,0] >= threshold:
            D_new = mm.diag(x_new) 
            # calculate p
            # calculate middle part first
            middle_part = self.A @ D_new @ D_new @ At
            middle_part = np.linalg.pinv(middle_part) # problem at the inverse method
            middle_part = D_new @ At @ middle_part
            middle_part = middle_part @ self.A @ D_new
            # calculate left_part
            left_part = np.matmul(ones_n, ones_n.T)
            left_part = (1/self.numVar) * left_part
            # calculate right_part
            right_part = np.identity(self.numVar)
            # find p from 3 parts
            p = right_part - middle_part
            p = p - left_part
            p = p @ D_new @ self.c
            p = p / np.linalg.norm(p)
            z_new = ones_n - (alpha*r) * p
            temp = ones_n.T @ D_new @ z_new
            x_new = (D_new / temp[0, 0]) @ z_new
        # need another step to convert to original variable
        if np.matmul(ct, x_new)[0, 0] == 0:
            self.fandb = True
            self.z = np.matmul(ct, x_new)[0][0]
            #self.res = mm.productNumVec([x_new[i][0] for i in range(len(self.res))] / self.M) how to deal with this
        else:
            self.fandb = False
    
    def getRes(self):
        if not self.fandb:
            print("The problem is either infeasible or unbounded.")
        else:
            for i in range(len(self.res)):
                print(self.x[i] + ": " + str(self.res[i]))
