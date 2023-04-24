from typing import List
from copy import deepcopy

def productNumVec(num: int, vec: List[int]) -> List[int]:
    return [num*vec[i] for i in range(len(vec))]

def addVec(vec1: List[int], vec2: List[int]) -> List[int]:
    try:
        return [vec1[i] + vec2[i] for i in range(len(vec1))]
    except:
        print("The size of two vector or the input of vector is having a problem.")
        return [0]
    
def minusVec(vec1: List[int], vec2: List[int]) -> List[int]:
    try:
        return [vec1[i] - vec2[i] for i in range(len(vec1))]
    except:
        print("The size of two vector or the input of vector is having a problem.")
        return [0]
    
def productNumMat(num: float, mat: List[List[int]]) -> List[float]:
    return [[num*mat[i][j] for j in range(len(mat[0]))] for i in range(len(mat))]

def addMat(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    try:
        return [[mat1[i][j] + mat2[i][j] for j in range(len(mat1[0]))] for i in range(len(mat1))]
    except:
        print("The size of two matrix or the input of matrix is having a problem.")
        return [[0]]
    
def minusMat(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    try:
        return [[mat1[i][j] - mat2[i][j] for j in range(len(mat1[0]))] for i in range(len(mat1))]
    except:
        print("The size of two matrix or the input of matrix is having a problem.")
        return [[0]]

def transpose(mat: List[List[int]]) -> List[List[int]]:
    return [[row[i] for row in mat] for i in range(len(mat[0]))]

def matmul(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    try:
        return [[sum(m1 * m2 for m1, m2 in zip(mat1_row, mat2_col)) for mat2_col in zip(*mat2)] for mat1_row in mat1]
    except:
        print("Error in matrix multiplication, need to check again the shape and input of two matrices")
        return [[0]]
    
def vert_concat(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    if len(mat1[0]) != len(mat2[0]):
        print("The two matrix does not have the same number of columns. ")
        return [[0]]
    else:
        return mat1 + mat2 # use list comprehension directly
    
def horz_concat(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    if len(mat1) != len(mat2):
        print("The two matrix does not have the same number of rows. ")
        return [[0]]
    else:
        m = len(mat1)
        res = [[] for i in range(m)]
        for i in range(m):
            res[i] = mat1[i] + mat2[i]
        return res

def identity(n: int) -> List[List[int]]:
    res = [[0 for j in range(n)] for i in range(n)]
    for i in range(n):
        res[i][i] = 1
    return res

def diag(vec: List[int]) -> List[List[int]]:
    res = [[0 for j in range(len(vec))] for i in range(len(vec))]
    for i in range(len(vec)):
        res[i][i] = vec[i]
    return res

def normalize(colVec: List[List[int]]):
    ss = 0
    for i in range(len(colVec)):
        for num in colVec[i]:
            ss += num**2
    for i in range(len(colVec)):
        for j in range(len(colVec[0])):
            colVec[i][j] /= ss
    return colVec

def inverse(mat: List[List[int]]) -> List[List[int]]:
    if len(mat) != len(mat[0]):
        print("Matrix need to be a square matrix.")
        return [[0]]
    else:
        n = len(mat)
        new_mat = deepcopy(mat)
        In = identity(n)
        indices = list(range(n)) # to allow flexible row referencing ***
        for fd in range(n): # fd stands for focus diagonal
            fdScaler = 1.0 / mat[fd][fd]
            # FIRST: scale fd row with fd inverse. 
            for j in range(n): # Use j to indicate column looping.
                new_mat[fd][j] *= fdScaler
                In[fd][j] *= fdScaler
            # SECOND: operate on all rows except fd row as follows:
            for i in indices[0:fd] + indices[fd+1:]: 
                # *** skip row with fd in it.
                crScaler = new_mat[i][fd] # cr stands for "current row".
                for j in range(n): 
                    # cr - crScaler * fdRow, but one element at a time.
                    new_mat[i][j] = new_mat[i][j] - crScaler * new_mat[fd][j]
                    In[i][j] = In[i][j] - crScaler * In[fd][j]
        if new_mat != identity(n):
            print("The matrix does not have inverse.")
            return [[0]]
        else:
            return new_mat
        
def shape(mat: List[List[int]]):
    print(str(len(mat)) + "*" + str(len(mat[0])))
