from typing import List

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
    
def productNumMat(num: int, mat: List[List[int]]) -> List[int]:
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

def matrixMultiplication(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    try:
        return [[sum(m1 * m2 for m1, m2 in zip(mat1_row, mat2_col)) for mat2_col in zip(*mat2)] for mat1_row in mat1]
    except:
        print("Error in matrix multiplication, need to check again the shape and input of two matrices")
        return [[0]]
    
def vert_concat(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    if len(mat1) != len(mat2):
        print("The two matrix does not have the same number of rows. ")
        return [[0]]
    else:
        m = len(mat1)
        res = [[] for i in range(m)]
        for i in range(m):
            res[i] = mat1[i] + mat2[i]
        return res
    
def horz_concat(mat1: List[List[int]], mat2: List[List[int]]) -> List[List[int]]:
    if len(mat1[0]) != len(mat2[0]):
        print("The two matrix does not have the same number of columns. ")
        return [[0]]
    else:
        return mat1 + mat2 # use list comprehension directly