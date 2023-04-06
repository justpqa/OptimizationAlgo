# This file will help convert inequlity to equality
import string
from InteriorMat import InteriorMat

# the objective function is only a linear combination of all variables
# bugs regarding the artificial variable and the big M method: we can only deal with coefficient smaller than 10^7, 
# so that we can use a value of M = 10^9

def TexttoDict(s: str) -> dict:
    # how to separate into same list as lst
    start = end = 0
    lst = []
    while end < len(s):
        if s[end] in ["+", "-", "="]:
            temp1 = s[start:end]
            temp2 = s[end]
            temp1 = temp1.replace(" ", "")
            lst.append(temp1)
            lst.append(temp2)
            end += 1
            start = end
        elif s[end] in [">", "<"]:
            temp1 = s[start:end]
            temp2 = s[end:end+2]
            temp1 = temp1.replace(" ", "")
            lst.append(temp1)
            lst.append(temp2)
            end += 2
            start = end
        elif end == len(s) - 1:
            temp = s[start:]
            temp = temp.replace(" ", "")
            lst.append(temp)
            end += 1
        else:
            end += 1
    d = {}
    if (">=" in lst) or ("<=" in lst) or ("=" in lst): 
        d["lhs"] = int(lst.pop())
        d["sign"] = lst.pop()
    for i in range(len(lst)):
        if lst[i] not in ["+", "-"]:
            for j in range(len(lst[i])):
                if lst[i][j] in string.ascii_letters:
                    var = lst[i][j:]
                    val = lst[i][:j]
                    if val == "":
                        d[var] = 1
                    elif val == "-":
                        d[var] = -1
                    else:
                        d[var] = int(val)
        elif lst[i] == "-":
            lst[i+1] = "-" + lst[i+1]
        else:
            continue
    return d

def InputtoLPMatrix() -> InteriorMat:
    # first we get the dictionary for objective function
    objDict = TexttoDict(input("Enter the objective function: "))
    
    # get the problem is max or min
    max = input("Are you doing a max/min problem: ") == "max"

    # start to create vector of variables
    vecx = set(objDict.keys())

    # get the number of constraint and create the list of constraint
    n =  int(input("Enter number of constraints: "))
    constraint = [""] * n

    # initialize sign vector
    signvec = [""] * n

    curr_n = n
    # get n constraints
    for i in range(curr_n):
        constraintDict = TexttoDict(input("Enter constraint number {}: ".format(i+1)))
        if constraintDict["sign"] == "=":
            constraintDict1 = constraintDict2 = constraintDict
            signvec[i] = ">="
            del constraintDict1["sign"]
            constraint[i] = constraintDict1
            signvec.append("<=")
            del constraintDict2["sign"]
            constraint.append(constraintDict2)
            n += 1
            vecx = set(list(vecx) + list(constraintDict1.keys()))
        else:
            signvec[i] = constraintDict["sign"]
            del constraintDict["sign"]
            constraint[i] = constraintDict
            vecx = set(list(vecx) + list(constraintDict.keys()))

    # some format to make vecx the right type
    vecx.remove("lhs")
    vecx = list(vecx)
    numVar = len(vecx)

    # add artificial variable to objective function in big M method:
    for var in vecx:
        if "a" in var:
            objDict[var] = -1000000000 if max else 1000000000

    # with the number of constraints and the number of variables, we create the component for the LPMatrix
    vecc = [0] * numVar
    vecb = [0] * n
    matA = [[0 for x in range(numVar)] for y in range(n)]

    # fill vecc
    for key in objDict:
        inx = vecx.index(key)
        vecc[inx] = objDict[key]
    
    # fill matA and vecb
    for i in range(n):
        for key in constraint[i]:
            if key == "lhs":
                vecb[i] = int(constraint[i][key])
            else:
                inx = vecx.index(key)
                matA[i][inx] = constraint[i][key]

    return InteriorMat(matA, vecb, vecc, vecx, signvec, max)
