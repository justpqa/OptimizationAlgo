# This file will help convert inequlity to equality
import string
from LPMatrix import LPMatrix

# the objective function is only a linear combination of all variables
# bugs regarding the artificial variable and the big M method: we can only deal with coefficient smaller than 10^7, 
# so that we can use a value of M = 10^9

def TexttoDict(s):
    lst = s.split(" ")
    for i in range(len(lst)):
        lst[i] = lst[i].replace(" ", "")
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

def getLPMatrix():
    # first we get the dictionary for objective function
    objDict = TexttoDict(input("Enter the objective function: "))
    
    # get the problem is max or min
    max = input("Are you doing a max/min problem: ") == "max"

    # start to create vector of variables
    vecx = set(objDict.keys())

    # get the number of constraint and create the list of constraint
    n =  int(input("Enter number of constraints: "))
    constraint = [""] * n

    # create number of surplus and excess variables in order to start create surplus and excess variable
    surplusInx = 1
    excessInx = 1
    artificialInx = 1

    # get n constraints
    curr_n = n # what if n is updated? we use this one instead
    for i in range(curr_n):
        constraintDict = TexttoDict(input("Enter constraint number {}: ".format(i+1)))
        if constraintDict["sign"] == "<=":
            varName = "s" + str(surplusInx)
            constraintDict[varName] = 1
            del constraintDict["sign"]
            constraint[i] = constraintDict
            vecx = set(list(vecx) + list(constraintDict.keys()))
            surplusInx += 1
        elif constraintDict["sign"] == ">=":
            varName = "e" + str(excessInx)
            constraintDict[varName] = -1
            varName = "a" + str(artificialInx)
            constraintDict[varName] = -10**9 if max else 10**9
            del constraintDict["sign"]
            constraint[i] = constraintDict
            vecx = set(list(vecx) + list(constraintDict.keys()))
            excessInx += 1
            artificialInx += 1
        else:
            constraintDict1, constraintDict2 = constraintDict, constraintDict
            varName = "s" + str(surplusInx)
            constraintDict1[varName] = 1
            del constraintDict1["sign"]
            varName = "e" + str(excessInx)
            constraintDict2[varName] = -1
            varName = "a" + str(artificialInx)
            constraintDict[varName] = -10**9 if max else 10**9
            del constraintDict2["sign"]
            # need to add 2 constraint instead of 1
            constraint[i] = constraintDict1
            constraint.append(constraintDict2) # add this to the end
            n += 1 # update n to get the latest number of constraint
            vecx = set(list(vecx) + list(constraintDict1.keys()) + list(constraintDict2.keys()))
            surplusInx += 1
            excessInx += 1
            artificialInx += 1
    
    # some format to make vecx the right type
    vecx.remove("lhs")
    vecx = list(vecx)
    numVar = len(vecx)

    # add artificial variable to objective function in big M method:
    for var in vecx:
        if "a" in var:
            objDict[var] = 1000000000

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

    return LPMatrix(matA, vecb, vecc, vecx, max)

if __name__ == "__main__":
    mat = getLPMatrix()
    mat.simplexAlg()
    mat.getRes()