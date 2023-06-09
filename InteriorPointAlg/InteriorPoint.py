from InteriorMat import InteriorMat
from Input import InputtoLPMatrix

def main():
    mat = InputtoLPMatrix()
    mat.standardize()
    mat = mat.Karmarkarize()
    mat.solveK()
    mat.getRes()

if __name__ == "__main__":
    main() 