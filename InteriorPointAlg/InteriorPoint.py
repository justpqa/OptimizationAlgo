from InteriorMat import InteriorMat
from Input import InputtoLPMatrix

def main():
    mat = InputtoLPMatrix()
    print("1st checkpoint")
    mat.standardize()
    print("2nd checkpoint")
    mat = mat.Karmarkarize()
    print("3rd checkpoint")
    mat.solveK()
    print("4th checkpoint")
    mat.getRes()
    print("5th checkpoint")

if __name__ == "__main__":
    main() 