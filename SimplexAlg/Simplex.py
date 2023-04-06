from SimplexMat import SimplexMat
from Input import InputtoLPMatrix

def main():
    mat = InputtoLPMatrix()
    mat.simplexAlg()
    mat.getRes()

if __name__ == "__main__":
    main()