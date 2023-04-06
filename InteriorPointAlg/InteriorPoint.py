from InteriorMat import InteriorMat
from Input import InputtoLPMatrix

def main():
    mat = InputtoLPMatrix()
    mat.interiorPoint()
    mat.getRes()

if __name__ == "__main__":
    main() 