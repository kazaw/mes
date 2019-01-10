from skrypt.Solver import Solver
import numpy as np
import shutil
import os

def removeOutput():
    shutil.rmtree('output')
    path = os.path.dirname(__file__)
    path = os.path.join(path, 'output')
    os.mkdir(path)

def main():
    removeOutput()
    solver = Solver("data.json")
    solver.execute(10)

if __name__ == "__main__":
    main()