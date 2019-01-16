import numpy as np
from skrypt.InitialData import InitialData

def N1(ksi, eta):
    return 0.25 * (1 - ksi) * (1 - eta)
def N2(ksi, eta):
    return 0.25 * (1 + ksi) * (1 - eta)
def N3(ksi, eta):
    return 0.25 * (1 + ksi) * (1 + eta)
def N4(ksi, eta):
    return 0.25 * (1 - ksi) * (1 + eta)

def dN_dski1(eta):
    return -1 * 0.25 * (1.0 - eta)
def dN_dski2(eta):
    return 0.25 * (1.0 - eta)
def dN_dski3(eta):
    return 0.25 * (1.0 + eta)
def dN_dski4(eta):
    return -1 * 0.25 * (1.0 + eta)

def dN_deta1(ksi):
    return -1 * 0.25 * (1.0 - ksi)
def dN_deta2(ksi):
    return -1 * 0.25 * (1.0 + ksi)
def dN_deta3(ksi):
    return 0.25 * (1.0 + ksi)
def dN_deta4(ksi):
    return 0.25 * (1.0 - ksi)

shapeFunctions = [N1, N2, N3, N4]
dN_dksiFunctions = [dN_dski1, dN_dski2, dN_dski3, dN_dski4]
dN_detaFunctions = [dN_deta1, dN_deta2, dN_deta3, dN_deta4]

class UniversalElement:
    def __init__(self):
        self.ip = 0.57735 #Wspólrzedne punktu calkowania
        #assign values to ksi & eta
        self.ksi = np.zeros((4))
        self.ksi[0] = -0.57735
        self.ksi[1] = 0.57735
        self.ksi[2] = 0.57735
        self.ksi[3] = -0.57735
        self.eta = np.zeros((4))
        self.eta[0] = -0.57735
        self.eta[1] = -0.57735
        self.eta[2] = 0.57735
        self.eta[3] = 0.57735

        self.edgeIntPoints = np.zeros((2))
        self.edgeIntPoints[0] = - 0.57735
        self.edgeIntPoints[1] = 0.57735

        self.dN_dksi = np.zeros((4,4))
        self.dN_deta = np.zeros((4,4))
        self.N = np.zeros((4,4))
        self.setN()
        self.setDerivatives()

        self.edgesN = np.zeros((4), object)
        self.setEdgesN()

    def setN(self):
        for i in range(4):
            for j in range(4):
                self.N[i, j] = shapeFunctions[j](self.ksi[i], self.eta[i])
    def setDerivatives(self):
        for i in range(4):
            for j in range(4):
                self.dN_dksi[i, j] = dN_dksiFunctions[j](self.eta[i])
                self.dN_deta[i, j] = dN_detaFunctions[j](self.ksi[i])
    def setEdgesN(self):
        bottomEdge = np.zeros((2, 4))
        rightEdge = np.zeros((2, 4))
        topEdge = np.zeros((2, 4))
        leftEdge = np.zeros((2, 4))
        for i in range(2):
            for j in range(4):
                bottomEdge[i, j] = shapeFunctions[j](self.edgeIntPoints[i], -1)
                rightEdge[i, j] = shapeFunctions[j](1, self.edgeIntPoints[i])
                topEdge[i, j] = shapeFunctions[j](self.edgeIntPoints[i], 1)
                leftEdge[i, j] = shapeFunctions[j](-1, self.edgeIntPoints[i])
        self.edgesN[0] = bottomEdge
        self.edgesN[1] = rightEdge
        self.edgesN[2] = topEdge
        self.edgesN[3] = leftEdge
        #print(self.edgesN)
    #----------------------------------PRINTING--------------------------------------
    def printN(self):
        i = 1
        for row in self.N:
            s = "N" + str(i)
            print(s, row)
            i = i + 1
        print()
    def printdN_ksi(self):
        i = 1
        for row in self.dN_dksi:
            s = "dN_ksi" + str(i)
            print(s, row)
            i = i + 1
        print()
    def printdN_eta(self):
        i = 1
        for row in self.dN_deta:
            s = "dN_eta" + str(i)
            print(s, row)
            i = i + 1
        print()





