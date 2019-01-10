from skrypt.UniversalElement import UniversalElement
from skrypt.Grid import Grid
import numpy as np
import math
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.collections



class Solver:
    def __init__(self, filename):
        self.ipa = 4 #integralPointsAmount
        self.grid = Grid(filename)
        self.universalElement = UniversalElement()
        self.globalVectorP = np.zeros(len(self.grid.nodes))
        self.globalMatrixH = np.zeros((len(self.grid.nodes),len(self.grid.nodes)))
        self.globalMatrixC = np.zeros((len(self.grid.nodes),len(self.grid.nodes)))
        self.MatrixMinMaxTemp = []
        #self.fig, self.axs = plt.subplots(len(self.grid.nodes), 1, figsize=(5, 5)) #do stworzenia subplots
    def execute(self, howManyIterations):
        for i in range(howManyIterations):
            self.clearGlobal()
            print("Iteracja:", i, "Sekundy:", (i+1) * self.grid.initData.simulationTimeStep)
            for element in self.grid.elements:
                #print(element)
                self.solveElement(element)
            self.printGlobalMatrixH()
            #--------------------------------
            self.globalMatrixC = np.dot(self.globalMatrixC, (1 / self.grid.initData.simulationTimeStep))
            self.globalMatrixH = np.add(self.globalMatrixH, self.globalMatrixC)
            t0 = self.grid.getTempMesh()
            tmp = np.dot(self.globalMatrixC, t0)
            self.globalVectorP = np.add(tmp, self.globalVectorP)
            self.updateTemperature()
            tmp = np.array([(i+1) * self.grid.initData.simulationTimeStep, self.grid.getMinTemp(), self.grid.getMaxTemp()])
            self.MatrixMinMaxTemp.append(tmp)
            self.printEverythingAndSaveFile(i)
            self.printPlot((i+1) * self.grid.initData.simulationTimeStep)
        print("Koniec")
    def solveElement(self, element):
        jacobian = self.solveJacobian(element)
        #print("jacobian_element\n", jacobian)
        dN_dx = np.zeros((4, 4))
        dN_dy = np.zeros((4, 4))
        for i in range(self.ipa):
            tmp1_DetJ = (1 / self.detJ(jacobian[i]))
            for j in range(self.ipa):
                dN_dx[i, j] = tmp1_DetJ * (jacobian[i, 1, 1] * self.universalElement.dN_dksi[i,j] - jacobian[i, 0, 1] * self.universalElement.dN_deta[i,j])
                dN_dy[i, j] = tmp1_DetJ * (jacobian[i, 1, 0] * self.universalElement.dN_dksi[i,j] - jacobian[i, 0, 0] * self.universalElement.dN_deta[i,j])
        #rozwiazanie calki integral(k * ({N/dx}*{N/dx}^T} + {N/dy}*{N/dy}^T) dV) 6.9 FEM_transient_2d
        localMatrixHfirst = np.zeros((4, 4))
        localMatrixHsecond = np.zeros((4, 4))
        localMatrixC = np.zeros((4, 4))
        localVectorP = np.zeros((4))
        for i in range(self.ipa):
            tmpH = np.zeros((4, 4))
            tmpC = np.zeros((4, 4))
            detJ = self.detJ(jacobian[i])
            #-----------------------------MatrixH-----------------------------
            tmp_dN_dx = np.outer(dN_dx[i], dN_dx[i])
            tmp_dN_dy = np.outer(dN_dy[i], dN_dy[i])
            tmpH = np.add(tmpH, tmp_dN_dx)
            tmpH = np.add(tmpH, tmp_dN_dy)
            tmpH = np.dot(tmpH, detJ)
            localMatrixHfirst = np.add(localMatrixHfirst, tmpH)
            #-----------------------------matrixC-----------------------------
            tmpN_x_N = np.outer(self.universalElement.N[i], self.universalElement.N[i])
            tmpC = np.add(tmpC, tmpN_x_N)
            tmpC = np.dot(tmpC, detJ)
            localMatrixC = np.add(localMatrixC, tmpC)

        tmp_specificHeat_x_density = self.grid.initData.specificHeat * self.grid.initData.density # *c*ro
        localMatrixC = np.dot(localMatrixC, tmp_specificHeat_x_density)
        localMatrixHfirst = np.dot(localMatrixHfirst, self.grid.initData.conductivity)
        #-----------------------------Boundary Condidions:-----------------------------
        edgesN = self.universalElement.edgesN
        liczba = 0
        for i in range(self.ipa):
            id1 = element.ID[i]
            k = i + 1
            if(i == 3):
                k = 0
            id2 = element.ID[k]
            node1 = self.grid.nodes[id1]
            node2 = self.grid.nodes[id2]
            if(self.checkBoundaryCondition(node1, node2)):
                for j in range(2):
                    tmp1 = np.outer(edgesN[i][j], edgesN[i][j])
                    tmpH = np.dot(tmp1, self.det1J(node1, node2))
                    localMatrixHsecond = np.add(localMatrixHsecond, tmpH)
                    # -----------------------------vectorP-----------------------------
                    tmpP = np.dot(edgesN[i][j], self.det1J(node1, node2))
                    localVectorP = np.add(localVectorP, tmpP)
        tmp_alfa_ambTemp = (self.grid.initData.alfa * self.grid.initData.ambientTemperature)
        localVectorP = np.dot(localVectorP, tmp_alfa_ambTemp)
        localMatrixHsecond = np.dot(localMatrixHsecond, self.grid.initData.alfa)
        localMatrixH = np.add(localMatrixHfirst, localMatrixHsecond)
        self.mergeToGlobal(localMatrixH, localMatrixC, localVectorP, element)
    def updateTemperature(self):
        vectorTemp = []
        for node in self.grid.nodes:
            vectorTemp.append(node.temp)
        vectorTemp = np.asarray(vectorTemp).reshape(self.grid.initData.nH * self.grid.initData.nL, 1)
        A = np.linalg.inv(self.globalMatrixH)
        vectorTempSolved = A.dot(self.globalVectorP)
        for i in range(len(self.grid.nodes)):
            self.grid.nodes[i].temp = vectorTempSolved[i]
    def checkBoundaryCondition(self, node1, node2):
        #print(node1.border and node2.border)
        return (node1.border and node2.border)
    def detJ(self, jacobian):
        return jacobian[0, 0] * jacobian[1, 1] - jacobian[1, 0] * jacobian[0, 1]
    def det1J(self, node1, node2):
        result = math.pow((node1.x - node2.x), 2) + math.pow((node1.y - node2.y), 2)
        result = math.sqrt(result) / 2
        return result
    def solveJacobian(self, element):
        jacobian = np.zeros((self.ipa, 2, 2))
        for i in range(4):
            for j in range(4):
                id = element.ID[j] # Zamienieno kolejnosc
                jacobian[i, 1, 0] += (self.universalElement.dN_deta[i, j] * self.grid.nodes[id].x) #dx_deta
                jacobian[i, 1, 1] += (self.universalElement.dN_deta[i, j] * self.grid.nodes[id].y) #dy_deta
                jacobian[i, 0, 0] += (self.universalElement.dN_dksi[i, j] * self.grid.nodes[id].x) #dx_dksi
                jacobian[i, 0, 1] += (self.universalElement.dN_dksi[i, j] * self.grid.nodes[id].y) #dy_dksi
        #print(element)
        #print(jacobian[0])
        return jacobian
    def mergeToGlobal(self, matrixH, matrixC, vectorP, element):
        for i in range(len(matrixH)):
            for j in range(len(matrixH[i])):
                self.globalMatrixH[element.ID[i], element.ID[j]] += matrixH[i, j]
                self.globalMatrixC[element.ID[i], element.ID[j]] += matrixC[i, j]
            self.globalVectorP[element.ID[i]] += vectorP[i]
    def clearGlobal(self):
        self.globalVectorP = np.zeros(len(self.grid.nodes))
        self.globalMatrixH = np.zeros((len(self.grid.nodes),len(self.grid.nodes)))
        self.globalMatrixC = np.zeros((len(self.grid.nodes),len(self.grid.nodes)))
    # ----------------------------------PRINTING--------------------------------------
    def printGlobalVectorP(self):
        print("{P} + {[C]/dt}*{T0}:")
        print(self.globalVectorP)
        print()
    def printGlobalMatrixH(self):
        print("Matrix H:")
        for row in self.globalMatrixH:
            print(row)
        print()
    def printGlobalMatrixH_Everything(self):
        print("-------------------------------------------------------------------------------------------------------")
        print("[H] + [C]/dt:")
        for row in self.globalMatrixH:
            print(row)
        print()
    def printGlobalMatrixC(self):
        print("Matrix C:")
        for row in self.globalMatrixC:
            print(row)
        print()
    def printTemp(self):
        print("Min temp:" , self.grid.getMinTemp())
        print("Max temp:", self.grid.getMaxTemp())
    def printPlot(self, sec):
        x = np.linspace(0, self.grid.initData.L, self.grid.initData.nL)
        y = np.linspace(0, self.grid.initData.H, self.grid.initData.nH)
        X, Y = np.meshgrid(x, y)
        z = np.asarray(self.grid.getTempMesh())
        Z = np.zeros((self.grid.initData.nH,self.grid.initData.nL))
        index = 0
        for i in range(self.grid.initData.nH):
            for j in range(self.grid.initData.nL):
                Z[i, j ] = z[index]
                index = index + 1
        plt.contourf(X, Y, Z, 20, cmap='jet');
        plt.title(str(int(sec)) + " sekund")
        cb = plt.colorbar()
        plt.savefig('plot' + str(sec) + '.png')
        plt.show()
    def toFile(self, iteration):
        path = os.path.dirname(__file__)
        path = os.path.join(path, 'output')
        df = pd.DataFrame(self.globalMatrixH)
        df.to_csv(path + "\\matrixH" + str(iteration) + ".csv")
        df = pd.DataFrame(self.globalVectorP)
        df.to_csv(path + "\\vectorP" + str(iteration) + ".csv")
        df = pd.DataFrame(self.globalMatrixC)
        df.to_csv(path + "\\matrixC" + str(iteration) + ".csv")
        df = pd.DataFrame(self.MatrixMinMaxTemp)
        df.to_csv(path + "\\MatrixMinMaxTemp.csv")
    def printEverythingAndSaveFile(self, iteration):
        self.printGlobalMatrixH_Everything()
        self.printGlobalVectorP()
        self.printTemp()
        self.toFile(iteration)
    def setSubPlot(self, iteration):
        seconds = (iteration + 1) * self.grid.initData.simulationTimeStep
        #self.axs[iteration, 0) = plt.contourf(X, Y, Z, 20, cmap='jet');
        pass


