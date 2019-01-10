import numpy as np
from skrypt.Element import Element
from skrypt.Node import Node
from skrypt.InitialData import InitialData

class Grid:
    def __init__(self,filename):
        self.initData = InitialData(filename)
        self.nodes = []
        self.elements = []
        self.genGrid()
    def genGrid(self):
        #Nodes:
        dH = self.initData.H / (self.initData.nH - 1)
        dL = self.initData.L / (self.initData.nL - 1)
        x = 0
        y = 0
        for i in range(self.initData.nL):
            y = 0
            for j in range(self.initData.nH):
                border = self.isBorder(i,j)
                self.nodes.append(Node(x,y, self.initData.initialTemperature,border))
                y = y + dH
            x = x + dL
        #Elements:
        index = 0
        for i in range(self.initData.nL - 1):
            for j in range(self.initData.nH - 1):
                a = index + i
                b = a + self.initData.nH
                c = b + 1
                d = a + 1
                element = Element(a,b,c,d)
                self.elements.append(element)
                index = index + 1
    def isBorder(self, x, y):
        if x == (self.initData.nL - 1) or y == (self.initData.nH - 1) or x == 0 or y == 0:
            return 1
        return 0
    #------------------------------------GET-----------------------------------------
    def getMinTemp(self):
        min = self.nodes[0].temp
        for node in self.nodes:
            if(min > node.temp):
                min = node.temp
        return min
    def getMaxTemp(self):
        max = self.nodes[0].temp
        for node in self.nodes:
            if(max < node.temp):
                max = node.temp
        return max
    def getTempMesh(self):
        tempMesh = []
        for node in self.nodes:
            tempMesh.append(node.temp)
        return tempMesh
    #----------------------------------PRINTING--------------------------------------
    def printNodes(self):
        for node in self.nodes:
            print(node)
    def printElements(self):
        for element in self.elements:
            print(element)
        
