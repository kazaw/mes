import json
import io
import sys


def readJSON(filename):
    with open(filename, "r") as read_file:
        data = json.load(read_file)
        return data
class InitialData:
    def __init__(self,filename):
        dane = readJSON(filename)
        self.initialTemperature = float(dane["initialTemperature"])
        self.simulationTime = float(dane["simulationTime"])
        self.simulationTimeStep = float(dane["simulationTimeStep"])
        self.ambientTemperature = float(dane["ambientTemperature"])
        self.alfa = float(dane["alfa"])
        self.H = float(dane["H"])
        self.L = float(dane["L"])
        self.nH = int(dane["nH"])
        self.nL = int(dane["nL"])
        self.specificHeat = float(dane["specificHeat"])
        self.conductivity = float(dane["conductivity"])
        self.density = float(dane["density"])





