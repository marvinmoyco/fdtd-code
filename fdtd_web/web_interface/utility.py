from . import simulation as fdtd

import numpy as np


def loadJSONdata(json_data):
    
    # Create a new simulation object and return it
    sim = fdtd.Simulation(json_data=json_data)

    

    return sim


def parseSimParam(data):
    ret = []
    for x in range(len(data)):
        ret.append(float(data[x]))

    return np.array(ret)
        

    


def createCSVfile(json_data):
    pass