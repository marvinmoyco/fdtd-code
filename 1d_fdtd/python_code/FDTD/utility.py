import simulation as fdtd




def loadJSONdata(json_data):
    
    # Create a new simulation object and return it
    sim = fdtd.Simulation(json_data=json_data)


    return sim


def createCSVfile(json_data=json_data):
    pass