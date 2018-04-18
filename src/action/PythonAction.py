mySourceData = 0
myTargetData = 0

def performAction(time, dt, sourceData, targetData):
    ''' This function is called first at configured timing. It can be omitted, if not 
    needed. Its parameters are time, (global) timestep size, the source data, followed by the target data. 
    Source and target data can be omitted (selectively or both) by not mentioning 
    them in the preCICE XML configuration (see the configuration reference).'''

    # Usage example 1:
    global mySourceData
    global myTargetData
    mySourceData = sourceData # store (reference to) sourceData for later use
    myTargetData = targetData # store (reference to) targetData for later use
    # Usage example 2:
    # for i in range(data.size):
    #     data[i] = data[i] + 1 # Add 1 to each data component
    #     i = i + 1
    
def vertexCallback(id, coords, normal):
    '''This function is called for every vertex in the configured mesh. It is called
    after performAction, and can also be omitted.'''

    # Usage example:
    global mySourceData # Make global data set in performAction visible
    global myTargetData
    # myTargetData[id] += coords[0] + mySourceData[id] # Add data to vertex coords
    
def postAction():
    '''This function is called at last, if not omitted.'''
    
    global mySourceData # Make global data set in performAction visible
    global myTargetData
    # Do something ...
    
