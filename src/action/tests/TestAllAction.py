mySourceData = 0
myTargetData = 0
myIteration = 0

#
# This function is called first. It can be omitted, if not needed. Its
# parameters are the source data, followed by the target data, which are
# omitted (selectively or both) if not mentioned in the preCICE configuration.
#
def performAction(time, dt, sourceData, targetData):
    global mySourceData
    global myTargetData
    mySourceData = sourceData # store (reference to) sourceData for later use
    myTargetData = targetData # store (reference to) targetData for later use
    # Usage example:
    for i in range(myTargetData.size):
        myTargetData[i] = mySourceData[i] + 1
    
#
# This function is called for every vertex in the configured mesh. It is called
# after performAction, and can also be omitted.
#
def vertexCallback(id, coords, normal):
    global mySourceData
    global myTargetData
    # Usage example:
    myTargetData[id] += coords[0]
    
#
# This function is called at last, if not omitted.
#
def postAction():
    global mySourceData
    global myTargetData
    global myIteration
    for i in range(myTargetData.size):
        myTargetData[i] -= myIteration
    myIteration += 1
    
