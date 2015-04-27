mySourceData = 0
myTargetData = 0

#
# This function is called first. It can be omitted, if not needed. Its
# parameters are the source data, followed by the target data, which are
# omitted (selectively or both) if not mentioned in the preCICE configuration.
#
def performAction(time, sourceData, targetData):
    global mySourceData
    global myTargetData
    mySourceData = sourceData # store (reference to) sourceData for later use
    myTargetData = targetData # store (reference to) targetData for later use
    # Usage example:
    for i in range(myTargetData.size):
        myTargetData[i] = mySourceData[i] - 1
   
    
