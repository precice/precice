myIteration = 0


#
# This function is called first. Its parameters are the source data, followed
# by the target data, which are omitted (selectively or both) if not mentioned
# in the preCICE configuration.
def performAction(time, sourceData, targetData):

    for i in range(targetData.size):
        targetData[i] = sourceData[i] + 1

    global myIteration
    for i in range(targetData.size):
        targetData[i] += myIteration
    myIteration += 1
