#
# Inverts direction of target data values 
#
def performAction(time, targetData):
    for i in range(targetData.size):
        targetData[i] = -targetData[i]