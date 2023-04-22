def performAction(time, sourceData, targetData):
    ''' This function is called at the configured timing.
    Its parameters are time, the source data, followed by the target data.
    Source and target data can be omitted (selectively or both) by not mentioning
    them in the preCICE XML configuration (see the configuration reference).'''

    # Usage example:
    # for i in range(sourceData.size):
    #     targetData[i] = sourceData[i] + 1 # Add 1 to each data component
    #     i = i + 1
