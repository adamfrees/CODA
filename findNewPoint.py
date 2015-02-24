from numpy import *
from sklearn import linear_model
import csv
import sys
set_printoptions(precision=5,linewidth=1000,suppress=True)
from regularize import regularize
from numpy.linalg import norm


if __name__ == '__main__':
    filename = str(sys.argv[1])
    acceptanceError = str(sys.argv[2])
    try:
        acceptanceError = float(acceptanceError)
    except ValueError:
        acceptanceError = 0.5
    target = [float(x) for x in str(sys.argv[2]).split(',')]
    data = regularize(filename,15,[4,5,6])
    data.makeRelativeData()
    data.declareTarget(target)
    #data.scaleTargets([4,5,6],10.0)
    data.useTargetIndices([0,1,4,5,6])
    data.constructModel()
    data.compressedSense()
    dataListOccupation = data.createNormList()
    includeInError=target[0:4] 
    for i,x in enumerate(data.results):
        #print x
        if x[1]>acceptanceError*norm(includeInError):
        	index = i
        	break
    voltages = map(lambda x: x*(abs(x)>1e-4),data.correspondingVecs)
    toPrint = ""
    
    #for i,x in enumerate(voltages):
    #    if sum(x)!=0.:
    #        index = i
    #        break
    #index = index #+ (len(voltages)-index)/2
    for x in data.controlVecs[index].tolist():
        toPrint += "%.4f" % x + ","
    print toPrint[0:-1]
