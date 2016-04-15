from numpy import *
from regularize import regularize
import os
from analyze_data import analyzeDataFile

def l1ErrorPlot(filename,target):
    data = regularize(filename,15,[4,5,6])
    data.makeRelativeData()
    data.declareTarget(target)
    data.useTargetIndices([0,1,4,5,6])
    data.constructModel()
    data.compressedSense(alphaRange=[1.,100.],numSamples=100)
    dataListOccupation = data.createNormList()
    distance = sum(abs(data.target))
    toReturn = [[],[],[],[]]
    voltages = []
    for index in range(len(data.results)-2,1,-1):
        if data.results[index][1]>data.results[index+1][1]:
            data.results[index][1]=data.results[index+1][1]
            data.controlVecs[index] = data.controlVecs[index+1]
    for index,x in enumerate(data.results):
        if True:#sum(abs(data.controlVecs[index]))<0.004 and sum(abs(data.controlVecs[index]))>0.0001:
            toReturn[0]+= [x[2]]
            toReturn[1] += [x[1]]
            toReturn[2] += [x[1]/distance]
            toReturn[3] += [x[0]]
            voltages += [data.controlVecs[index]]
    return toReturn,voltages

codaDataFolder = '/home/frees/Dropbox/_UW/CODA/version3/Autotune/attempt_2/autoTuneIteration_4/CODA_data'

target = [1.,2.,0,0,0.01,0.01,0.01]

plot,voltageChanges = l1ErrorPlot(codaDataFolder+'/CODA_data.txt',target)
#print plot
for i in range(len(voltageChanges)):
    print sum(abs(voltageChanges[i])),plot[2][i]
