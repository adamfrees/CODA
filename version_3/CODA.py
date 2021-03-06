from numpy import *
import pymatlab
from model import model
from ReadOutput import processData
#from regularize import regularize
from numpy.linalg import norm
from operator import add

#dataFolder = '/home/frees/Dropbox/_UW/CODA/Var_Height/Data'
dataFolder = '/home/frees/Dropbox/_UW/CODA/8Dot/Data'

class TransmissionCoefficientError(Exception):
     pass

def analyzeFile(filename,modelType=0):
    procDat = processData(filename,modelType=modelType)
    procDat.processPotential(None,None)
    for i,transmission in enumerate(procDat.Ts):
        if transmission <=1.e-8 or transmission >= 1.:
            raise TransmissionCoefficientError('Transmission coefficient outside acceptable range. Coupling '+str(i)+' is '+str(transmission))
    return procDat.densities, procDat.Ts


def updateFiles(model, voltages, heightNum, runNum, sweepNum,modelType=0):
    model.exportData(dataFolder+'/height_'+str(heightNum)+'_run_'+str(runNum)+'.txt')
    procDat = processData(dataFolder+'/height_'+str(heightNum)+'_run_'+str(runNum)+'.txt',modelType=modelType)
    procDat.processPotential(None,None)
    for i,transmission in enumerate(procDat.Ts):
        if transmission <=1.e-8 or transmission >= 1.:
            raise TransmissionCoefficientError('Transmission coefficient outside acceptable range. Coupling '+str(i)+' is '+str(transmission))
    toWrite = ""
    for volts in voltages:
        toWrite += str(volts)+" "
    for density in procDat.densities:
        toWrite += str(density)+" "
    for transmission in procDat.Ts:
        toWrite += str(transmission)+" "
    toWrite = toWrite+'\n'
    f = open(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt','a')
    f.write(toWrite)
    f.close()

def findNewPoint(filename,target,acceptanceError=0.5):
    data = regularize(filename,15,[4,5,6])
    data.makeRelativeData()
    data.declareTarget(target)
    #data.scaleTargets([4,5,6],10.0)
    data.useTargetIndices([0,1,4,5,6])
    data.constructModel()
    data.compressedSense()
    dataListOccupation = data.createNormList() 
    for i,x in enumerate(data.results):
        #print x
        if x[1]>acceptanceError*norm(data.target):
        	index = i
        	break
    #voltages = map(lambda x: x*(abs(x)>1e-4),data.correspondingVecs)
    
    #for i,x in enumerate(voltages):
    #    if sum(x)!=0.:
    #        index = i
    #        break
    #index = index #+ (len(voltages)-index)/2
    return  data.controlVecs[index]

def findNewPoint8Dot(filename,target,acceptanceError=0.5):
    data = regularize(filename,41,[8,9,10,11])
    data.makeRelativeData()
    data.declareTarget(target)
    #data.scaleTargets([4,5,6],10.0)
    data.useTargetIndices([0,1,2,3,4,5,6,7,8,9,10,11])
    data.constructModel()
    data.compressedSense()
    dataListOccupation = data.createNormList() 
    for i,x in enumerate(data.results):
        #print x
        if x[1]>acceptanceError*norm(data.target):
        	index = i
        	break
    #voltages = map(lambda x: x*(abs(x)>1e-4),data.correspondingVecs)
    
    #for i,x in enumerate(voltages):
    #    if sum(x)!=0.:
    #        index = i
    #        break
    #index = index #+ (len(voltages)-index)/2
    return  data.controlVecs[index]

def isConverged(filename,target,convergenceTest=0.2):
    data = regularize(filename,15,[4,5,6])
    data.makeRelativeData()
    data.declareTarget(target)
    data.useTargetIndices([0,1,4,5,6])
    converged = False
    if norm(data.target)<convergenceTest:
        converged = True
    return converged

def isConverged8Dot(filename,target,convergenceTest=0.2):
    data = regularize(filename,41,[8,9,10,11])
    data.makeRelativeData()
    data.declareTarget(target)
    data.useTargetIndices([0,1,2,3,4,5,6,7,8,9,10,11])
    converged = False
    if norm(data.target)<convergenceTest:
        converged = True
    return converged
    
