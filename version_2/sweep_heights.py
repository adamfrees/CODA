from numpy import *
import pymatlab
from model import model
from ReadOutput import processData
from regularize import regularize
from numpy.linalg import norm
from operator import add

class TransmissionCoefficientError(Exception):
     pass

def updateFiles(model, voltages, heightNum, runNum, sweepNum):
    model.exportData(dataFolder+'/height_'+str(heightNum)+'_run_'+str(runNum)+'.txt')
    procDat = processData(dataFolder+'/height_'+str(heightNum)+'_run_'+str(runNum)+'.txt',modelType=0)
    procDat.processPotential(None,None)
    for transmission in procDat.Ts:
        if transmission <=0. or transmission >= 1.:
            print "Error:", voltages
            print procDat.density, procDat.Ts
            raise TransmissionCoefficientError('Transmission coefficient outside acceptable range.')
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
    return data.controlVecs[index]

dataFolder = '/home/frees/Dropbox/_UW/CODA/Var_Height/Data'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height.mph'

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

target = [1.,1.,0,0,0.01,0.01,0.01]
numHeights = 11
heightList = linspace(0.08,0.09,numHeights)

dDot = model()
dDot.loadModelFile(modelFile)
oldStep = []
for heightNum,height in enumerate(heightList):
    dDot.changeHeight('wp4',height)
    #looping =True
    runNum = 0
    sweepNum = 0
    for myCount in range(5):
        voltages = []
        for gate in order:
            voltages += [dDot.getVoltage('pot'+str(gate))]
        numTries = 0
        maxTries = 10
        while True:
            acceptanceError = 0.5
            if numTries>maxTries:
                raise TransmissionCoefficientError('Transmission coefficient failed to converge - failure in finding new operating point.')
            numTries +=1
            try:
                dDot.runSimulation()
                updateFiles(dDot, voltages, heightNum, runNum, sweepNum)
                break
            except TransmissionCoefficientError:
                if sweepNum==0:
                    raise TransmissionCoefficientError('Transmission coefficient failed to converge - bad starting operating point.')
                acceptanceError = 1.- 0.8*(1.-acceptanceError) #Bringing the acceptable error closer to 1 limits the L1 Norm of the solution, making the step size smaller.
                voltages = map(lambda x,y: x-y, voltages,oldStep)
                newStep = findNewPoint(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum-1)+'.txt',target,acceptanceError=acceptanceError)
                voltages = map(add,voltages,newStep)
                oldStep = newStep
                 
        for gateNum, gate in enumerate(order):
            voltChange = -0.0016
            numTries = 0
            maxTries = 10
            while True:
                if numTries>maxTries:
                    raise TransmissionCoefficientError('Transmission coefficient failed to converge - failure in modulating gate voltage.')
                numTries +=1
                voltages[gateNum] += voltChange
                dDot.setVoltage('pot'+str(gate),voltages[gateNum])
                dDot.runSimulation()
                runNum += 1
                try:
                    updateFiles(dDot, voltages, heightNum, runNum, sweepNum)
                    voltages[gateNum] -= voltChange
                    dDot.setVoltage('pot'+str(gate),voltages[gateNum])
                    break
                except TransmissionCoefficientError:
                    voltages[gateNum] -= voltChange
                    voltChange = voltChange/2.
        
        newStep = findNewPoint(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt',target)
        voltages = map(add,voltages,newStep)
        for gateNum, gate in enumerate(order):
            dDot.setVoltage('pot'+str(gate),voltages[gateNum])
        oldStep = newStep
        sweepNum+=1
