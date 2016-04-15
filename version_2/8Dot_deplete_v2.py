from numpy import *
import pymatlab
from model import model
from ReadOutput import processData
from regularize import regularize
from numpy.linalg import norm
from operator import add
from CODA import *
import os
from time import time

heightNum=0
reloadData = False

dataFolder = '/home/frees/Dropbox/_UW/CODA/8Dot/Data'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/8Dots/COMSOL/8Dots_v2.mph'

order = range(19,60)

target = [1.,1.,1.,1.,1.,1.,1.,1.,0.01,0.01,0.01,0.01]

sweepNum = 0
runNum = 0

octoDot = model()
octoDot.loadModelFile(modelFile)

voltages = array([0.4,0.4,0.4,0.4,-0.0595889752117,-0.02,0.00105537163149,-0.02,-0.03,0.00190583399534,-0.0195348709,-0.0182482415134,0.00100138261228,-0.02,-0.02,-0.00235333612439,-0.02,-0.06,-0.1,0.04,-0.1,-0.06,0.0345949672212,-0.05,0.0242011620479,-0.1,0.024055166628,-0.0484,0.0240067069583,-0.1,0.0242025279056,-0.05,0.0252998181465,-0.1,0.0244078867626,-0.05,0.0225339476475,-0.06,-0.1,0.04,-0.1])
for gateNum,gate in enumerate(order):
    octoDot.setVoltage('pot'+str(gate),voltages[gateNum])

iterateSweep = True

while iterateSweep:
    #If not the first sweep of the height, find new operating point.
    if sweepNum != 0 and not reloadData:
        acceptableL2Error = 0.5
        lookingForNewSweepOperatingPoint = True
        numAttemptsNewSweepOperatingPoint = 0
        while lookingForNewSweepOperatingPoint:
            numAttemptsNewSweepOperatingPoint +=1
            newStep = findNewPoint8Dot(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum-1)+'.txt',target,acceptanceError=acceptableL2Error)
            voltages = map(add,voltages,newStep)
            for gateNum,gate in enumerate(order):
                octoDot.setVoltage('pot'+str(gate),voltages[gateNum])
            try:
                octoDot.runSimulation()
                updateFiles(octoDot, voltages, heightNum, runNum, sweepNum, modelType=2)
                lookingForNewSweepOperatingPoint = False
            except TransmissionCoefficientError:
                acceptableL2Error = 1.- 0.8*(1.-acceptableL2Error) #Bringing the acceptable error closer to 1 limits the L1 Norm of the solution, making the step size smaller.
                voltages = map(lambda x,y: x-y, voltages,newStep)
            if numAttemptsNewSweepOperatingPoint>10:
                raise TransmissionCoefficientError('Transmission coefficient failed to converge - failure in finding new sweep Operating point.')
        
    #Iterate through runs (perturb voltage on each gate in turn).
    for gateNum, gate in enumerate(order):
        voltChange = 0.0008
        numTries = 0
        maxTries = 10
        while True:
            if numTries>maxTries:
                raise TransmissionCoefficientError('Transmission coefficient failed to converge - failure in modulating gate voltage.')
            numTries +=1
            voltages[gateNum] += voltChange
            octoDot.setVoltage('pot'+str(gate),voltages[gateNum])
            startTime = time()
            octoDot.runSimulation()
            print time()-startTime
            runNum += 1
            try:
                updateFiles(octoDot, voltages, heightNum, runNum, sweepNum, modelType=2)
                voltages[gateNum] -= voltChange
                octoDot.setVoltage('pot'+str(gate),voltages[gateNum])
                break
            except TransmissionCoefficientError:
                voltages[gateNum] -= voltChange
                voltChange = voltChange/2.
    
    reloadData = False
    if isConverged8Dot(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt',target) or sweepNum>5:
        iterateSweep = False
    sweepNum += 1
