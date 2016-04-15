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

target = [1,1,1,1,1,1,1,1,0.01,0.01,0.01,0.01]

sweepNum = 0
runNum = 0

octoDot = model()
octoDot.loadModelFile(modelFile)

voltages = zeros(len(order))
for gateNum,gate in enumerate(order):
    voltages[gateNum] = octoDot.getVoltage('pot'+str(gate))

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
