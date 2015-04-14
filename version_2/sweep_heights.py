from numpy import *
import pymatlab
from model import model
from ReadOutput import processData
from regularize import regularize
from numpy.linalg import norm
from operator import add
from CODA import *
import os


reloadData = True

dataFolder = '/home/frees/Dropbox/_UW/CODA/Var_Height/Data'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height.mph'

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

target = [1.,1.,0,0,0.01,0.01,0.01]

maxHeight=0.09

height=0.08

heightNum = 0
sweepNum = 0
runNum = 0

dDot = model()
dDot.loadModelFile(modelFile)

if reloadData:
    f = open(dataFolder+'/heightList.txt','r')
    heightData = f.readlines()[-1]
    f.close()
    heightNum = int(heightData[0:heightData.find(':')])
    height = float(heightData[heightData.find(':')+1:-1])
    sweepList = []
    runList = []
    for file in os.listdir(dataFolder):
        if file.find('height_'+str(heightNum))!=-1 and file.find('sweep')!=-1:
            sweepList += [file]
    sweepNum = len(sweepList)-1
    f = open(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt','r')
    voltData = f.readline()
    f.close()
    voltages = [float(x) for x in voltData.split(' ')[0:len(order)]]
    for gateNum,gate in enumerate(order):
        dDot.setVoltage('pot'+str(gate),voltages[gateNum])
    runNum = (len(order)-1)*(sweepNum-1)
    os.remove(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt')
    dDot.changeHeight('wp4',height)


iterateHeight = True

while iterateHeight:
    #If not the first height, find optimal height
    if heightNum != 0 and not reloadData:
        heightStep=0.001
        lookingForNewHeight = True
        numAttemptsNewHeight =0
        while lookingForNewHeight:
            numAttemptsNewHeight +=1
            height += heightStep
            dDot.changeHeight('wp4',height)
            try:
                dDot.runSimulation()
                updateFiles(dDot, voltages, heightNum, runNum, sweepNum)
                lookingForNewHeight = False
            except TransmissionCoefficientError:
                height -=heightStep
                heightStep = heightStep*0.8
            if numAttemptsNewHeight>10:
                raise TransmissionCoefficientError('Transmission coefficient failed to converge - failure in finding new height.')
        sweepNum = 0
        runNum = 0

    toWrite = str(heightNum)+': '+str(height)+'\n'        
    f = open(dataFolder+'/heightList.txt','a')
    f.write(toWrite)
    f.close()
    
    iterateSweep = True
    
    while iterateSweep:
        #If not the first sweep of the height, find new operating point.
        if sweepNum != 0 and not reloadData:
            acceptableL2Error = 0.5
            lookingForNewSweepOperatingPoint = True
            numAttemptsNewSweepOperatingPoint = 0
            while lookingForNewSweepOperatingPoint:
                numAttemptsNewSweepOperatingPoint +=1
                newStep = findNewPoint(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum-1)+'.txt',target,acceptanceError=acceptableL2Error)
                voltages = map(add,voltages,newStep)
                for gateNum,gate in enumerate(order):
                    dDot.setVoltage('pot'+str(gate),voltages[gateNum])
                try:
                    dDot.runSimulation()
                    updateFiles(dDot, voltages, heightNum, runNum, sweepNum)
                    lookingForNewSweepOperatingPoint = False
                except TransmissionCoefficientError:
                    acceptableL2Error = 1.- 0.8*(1.-acceptableL2Error) #Bringing the acceptable error closer to 1 limits the L1 Norm of the solution, making the step size smaller.
                    voltages = map(lambda x,y: x-y, voltages,newStep)
                if numAttemptsNewSweepOperatingPoint>10:
                    raise TransmissionCoefficientError('Transmission coefficient failed to converge - failure in finding new sweep Operating point.')
            
        #Iterate through runs (perturb voltage on each gate in turn).
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
        
        reloadData = False
        if isConverged(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt',target) or sweepNum>5:
            iterateSweep = False
        sweepNum += 1
    heightNum +=1
    if height>maxHeight:
        iterateHeight = False