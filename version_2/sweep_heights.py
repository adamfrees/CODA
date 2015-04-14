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
numHeights = 16
#heightList = linspace(0.08,0.09,numHeights)
height=0.08
maxHeight=0.09
heightStep=0.001
heightNum=0

dDot = model()
dDot.loadModelFile(modelFile)
oldStep = []

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
    voltages = voltData.split(' ')[0:len(order)]
    for gateNum,gate in enumerate(order):
        dDot.setVoltage('pot'+str(gate),voltages[gateNum])
    runNum = (len(order)-1)*(sweepNum-1)
    os.remove(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt')

while height<maxHeight:
    dDot.changeHeight('wp4',height)
    looping =True
    if not reloadData:
        runNum = 0
        sweepNum = 0
    reloadData = False
    while looping:
        voltages = []
        for gate in order:
            voltages += [dDot.getVoltage('pot'+str(gate))]
        numTries = 1
        maxTries = 10
        lookingForNextStep = True
        needNewHeight=True
        while lookingForNextStep:
            acceptanceError = 0.5
            if numTries>maxTries:
                raise TransmissionCoefficientError('Transmission coefficient failed to converge - failure in finding new operating point.')
            numTries +=1
            try:
                dDot.runSimulation()
                updateFiles(dDot, voltages, heightNum, runNum, sweepNum)
                lookingForNextStep = False
            except TransmissionCoefficientError:
                if sweepNum==0:
                    break
                acceptanceError = 1.- 0.8*(1.-acceptanceError) #Bringing the acceptable error closer to 1 limits the L1 Norm of the solution, making the step size smaller.
                voltages = map(lambda x,y: x-y, voltages,oldStep)
                converged,newStep = findNewPoint(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum-1)+'.txt',target,acceptanceError=acceptanceError)
                voltages = map(add,voltages,newStep)
                oldStep = newStep
        else:
            needNewHeight=False
        if needNewHeight:#this runs if there is a TransmissionCoefficientError on the first sweep on a height
            if heightNum==0:
                raise TransmissionCoefficientError('Transmission coefficient failed to converge - bad starting operating point.')
            else:
                height -= heightStep
                heightStep = heightStep*0.6
                height += heightStep
                break
        toWrite = str(heightNum)+': '+str(height)+'\n'        
        f = open(dataFolder+'/heightList.txt','a')
        f.write(toWrite)
        f.close()    
                 
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
        
        converged,newStep = findNewPoint(dataFolder+'/height_'+str(heightNum)+'_sweep_'+str(sweepNum)+'.txt',target)
        if converged or sweepNum>5:
            looping = False
            height +=heightStep
            heightNum += 1
        else:
            voltages = map(add,voltages,newStep)
            for gateNum, gate in enumerate(order):
                dDot.setVoltage('pot'+str(gate),voltages[gateNum])
            oldStep = newStep
            sweepNum+=1
