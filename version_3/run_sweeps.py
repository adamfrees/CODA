from numpy import *
from model import model
from time import time
import os

dataFolder = '/home/frees/Dropbox/_UW/CODA/version3/Data/0_005'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height_3.mph'

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

params = [
    [0.03,[0.4,0.45,0.4,-0.4,0.0722,-0.175,-0.0087,-0.012,-0.094,0.0642,-0.4,-0.2,-0.25,-0.264,-0.2]],
    [0.032,[0.4,0.45,0.4,-0.4,0.075,-0.18,-0.01,-0.012,-0.1,0.0685,-0.4,-0.2,-0.26,-0.26,-0.2]],
    [0.034,[0.4,0.45,0.4,-0.4,0.077,-0.177,-0.012,-0.012,-0.105,0.073,-0.4,-0.2,-0.27,-0.256,-0.2]],
    [0.036,[0.4,0.45,0.4,-0.4,0.079,-0.1765,-0.014,-0.012,-0.109,0.0775,-0.4,-0.2,-0.285,-0.252,-0.2]],
    [0.038,[0.4,0.45,0.4,-0.4,0.082,-0.179,-0.016,-0.012,-0.114,0.081,-0.4,-0.2,-0.295,-0.248,-0.2]],
    [0.04,[0.4,0.45,0.4,-0.4,0.084,-0.183,-0.0177,-0.012,-0.12,0.0855,-0.4,-0.2,-0.305,-0.24,-0.2]]
    ]

dDot = model()
dDot.loadModelFile(modelFile)

for heightParams in params:
    height = heightParams[0]
    startingVolts = heightParams[1]
    dDot.setParam('SiGe',height)
    for i,gateNum in enumerate(order):
        dDot.setParam('pot'+str(gateNum),startingVolts[i])
    
    runNum = 0
    
    startTime = time()
    dDot.runSimulation('std1')
    os.makedirs(dataFolder+'/Height_'+str(height))
    dDot.exportData(dataFolder+'/Height_'+str(height)+'/Run_'+str(runNum)+'.txt','data1')
    runNum +=1
    print "First run: ", time()-startTime
    
    gateList = []
    paramList = []
    for gateNum in order:
        paramList += [dDot.getParam('pot'+str(gateNum))]
        gateList += ['pot'+str(gateNum)]
    paramMat = array([paramList]*len(order)).transpose()
    
    voltageStepSize = 0.005
    
    paramMat += voltageStepSize*eye(len(order))
    
    dDot.initializeSweep(gateList,paramMat)
    
    startTime = time()
    dDot.runSimulation('std2')
    
    for i in range(len(order)):
        dDot.exportData(dataFolder+'/Height_'+str(height)+'/Run_'+str(runNum)+'.txt','data3',i+1)
        runNum +=1
    
    print "Subsequent 15 runs: ",time()-startTime
