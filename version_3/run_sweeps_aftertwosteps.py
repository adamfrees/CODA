from numpy import *
from model import model
from time import time
import os

dataFolder = '/home/frees/Dropbox/_UW/CODA/version3/Data'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height_3.mph'

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

params = [
    [0.03,[0.4,0.45408,0.4,-0.4,0.0722,-0.175,-0.01096,-0.00903,-0.11769,0.07804,-0.4,-0.2,-0.25,-0.264,-0.2]],
    [0.032,[0.4,0.45,0.4,-0.4,0.07567,-0.19354,-0.01226,-0.012,-0.1222,0.08368,-0.4,-0.2,-0.26,-0.26,-0.2]],
    [0.034,[0.4,0.45,0.4,-0.4,0.07671,-0.16176,-0.01628,-0.012,-0.12483,0.0894,-0.4,-0.2,-0.27,-0.256,-0.2]],
    [0.036,[0.4,0.45,0.4,-0.4,0.07801,-0.1573,-0.01934,-0.012,-0.11582,0.0927,-0.4,-0.2,-0.285,-0.252,-0.2]],
    [0.038,[0.4,0.45,0.4,-0.4,0.07939,-0.15473,-0.02061,-0.012,-0.1591,0.10481,-0.4,-0.2,-0.295,-0.248,-0.2]],
    [0.04,[0.4,0.45,0.4,-0.4,0.08271,-0.16108,-0.02436,-0.012,-0.1268,0.10201,-0.4,-0.2,-0.305,-0.24,-0.2]]
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
    
    voltageStepSize = 0.001
    
    paramMat += voltageStepSize*eye(len(order))
    
    dDot.initializeSweep(gateList,paramMat)
    
    startTime = time()
    dDot.runSimulation('std2')
    
    for i in range(len(order)):
        dDot.exportData(dataFolder+'/Height_'+str(height)+'/Run_'+str(runNum)+'.txt','data3',i+1)
        runNum +=1
    
    print "Subsequent 15 runs: ",time()-startTime
