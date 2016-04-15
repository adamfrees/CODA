from numpy import *
from model import model
from time import time
from regularize import regularize
import os
from analyze_data import analyzeDataFile 

def l1ErrorPlot(filename,target):
    data = regularize(filename,15,[4,5,6])
    data.makeRelativeData()
    data.declareTarget(target)
    data.useTargetIndices([0,1,4,5,6])
    data.constructModel()
    data.compressedSense(alphaRange=[1.e-5,100000.],numSamples=5000)
    dataListOccupation = data.createNormList()
    distance = sum(abs(data.target)) 
    toReturn = [[],[],[]]
    voltages = []
    for index in range(len(data.results)-2,1,-1):
        if data.results[index][1]>data.results[index+1][1]:
            data.results[index][1]=data.results[index+1][1]
            data.controlVecs[index] = data.controlVecs[index+1]
    for index,x in enumerate(data.results):
        if sum(abs(data.controlVecs[index]))<0.004 and sum(abs(data.controlVecs[index]))>0.0001:
            toReturn[0]+= [x[2]]
            toReturn[1] += [x[1]]
            toReturn[2] += [x[1]/distance]
            voltages += [data.controlVecs[index]]
    return toReturn,voltages

def getColor(x):
    #outputs rgb color from rainbow, takes in x:0-1
    case = x*4.0
    if case <1.:
        return (1.,case,0.)
    if case <2.:
        return (2.-case,1.,0.)
    if case<3.:
        return (0.,1.,case-2.)
    return (0.,4-case,1.)

plotL1 = True

autoTuneRunning = True
autoTuneIteration = 0

parentFolder = '/home/frees/Dropbox/_UW/CODA/version3/Autotune/attempt_2'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height_3.mph'

allHeightsStartingParams = [
    [0.03,[0.4,0.45,0.4,-0.4,0.0722,-0.175,-0.0087,-0.012,-0.094,0.0642,-0.4,-0.2,-0.25,-0.264,-0.2]],
    [0.032,[0.4,0.45,0.4,-0.4,0.075,-0.18,-0.01,-0.012,-0.1,0.0685,-0.4,-0.2,-0.26,-0.26,-0.2]],
    [0.034,[0.4,0.45,0.4,-0.4,0.077,-0.177,-0.012,-0.012,-0.105,0.073,-0.4,-0.2,-0.27,-0.256,-0.2]],
    [0.036,[0.4,0.45,0.4,-0.4,0.079,-0.1765,-0.014,-0.012,-0.109,0.0775,-0.4,-0.2,-0.285,-0.252,-0.2]],
    [0.038,[0.4,0.45,0.4,-0.4,0.082,-0.179,-0.016,-0.012,-0.114,0.081,-0.4,-0.2,-0.295,-0.248,-0.2]],
    [0.04,[0.4,0.45,0.4,-0.4,0.084,-0.183,-0.0177,-0.012,-0.12,0.0855,-0.4,-0.2,-0.305,-0.24,-0.2]]
    ]

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

target = [1.,2.,0,0,0.01,0.01,0.01]


dDot = model()
dDot.loadModelFile(modelFile)

params = allHeightsStartingParams[0]

#######################################
height = params[0]
dDot.setParam('SiGe',height)
currentVolts = params[1]

while autoTuneRunning:
    currentIterationFolder = parentFolder + '/autoTuneIteration_'+str(autoTuneIteration)
    if not os.path.exists(currentIterationFolder):
        os.makedirs(currentIterationFolder)
    
    
    #Generate base COMSOL data
    dataFolder = currentIterationFolder + '/COMSOL_data'
    
    if not os.path.exists(dataFolder):
        os.makedirs(dataFolder)
    
    for i,gateNum in enumerate(order):
        dDot.setParam('pot'+str(gateNum),currentVolts[i])
    
    runNum = 0
    
    startTime = time()
    dDot.runSimulation('std1')
    
    dDot.exportData(dataFolder+'/Run_'+str(runNum)+'.txt','data1')
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
        dDot.exportData(dataFolder+'/Run_'+str(runNum)+'.txt','data3',i+1)
        runNum +=1
    
    print "Subsequent "+str(len(order))+" runs: ",time()-startTime
    

    #Generate CODA capacitive Matrix
    
    codaDataFolder = currentIterationFolder + '/CODA_data'
    
    if not os.path.exists(codaDataFolder):
        os.makedirs(codaDataFolder)
    outputFile = open(codaDataFolder+'/CODA_data.txt','w')
    
    for i in range(len(order)+1):
        density,ts = analyzeDataFile(dataFolder+'/Run_'+str(i)+'.txt')
        outputVolts = currentVolts
        if i>0:
            outputVolts[i-1] += 0.001
        outputFile.write(str(outputVolts+density.tolist()+ ts.tolist()).replace(',','').replace(']','').replace('[','')+' \n')
        if i>0:
            outputVolts[i-1] -= 0.001
    outputFile.close()

    #Run CODA, select certain solutions
    plot,voltageChanges = l1ErrorPlot(codaDataFolder+'/CODA_data.txt',target)
    numTries = 10
    
    lastError = 0.
    paramMat = []
    dataIndecesUsed = []
    for i in range(numTries):
        index = int(i*1./numTries*len(plot[0]))
        if plot[1][index]==lastError:
            continue
        print plot[0][index],plot[1][index],plot[2][index]
        paramMat += [array(currentVolts)]
        paramMat[-1]+= voltageChanges[index]
        lastError = plot[1][index]
        dataIndecesUsed+= [index]
    
    numTries = len(paramMat)
    paramMat = array(paramMat)
    paramMat = paramMat.transpose()
    
    print paramMat
    
    #Test actual convergence of selected solutions

    gateList = []
    for gateNum in order:
        gateList += ['pot'+str(gateNum)]

    dDot.initializeSweep(gateList,paramMat)
    
    startTime = time()
    dDot.runSimulation('std2')
    
    for i in range(numTries):
        dDot.exportData(dataFolder+'/Run_'+str(runNum)+'.txt','data3',i+1)
        runNum +=1
    
    print "Subsequent "+str(numTries)+" runs: ",time()-startTime
    

    #Select optimal step and update parameters
    convergenceData = [[],[]]
    minDistance = inf
    minDistanceDataIndex = 0
    minDistanceSystemParams = None
    for i,dataIndex in enumerate(dataIndecesUsed):
        convergenceData[0] += [plot[0][dataIndex]]
        density,ts = analyzeDataFile(dataFolder+'/Run_'+str(len(order)+1+i)+'.txt')
        pointDistance = sum(abs(density[0:2]-target[0:2]))
        pointDistance += sum(abs(array(map(log,ts))-array(map(log,target[4:]))))
        convergenceData[1] += [pointDistance*plot[2][dataIndex]/plot[1][dataIndex]]
        if pointDistance <minDistance:
            minDistanceDataIndex = dataIndex
            minDistanceSystemParams = (density,ts)
            minDistance = pointDistance
    
    
    if plotL1:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.cla()
        plt.clf()
        plt.plot(plot[0],plot[2],linewidth=2.)
        plt.plot(convergenceData[0],convergenceData[1],'s--',linewidth=2.)
        labels = ['Predicted Error','Actual Error']
        plt.legend(labels,loc='upper left')
        plt.xlabel('L1 Norm (V)')
        plt.ylabel('Percent Error')
        plt.title('Predicted Error and Actual Error for '+str(int(height*1000.))+'nm thick SiGe Buffer at Autotune step '+str(autoTuneIteration))
        plt.savefig(codaDataFolder+'/predictedError_vs_actual.png')
    
    currentVolts = array(currentVolts)
    currentVolts += voltageChanges[minDistanceDataIndex]
    currentVolts = currentVolts.tolist()
    print 'Completed step '+str(autoTuneIteration)
    print 'The voltages:'
    print currentVolts
    print 'yield the system state:'
    print minDistanceSystemParams
    
    print 'Current distance from target: '+ str(minDistance)
    
    autoTuneIteration += 1


    if minDistance < 0.1:
        print 'Finished AutoTune!'
        print 'The voltages:'
        print currentVolts
        print 'yield the system state:'
        print minDistanceSystemParams
        autoTuneRunning = False
        
        

