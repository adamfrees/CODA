from numpy import *
from model import model
from time import time
from regularize import regularize
import os
from analyze_data import analyzeDataFile 

dataFolder = '/home/frees/Dropbox/_UW/CODA/version3/Convergent_Data/0_005'
codaDataFolder = '/home/frees/Dropbox/_UW/CODA/version3/SiGeSweep'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height_3.mph'

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

params = [
    [0.03,[0.4,0.45,0.4,-0.4,0.0722,-0.175,-0.0087,-0.012,-0.094,0.0642,-0.4,-0.2,-0.25,-0.264,-0.2]],
    [0.032,[0.4,0.45,0.4,-0.4,0.075,-0.18,-0.01,-0.012,-0.1,0.0685,-0.4,-0.2,-0.26,-0.26,-0.2]],
    [0.034,[0.4,0.45,0.4,-0.4,0.077,-0.177,-0.012,-0.012,-0.105,0.073,-0.4,-0.2,-0.27,-0.256,-0.2]]#,
    #[0.036,[0.4,0.45,0.4,-0.4,0.079,-0.1765,-0.014,-0.012,-0.109,0.0775,-0.4,-0.2,-0.285,-0.252,-0.2]],
    #[0.038,[0.4,0.45,0.4,-0.4,0.082,-0.179,-0.016,-0.012,-0.114,0.081,-0.4,-0.2,-0.295,-0.248,-0.2]],
    #[0.04,[0.4,0.45,0.4,-0.4,0.084,-0.183,-0.0177,-0.012,-0.12,0.0855,-0.4,-0.2,-0.305,-0.24,-0.2]]
    ]

def l1ErrorPlot(filename,target):
    data = regularize(filename,15,[4,5,6])
    data.makeRelativeData()
    data.declareTarget(target)
    #data.scaleTargets([4,5,6],10.0)
    data.useTargetIndices([0,1,4,5,6])
    data.constructModel()
    data.compressedSense(alphaRange=[1.e-4,100.],numSamples=100)
    dataListOccupation = data.createNormList()
    distance = sum(abs(data.target)) 
    toReturn = [[],[],[]]
    voltages = []
    for index in range(len(data.results)-2,1,-1):
        if data.results[index][1]>data.results[index+1][1]:
            data.results[index][1]=data.results[index+1][1]
            data.controlVecs[index] = data.controlVecs[index+1]
    for index,x in enumerate(data.results):
        if x[1]<distance and x[1]>0.5*distance:
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

target = [1.,2.,0,0,0.01,0.01,0.01]
data = []
voltageChanges = []

for heightParams in params:
    height = heightParams[0]
    fileName = codaDataFolder+'/height_'+str(height)+'.txt'
    plot = l1ErrorPlot(fileName,target)
    data += [plot[0]]
    voltageChanges += [plot[1]]

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

cycle = []
labels = []
for i,param in enumerate(params):
    if len(params)>1:
        cycle += [getColor(1./(len(params)-1.)*i)]
    else:
        cycle = [(0.,0.,1.)]
    labels += [str(int(param[0]*1000.))+"nm SiGe Buffer"]

plt.gca().set_color_cycle(cycle)

for x in data:
    plt.plot(x[0],map(lambda a: log(a),x[1]),linewidth=2.)

#plt.xlim(0.,0.15)
plt.legend(labels,loc='upper right')
plt.xlabel('L1 Norm (V)')
plt.ylabel('Log(Error)')
plt.title('L1-Error Plot')
plt.savefig(codaDataFolder+'/L1-Error.png')
plt.cla()
plt.clf()

plt.gca().set_color_cycle(cycle)

for x in data:
    plt.plot(x[0],x[2],linewidth=2.)

#plt.xlim(0.,0.15)
plt.legend(labels,loc='upper right')
plt.xlabel('L1 Norm (V)')
plt.ylabel('Percent Error')
plt.title('Percent Error Plot')
plt.savefig(codaDataFolder+'/percent_Error.png')

numTries = [10]*len(params)

dataIndecesUsed = [[] for i in range(len(params))]

paramMats = []
for dataIndex,x in enumerate(data):
    paramList = params[dataIndex][1]
    paramMat = []#array([paramList]*numTries)
    lastError = 0.
    for i in range(numTries[dataIndex]):
        index = int(i*1./numTries[dataIndex]*len(x[0]))
        if x[1][index]==lastError:
            continue
        print x[0][index],x[1][index],x[2][index]
        paramMat += [array(paramList)]
        paramMat[-1]+= voltageChanges[dataIndex][index]
        lastError = x[1][index]
        dataIndecesUsed[dataIndex]+= [index]
    numTries[dataIndex] = len(paramMat)
    paramMat = array(paramMat)
    paramMat = paramMat.transpose()
    paramMats += [paramMat]

print paramMats

for i,param in enumerate(params):
    if len(params)>1:
        cycle += [getColor(1./(len(params)-1.)*i),getColor(1./(len(params)-1.)*i)]
    else:
        cycle = [(0.,0.,1.),(0.,0.,1.)]

plt.gca().set_color_cycle(cycle)

dDot = model()
dDot.loadModelFile(modelFile)

for index,heightParams in enumerate(params):
    height = heightParams[0]
    startingVolts = heightParams[1]
    dDot.setParam('SiGe',height)
    for i,gateNum in enumerate(order):
        dDot.setParam('pot'+str(gateNum),startingVolts[i])
    
    runNum = 0
    
    startTime = time()
    dDot.runSimulation('std1')
    if not os.path.exists(dataFolder+'/Height_'+str(height)):
        os.makedirs(dataFolder+'/Height_'+str(height))
    dDot.exportData(dataFolder+'/Height_'+str(height)+'/Run_'+str(runNum)+'.txt','data1')
    runNum +=1
    print "First run: ", time()-startTime
    
    gateList = []
    paramList = []
    for gateNum in order:
        gateList += ['pot'+str(gateNum)]
    paramMat = paramMats[index]

    dDot.initializeSweep(gateList,paramMat)
    
    startTime = time()
    dDot.runSimulation('std2')
    
    for i in range(numTries[index]):
        dDot.exportData(dataFolder+'/Height_'+str(height)+'/Run_'+str(runNum)+'.txt','data3',i+1)
        runNum +=1
    
    print "Subsequent "+str(numTries[index])+" runs: ",time()-startTime

for index,heightParams in enumerate(params):
    convergenceData = [[],[]]
    for i,dataIndex in enumerate(dataIndecesUsed[index]):
        convergenceData[0] += [data[index][0][dataIndex]]
        density,ts = analyzeDataFile(dataFolder+'/Height_'+str(heightParams[0])+'/Run_'+str(i+1)+'.txt')
        print convergenceData[0][-1], density, ts
        pointDistance = sum(abs(density[0:2]-target[0:2]))
        pointDistance += sum(abs(array(map(log,ts))-array(map(log,target[4:]))))
        convergenceData[1] += [pointDistance*data[index][2][dataIndex]/data[index][1][dataIndex]]
    plt.cla()
    plt.clf()
    
    plt.plot(data[index][0],data[index][2],linewidth=2.)
    plt.plot(convergenceData[0],convergenceData[1],'s--',linewidth=2.)
    labels = ['Predicted Error','Actual Error']
    #plt.xlim(0.,0.15)
    plt.legend(labels,loc='upper left')
    plt.xlabel('L1 Norm (V)')
    plt.ylabel('Percent Error')
    plt.title('Predicted Error and Actual Error for '+str(int(heightParams[0]*1000.))+'nm thick SiGe Buffer')
    plt.savefig(codaDataFolder+'/predictedError_vs_actual_0_005_'+str(heightParams[0])+'.png')
