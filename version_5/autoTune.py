from numpy import *
import sys
sys.path.append('/home/frees/Code/comsol-python')
from comsolpy_matlab_engine import model
from time import time
from regularize import regularize
import os
from analyze_data import analyzeDataFile 
from numpy.linalg import norm

useCOMSOL = True
plotL1 = True

autoTuneRunning = True
autoTuneIteration = 0

modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height_3.mph'

allHeightsStartingParams = [
    [0.03,[0.4,0.45,0.4,-0.4,0.0722,-0.175,-0.0087,-0.012,-0.094,0.0642,-0.4,-0.2,-0.25,-0.264,-0.2],'L1'],
    [0.032,[0.4,0.45,0.4,-0.4,0.075,-0.18,-0.01,-0.012,-0.1,0.0685,-0.4,-0.2,-0.26,-0.26,-0.2],'L1'],
    [0.034,[0.4,0.45,0.4,-0.4,0.077,-0.177,-0.012,-0.012,-0.105,0.073,-0.4,-0.2,-0.27,-0.256,-0.2],'L1'],
    [0.036,[0.4,0.45,0.4,-0.4,0.079,-0.1765,-0.014,-0.012,-0.109,0.0775,-0.4,-0.2,-0.285,-0.252,-0.2],'L1'],
    [0.038,[0.4,0.45,0.4,-0.4,0.082,-0.179,-0.016,-0.012,-0.114,0.081,-0.4,-0.2,-0.295,-0.248,-0.2],'L1'],
    [0.04,[0.4,0.45,0.4,-0.4,0.084,-0.183,-0.0177,-0.012,-0.12,0.0855,-0.4,-0.2,-0.305,-0.24,-0.2],'L1']
    #[0.036,[0.4,0.45,0.4,-0.4,0.079,-0.1765,-0.014,-0.012,-0.109,0.0775,-0.4,-0.2,-0.285,-0.252,-0.2],'L1'],
    #[0.038,[0.4,0.45,0.4,-0.4,0.082,-0.179,-0.016,-0.012,-0.114,0.081,-0.4,-0.2,-0.295,-0.248,-0.2],'L1']
    #[0.04,[0.4,0.45,0.4,-0.4,0.084,-0.183,-0.0177,-0.012,-0.12,0.0855,-0.4,-0.2,-0.305,-0.24,-0.2],'L1']
    ]

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

target = [1.,2.,0,0,0.01,0.01,0.01]

if useCOMSOL:
    dDot = model(comsolBinLocation='')
    dDot.loadModelFile(modelFile)

#params = allHeightsStartingParams[0]
for params in allHeightsStartingParams:
    autoTuneRunning = True
    autoTuneIteration = 0
    #######################################
    method = params[2]
    height = params[0]
    parentFolder = '/home/frees/Dropbox/_UW/CODA/version_5/Autotune/'+str(int(height*1000))+'_'+method
    if not os.path.exists(parentFolder):
            os.makedirs(parentFolder)
    
    def findDirection(filename,target):
        data = regularize(filename,15,[4,5,6],method=method)
        data.makeRelativeData()
        data.declareTarget(target)
        data.useTargetIndices([0,1,4,5,6])
        data.scaleTargets([2,3,4],0.5)
        data.constructModel()
        data.findDirection()
        return data.voltages,data.distance
    
    def predictDeps(filename,changes,target=target):
        data = regularize(filename,15,[4,5,6],method=method)
        data.makeRelativeData()
        data.declareTarget(target)
        data.useTargetIndices([0,1,4,5,6])
        data.scaleTargets([2,3,4],0.5)
        data.constructModel()
        #print data.M,changes
        deps = dot(data.M,changes)
        #print deps
        for i,x in enumerate([0,1,4,5,6]):
            deps[i] += data.startDeps[x]
        return deps
    
    if useCOMSOL:
        dDot.setParam('SiGe',height)
    currentVolts = params[1]
    
    while autoTuneRunning:
        currentIterationFolder = parentFolder + '/autoTuneIteration_'+str(autoTuneIteration)
        
        if not os.path.exists(currentIterationFolder):
            os.makedirs(currentIterationFolder)
        
        codaDataFolder = currentIterationFolder + '/CODA_data' 
        
        #Generate base COMSOL data
        dataFolder = currentIterationFolder + '/COMSOL_data'
        
        if not os.path.exists(dataFolder):
            os.makedirs(dataFolder)
        if useCOMSOL: 
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
            
            voltageStepSize = 1.e-4
            
            paramMat += voltageStepSize*eye(len(order))
            
            dDot.initializeSweep(gateList,paramMat)
            
            startTime = time()
            dDot.runSimulation('std2')
            
            for i in range(len(order)):
                dDot.exportData(dataFolder+'/Run_'+str(runNum)+'.txt','data3',i+1)
                runNum +=1
            
            print "Subsequent "+str(len(order))+" runs: ",time()-startTime
            
            
            #Generate CODA capacitive Matrix
            
            if not os.path.exists(codaDataFolder):
                os.makedirs(codaDataFolder)
            outputFile = open(codaDataFolder+'/CODA_data.txt','w')
            
            for i in range(len(order)+1):
                density,ts = analyzeDataFile(dataFolder+'/Run_'+str(i)+'.txt')
                outputVolts = currentVolts
                if i>0:
                    outputVolts[i-1] += voltageStepSize
                outputFile.write(str(outputVolts+density.tolist()+ ts.tolist()).replace(',','').replace(']','').replace('[','')+' \n')
                if i>0:
                    outputVolts[i-1] -= voltageStepSize
            outputFile.close()
        
        #Run CODA, select certain solutions
        voltageChanges,totalDistance = findDirection(codaDataFolder+'/CODA_data.txt',target)
        solutionL2Norm = norm(voltageChanges,ord=2) 
        paramMat = []
        convergenceData = [[],[]]
        
        minDistance = 10.
        errorCounter = 0
        for alpha in logspace(-4,0,11):
            for i,gateNum in enumerate(order):
                print alpha*voltageChanges[i]
                dDot.setParam('pot'+str(gateNum),currentVolts[i]+alpha*voltageChanges[i])
            
            startTime = time()
            dDot.runSimulation('std1')
            
            dDot.exportData(dataFolder+'/Run_'+str(runNum)+'.txt','data1')
            print "Run number "+str(runNum)+": ", time()-startTime
            
            density,ts = analyzeDataFile(dataFolder+'/Run_'+str(runNum)+'.txt')
            predictedDeps = predictDeps(codaDataFolder+'/CODA_data.txt',alpha*voltageChanges)
            differenceFromLinear = abs(density[0:2]-predictedDeps[0:2]).tolist()
            differenceFromLinear += (0.5*abs(array(map(log,ts))-predictedDeps[2:])).tolist()
            differenceFromLinearNorm = norm(array(differenceFromLinear),ord=2)
            normalizedDistance = differenceFromLinearNorm/totalDistance
            if normalizedDistance>0.2:
                errorCounter += 1
            else:
                errorCounter = 0
                optimalSystemParams = [density,ts]
                optimalVoltageChanges = alpha*array(voltageChanges)
                convergenceData[0] += [alpha*solutionL2Norm]
                convergenceData[1] += [normalizedDistance]
                attemptDistance = abs(density[0:2]-target[0:2]).tolist()
                attemptDistance += (0.5*abs(array(map(log,ts))-array(map(log,target[4:])))).tolist()
                attemptDistanceNorm = norm(array(attemptDistance),ord=2)
            if errorCounter>1:
                break
            runNum +=1
            
        if plotL1:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            plt.cla()
            plt.clf()
            plt.plot(convergenceData[0],convergenceData[1],'s--',linewidth=2.)
            plt.xlabel('L2 Length (V)')
            plt.ylabel('Percent Error')
            plt.title('Error for '+str(int(height*1000.))+'nm thick SiGe Buffer at Autotune step '+str(autoTuneIteration))
            plt.savefig(codaDataFolder+'/error_along_direction.png')
        
        currentVolts = array(currentVolts)
        currentVolts += optimalVoltageChanges
        currentVolts = currentVolts.tolist()
        print 'Completed step '+str(autoTuneIteration)
        print 'The voltages:'
        print currentVolts
        print 'yield the system state:'
        print optimalSystemParams
        
        print 'Current distance from target: '+ str(attemptDistanceNorm)
        
        autoTuneIteration += 1
    
    
        if attemptDistanceNorm < 0.05:
            print 'Finished AutoTune!'
            print 'The voltages:'
            print currentVolts
            print 'yield the system state:'
            print optimalSystemParams
            autoTuneRunning = False
            
    
