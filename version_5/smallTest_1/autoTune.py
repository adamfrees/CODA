from numpy import *
import sys
sys.path.append('/Users/adamfrees/GitRepos/comsol-python')
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

modelFile = '/Users/adamfrees/Dropbox/_UW/Scale-up/Comsol_test_1/simpleSingleDot_w_high_Lead.mph'

target = [2.]

#allHeightsStartingParams=[[30.,[0.190050340468617,-0.2,-0.2],'L1'],[32.,[0.20027523864318902,-0.2,-0.2],'L1'],[34.,[0.21058880122396614,-0.2,-0.2],'L1'],[36.,[0.221065369580485,-0.2,-0.2],'L1'],[38.,[0.23163401295215844,-0.2,-0.2],'L1']]#simpleSingleDot.mph
#allHeightsStartingParams=[[30.,[0.13577273372369625,-0.2,-0.2],'L1'],[32.,[0.14137342528783453,-0.2,-0.2],'L1'],[34.,[0.14683575221584297,-0.2,-0.2],'L1'],[36.,[0.15213096912565122,-0.2,-0.2],'L1'],[38.,[0.15723684487136977,-0.2,-0.2],'L1']]#simpleSingleDot_w_Lead.mph#voltLead=0.6
#allHeightsStartingParams=[[30.,[0.12558932329196001,-0.2,-0.2],'L1'],[32.,[0.1298659398614813,-0.2,-0.2],'L1'],[34.,[0.13402961795887958,-0.2,-0.2],'L1'],[36.,[0.13824949682438456,-0.2,-0.2],'L1'],[38.,[0.1423095716861876,-0.2,-0.2],'L1']]#simpleSingleDot_w_Lead.mph#voltLead=1.
allHeightsStartingParams=[[30.,[0.1422856442898976,-0.2,-0.2],'L1'],[32.,[0.14873661066908628,-0.2,-0.2],'L1'],[34.,[0.1551211581300123,-0.2,-0.2],'L1'],[36.,[0.16143216556424986,-0.2,-0.2],'L1'],[38.,[0.16764987251535599,-0.2,-0.2],'L1']]#simpleSingleDot_w_high_Lead.mph#voltLead=0.6


order=range(1,4)


if useCOMSOL:
    dDot = model(comsolBinLocation='/Applications/COMSOL52/Multiphysics/bin/')
    dDot.loadModelFile(modelFile)

#params = allHeightsStartingParams[0]
for params in allHeightsStartingParams:
    autoTuneRunning = True
    autoTuneIteration = 0
    #######################################
    method = params[2]
    height = params[0]
    parentFolder = '/Users/adamfrees/Dropbox/_UW/CODA/version_5/SmallTest1_w_high_Lead/'+("%.0f" % height)+'nm'
    if not os.path.exists(parentFolder):
            os.makedirs(parentFolder)

    def findDirection(filename,target):
        data = regularize(filename,3,[],method=method)
        data.makeRelativeData()
        data.declareTarget(target)
        data.useTargetIndices([0])
        data.constructModel()
        data.findDirection()
        return data.voltages,data.distance

    def predictDeps(filename,changes,target=target):
        data = regularize(filename,3,[],method=method)
        data.makeRelativeData()
        data.declareTarget(target)
        data.useTargetIndices([0])
        data.constructModel()
        #print data.M,changes
        deps = dot(data.M,changes)
        #print deps
        for i,x in enumerate([0]):
            deps[i] += data.startDeps[x]
        return deps

    if useCOMSOL:
        dDot.setParam('SiGeHeight',height)
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
                dDot.setParam('volt'+str(gateNum),currentVolts[i])

            runNum = 0

            startTime = time()
            dDot.runSimulation('std1')
            densities = [[dDot.evaluateNumerical('int1')]]
            #dDot.exportData(dataFolder+'/Run_'+str(runNum)+'.txt','data1')
            runNum +=1
            print "First run: ", time()-startTime

            gateList = []
            paramList = []
            voltageStepSize = 1.e-4
            startTime = time()
            for gateNum in order:
                newVolt=dDot.getParam('volt'+str(gateNum))+voltageStepSize
                dDot.setParam('volt'+str(gateNum),newVolt)
                dDot.runSimulation('std1')
                densities += [[dDot.evaluateNumerical('int1')]]
                dDot.setParam('volt'+str(gateNum),newVolt-voltageStepSize)
                runNum +=1

            print "Subsequent "+str(len(order))+" runs: ",time()-startTime


            #Generate CODA capacitive Matrix

            if not os.path.exists(codaDataFolder):
                os.makedirs(codaDataFolder)
            outputFile = open(codaDataFolder+'/CODA_data.txt','w')

            for i in range(len(order)+1):
                density = array(densities[i])
                outputVolts = currentVolts
                if i>0:
                    outputVolts[i-1] += voltageStepSize
                outputFile.write(str(outputVolts+density.tolist()).replace(',','').replace(']','').replace('[','')+' \n')
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
                dDot.setParam('volt'+str(gateNum),currentVolts[i]+alpha*voltageChanges[i])

            startTime = time()
            dDot.runSimulation('std1')

            #dDot.exportData(dataFolder+'/Run_'+str(runNum)+'.txt','data1')
            print "Run number "+str(runNum)+": ", time()-startTime

            #density,ts = analyzeDataFile(dataFolder+'/Run_'+str(runNum)+'.txt')
            density=[dDot.evaluateNumerical('int1')]
            predictedDeps = predictDeps(codaDataFolder+'/CODA_data.txt',alpha*voltageChanges)
            differenceFromLinear = [abs(density[0]-predictedDeps[0])]
            #differenceFromLinear += (0.5*abs(array(map(log,ts))-predictedDeps[2:])).tolist()
            differenceFromLinearNorm = norm(array(differenceFromLinear),ord=2)
            normalizedDistance = differenceFromLinearNorm/totalDistance
            if normalizedDistance>0.2:
                errorCounter += 1
            else:
                errorCounter = 0
                optimalSystemParams = [density]
                optimalVoltageChanges = alpha*array(voltageChanges)
                convergenceData[0] += [alpha*solutionL2Norm]
                convergenceData[1] += [normalizedDistance]
                attemptDistance = [abs(density[0]-target[0])]
                #attemptDistance += (0.5*abs(array(map(log,ts))-array(map(log,target[4:])))).tolist()
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


        if attemptDistanceNorm < 0.01:
            print 'Finished AutoTune!'
            print 'The voltages:'
            print currentVolts
            print 'yield the system state:'
            print optimalSystemParams
            autoTuneRunning = False
