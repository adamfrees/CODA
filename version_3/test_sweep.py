from numpy import *
from model import model
from time import time

dataFolder = '/home/frees/Dropbox/_UW/CODA/version3/Data'
modelFile = '/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height_2.mph'

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

#startingVolts = [0.411818943243,0.341810107612,0.398860825154,-0.35,0.121025101459,-0.39,-0.099766291307,-0.0118656225221,-0.370517974419,0.120554018307,-0.35,-0.1,-0.0951015486058,-0.100766615989,-0.1]
startingVolts = [0.411818943243,0.301810107612,0.398860825154,-0.35,0.118525101459,-0.355,-0.097766291307,-0.0118656225221,-0.34,0.119554018307,-0.35,-0.1,-0.0951015486058,-0.100766615989,-0.1]

dDot = model()
dDot.loadModelFile(modelFile)

for i,gateNum in enumerate(order):
    dDot.setParam('pot'+str(gateNum),startingVolts[i])

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

print "Subsequent 15 runs: ",time()-startTime
