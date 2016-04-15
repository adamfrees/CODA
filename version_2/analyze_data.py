from numpy import *
from regularize import regularize
from numpy.linalg import norm
from operator import add
import os
import pickle
from random import random

def getMetric(filename,target,acceptanceError=0.5):
    data = regularize(filename,15,[4,5,6])
    data.makeRelativeData()
    data.declareTarget(target)
    #data.scaleTargets([4,5,6],10.0)
    data.useTargetIndices([0,1,4,5,6])
    data.target = array([0.,1.,0.,0.,0.])
    data.constructModel()
    data.compressedSense()
    #dataListOccupation = data.createNormList() 
    pickle.dump(data.results,open(filename+'_results_pickle.p','wb'))
    for i,x in enumerate(data.results):
        #print x
        if x[1]>acceptanceError*norm(data.target):
        	index = i
        	break
    #voltages = map(lambda x: x*(abs(x)>1e-4),data.correspondingVecs)
    
    #for i,x in enumerate(voltages):
    #    if sum(x)!=0.:
    #        index = i
    #        break
    #index = index #+ (len(voltages)-index)/2
    return  data.results[index][2]
    
target = [1.,2.,0,0,0.01,0.01,0.01]

acc = 0.95

path = "/home/frees/Dropbox/_UW/CODA/Var_Height"
folder = "/Condensed_data"

#numHeights = 0
listOfFiles = []

for subdir, dirs, files in os.walk(path+folder):
    for file in files:
        if file.startswith('height_'):
            listOfFiles += [os.path.join(subdir, file)]
#numHeights = max(numHeights,int(file.split('_')[1]))

data = zeros((len(listOfFiles),2))

for i,file in enumerate(listOfFiles):
    l1Norm = getMetric(file,target,acceptanceError=acc)
    data[i,0] = int(file.strip(path).strip(folder).strip('/').split('_')[1])
    data[i,1] = l1Norm
    break
'''
data = array(sorted(data, key=lambda a_entry: a_entry[0]))

heightFile = "/heightList.txt"

heightList = []

f = open(path+heightFile,'r')

for line in f:
    heightList += [float(line.split(' ')[1])]

f.close()

f = open(path+'/'+str(acc)+'_metric.txt','w')

for a in data:
    f.write(str(heightList[int(a[0])])+','+str(a[1])+'\n')
    
f.close()

pickle.dump(data,open(path+'/'+str(acc)+'_metric.p','wb'))'''
