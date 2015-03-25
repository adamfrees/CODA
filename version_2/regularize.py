from numpy import *
from sklearn import linear_model
import csv
import sys
set_printoptions(precision=5,linewidth=1000,suppress=True)

class regularize:
    def __init__(self,inputFile,ind,logIndexList = []):
        dataFile = open(inputFile,"rb")
        data = reader=csv.reader(dataFile,delimiter=' ')
        x=list(reader)
        data = array(x)
        data=data[:,0:-1].astype('float')
        dataFile.close()
        self.X = data[:,:ind] #voltage data
        self.Y = data[:,ind:] #target data 
        print self.X
        print self.Y
        self.ind = ind
        self.dep =self.Y.shape[1]
        self.controlVecs = []
        self.logIndexList = logIndexList
        self.startDeps = []
        for index in logIndexList:
            self.Y[:,index] = log(self.Y[:,index])
    
    def declareTarget(self,target):
        self.target = array(target).reshape((self.dep,1))#-self.startDeps
        for index in self.logIndexList:
            self.target[index,0] = log(self.target[index,0])#log((exp(self.startDeps[index,0])+self.target[index,0])/exp(self.startDeps[index,0]))
        self.target = self.target-self.startDeps
    
    def scaleTargets(self,indexList,multiplier):
        for index in indexList:
            self.target[index,0] *= multiplier
            self.Y[:,index] *= multiplier
    
    def useTargetIndices(self,indexList):
        newY = zeros((self.Y.shape[0],len(indexList)))
        newTarget = zeros((len(indexList),1))
        for i,index in enumerate(indexList):
            newY[:,i] = self.Y[:,index]
            newTarget[i,0] = self.target[index,0]
        self.Y = newY
        self.target = newTarget
          
    def makeRelativeData(self):
        self.startDeps = self.Y[0,:].reshape((self.dep,1))
        self.X = self.X[1:,:] - self.X[0,:]
        self.Y = self.Y[1:,:] - self.Y[0,:]
        
    def constructModel(self):
        for y in nditer(self.Y, op_flags=['readwrite']):
            if y==-inf:
                y[...] = -sys.float_info.max
        print self.Y
        self.M = linalg.lstsq(self.X,self.Y)[0].T
        print self.M

    def compressedSense(self,alphaRange=[1e-5,1000.0],numSamples=7000,max_iter=50000,tol=1e-6):
        resultList = []
        for alpha2 in logspace(log(alphaRange[0])/log(10),log(alphaRange[1])/log(10),numSamples):
            clf = linear_model.Lasso(alpha=alpha2,fit_intercept=False,max_iter=max_iter,tol=tol)
            clf.fit(self.M,self.target)    
            error = sum(abs(self.target-dot(array(self.M),clf.coef_.reshape((self.ind,1)))))
            oneNorm = sum(abs(clf.coef_))
            zeroNorm = linalg.norm(clf.coef_,ord=0)
            resultList += [[alpha2,error,oneNorm,zeroNorm]]    
            self.controlVecs += [clf.coef_[:]]
        self.controlVecs = array(self.controlVecs)
        self.results = array(resultList)
        

    def createNormList(self):
        errorSlice = self.results[:,1:4]
        sortedVecs = self.controlVecs[errorSlice[:,1].argsort()]
        errorSlice = errorSlice[errorSlice[:,1].argsort()]
        oneNormVec = logspace(log(1e-5)/log(10),log(5.0)/log(10),2000)    
        returnList = []        
        self.correspondingVecs = []
        for oneNormThreshold in oneNormVec:
            if amax(errorSlice[:,1])<=oneNormThreshold:
                indexThreshold = len(errorSlice[:,1])
            else:        
                indexThreshold = argmax(errorSlice[:,1]>oneNormThreshold)
            if indexThreshold==0:
                errorVal = 100
                zeroNormVal = 100
            else:
                minIndex = argmin(errorSlice[:indexThreshold,0])
                errorVal = errorSlice[minIndex,0]
                zeroNormVal = errorSlice[minIndex,2]
                
            returnList += [[oneNormThreshold*1000.0,errorVal,zeroNormVal]]
            try:
                self.correspondingVecs +=  [sortedVecs[minIndex]]
            except:
                self.correspondingVecs += [0.0*sortedVecs[0]]
        return returnList
    
    def computeErrors(self,resultList,controlVecs):
        newErrorList = []
        for i,result in enumerate(resultList):
            error = sum(abs(self.target-dot(array(self.M),controlVecs[i].reshape((self.ind,1)))))
            newErrorList += [error]
        return newErrorList
