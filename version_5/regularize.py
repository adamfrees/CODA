from numpy import *
from sklearn import linear_model
import csv
import sys
import cvxpy as cvx
from numpy.linalg import norm
set_printoptions(precision=5,linewidth=1000,suppress=True)

class regularize:
    def __init__(self,inputFile,ind,logIndexList = [],method='L1'):
        dataFile = open(inputFile,"rb")
        data = reader=csv.reader(dataFile,delimiter=' ')
        x=list(reader)
        data = array(x)
        data=data[:,0:-1].astype('float')
        dataFile.close()
        self.X = data[:,:ind] #voltage data
        #print sum(abs(self.X[0,:]))
        self.Y = data[:,ind:] #target data 
        #print self.X
        self.ind = ind
        self.dep =self.Y.shape[1]
        self.controlVecs = []
        self.logIndexList = logIndexList
        self.startDeps = []
        self.method = method
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
        self.M = linalg.lstsq(self.X,self.Y)[0].T
        self.distance = norm(self.target,ord=2)
        print self.target
        print self.distance
        print norm(self.target.T,ord=2)

    def findDirection(self):
        resultList = []
        if self.method=='L1':
            x = cvx.Variable(self.M.shape[1])
            print self.M.shape[1]            
            # Create two constraints.
            constraints = [self.M*x == self.target]
            
            # Form objective.
            obj = cvx.Minimize(cvx.norm1(x))
            
            # Form and solve problem.
            prob = cvx.Problem(obj, constraints)
            prob.solve()  # Returns the optimal value.
            print "status:", prob.status
            print "optimal value", prob.value
            print "optimal var", x.value
            self.voltages = array(x.value.T)[0]
        elif self.method=='L2':
            #self.voltages = array(linalg.lstsq(self.M,self.target)[0]).T[0]
            x = cvx.Variable(self.M.shape[1])
            print self.M.shape[1]
            # Create two constraints.
            constraints = [self.M*x == self.target]

            # Form objective.
            obj = cvx.Minimize(cvx.norm2(x))

            # Form and solve problem.
            prob = cvx.Problem(obj, constraints)
            prob.solve()  # Returns the optimal value.
            print "status:", prob.status
            print "optimal value", prob.value
            print "optimal var", x.value
            self.voltages = array(x.value.T)[0]
        else:
            print "Warning: method not recognized." 

