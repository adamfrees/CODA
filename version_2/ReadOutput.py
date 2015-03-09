import sys
import os
from numpy import *
import scipy.sparse as sps
from scipy.sparse.linalg import eigs
import scipy.ndimage as sndi
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize
import pickle
#import matplotlib as mpl
#import matplotlib.cm as cm
#import matplotlib.colors as co
#mpl.use("Agg")
#import matplotlib.pyplot as plt
import time
import csv

# The data processing class:

class processData:
    def __init__(self,potentialDataName):
        self.mt = 0.19*9.1e-31 # transverse EM in SI
        self.hbar = 1.0546e-34 # hbar in SI
        self.q = 1.602e-19 # elementary charge
        self.eVtoJ = 1.602e-19 #eV in SI
    
        xVec,yVec, potentialData, densityData = self.readDataFile(potentialDataName)
        
        self.xVec = xVec
        self.yVec = yVec
        self.Lx = xVec[-1]-xVec[0]
        self.Ly = yVec[-1]-yVec[0]
        self.effectivePotential = potentialData
        self.lateralDensity = densityData
        
        #implot = plt.imshow(transpose(densityData),origin='lower',aspect=self.Ly/self.Lx)
        #plt.savefig('density.png')
    	#plt.cla()
    	#plt.clf()
    	#implot = plt.imshow(transpose(potentialData),origin='lower',aspect=self.Ly/self.Lx)
        #plt.savefig('potential.png')
    	#plt.cla()
    	#plt.clf()
    	
        totalElectrons =trapz(trapz(densityData,xVec),yVec)       


 ##Function takes in data file, outputs xVec, yVec, in nanometers!!! potentialData in V, Density in number/nm^2
    def readDataFile(self,fileName):
		file = map(lambda x:x.split(),open(fileName,'r').readlines())
		startReadingVoltage=False
		startReadingDensity=False
		for i,line in enumerate(file):
			if len(line)<=1:
				continue
			if line[1]=='Grid':
				xHolder = array(map(float,file[i+1]))*1e3
				yHolder = array(map(float,file[i+2]))*1e3
				xDim = xHolder.shape[0]
				yDim = yHolder.shape[0]
				voltageHolder = zeros((xDim,yDim)); densityHolder = zeros((xDim,yDim))
			elif line[1]=='-(V+0.0)':
				startReadingVoltage=True
				yIndex=0; xIndex=0
			elif line[1]=='Data': continue
			#elif line[1]=='2*0.19*me_const/(pi*hbar_const^2)*(e_const*(int1(root.x[Mm],root.y[Mm])+0.0))*(-0.0<int1(root.x[Mm],root.y[Mm]))*(0>=y)+2*0.19*me_const/(pi*hbar_const^2)*(e_const*(V+0.0))*(-0.0<V)*(0<y)':
			elif line[1]=='2*0.19*me_const/(pi*hbar_const^2)*(e_const*(V+0.0))*(-0.0<V)':
				startReadingDensity=True
				yIndex=0; xIndex=0
			elif startReadingVoltage and not startReadingDensity:
				voltageHolder[:,yIndex]=array(map(float,line))
				yIndex +=1
			elif startReadingDensity:
				densityHolder[:,yIndex]=array(map(float,line))
				yIndex+=1
		outputLines = ""
		for i,y in enumerate(yHolder):
			for j,x in enumerate(xHolder):
				outputLines += str(x)+" "+str(y)+" "+str(voltageHolder[j,i])+"\n"
		#outputFile = open('Volts.txt','w')
		#outputFile.writelines(outputLines)
		#outputFile.close()
		return xHolder,yHolder,voltageHolder,densityHolder
    
    
    def processPotential(self,guesses,tunnelingPairs):
        # First, calculate the target region centers based on the guesses provided:
        self.findCenters()
        self.identifyCenters(guesses)

        # Next, specify the watershed partition:
        self.watershedDoubleDot()
        self.determineWatershedRegions()

        # Finally, compute the quantities of interest:
        self.computeDensities()
        self.computeTunnelings(tunnelingPairs)

########This all has to do with finding the centers of the dots and leads
    def addNeighbors(self,inputSet,inputArray):
        newSet = []
        xLen = len(self.xVec)
        yLen = len(self.yVec)
        for element in inputSet:
            i = element[0]; j=element[1]
            if i>0:
                if inputArray[i-1,j]!=0.0: newSet += [(i-1,j)]
            if i>0 and j<(yLen-1):
                if inputArray[i-1,j+1]!=0.0: newSet += [(i-1,j+1)]
            if i<(xLen-1):
                if inputArray[i+1,j]!=0.0: newSet += [(i+1,j)]
            if j<(yLen-1):
                if inputArray[i,j+1]!=0.0: newSet += [(i,j+1)]
            if j>0:
                if inputArray[i,j-1]!=0.0: newSet += [(i,j-1)]
            if i<(xLen-1) and j<(yLen-1):
                if inputArray[i+1,j+1]!=0.0: newSet += [(i+1,j+1)]
            if i<(xLen-1) and j>0:
                if inputArray[i+1,j-1]!=0.0: newSet += [(i+1,j-1)]
            if i>0 and j>0:
                if inputArray[i-1,j-1]!=0.0: newSet += [(i-1,j-1)]      
        return set(newSet+list(inputSet))   
    def findRegion(self,inputPoint,inputArray):
        testSet = set([inputPoint])
        while True:
            oldLen = len(testSet)
            testSet = self.addNeighbors(testSet,inputArray)
            newLen = len(testSet)
            if oldLen == newLen:
                break
        return testSet
    def createPartitions(self,inputArray):
        coveredPoints = []
        partitions = []
        for i,x in enumerate(self.xVec):
            for j,y in enumerate(self.yVec):
                if inputArray[i,j]==0:
                    coveredPoints +=[(i,j)]
                elif (i,j) in coveredPoints:
                    pass
                else:
                    newRegion = list(self.findRegion((i,j),inputArray))
                    partitions += [newRegion]
                    coveredPoints += newRegion
        partitionPlot =  zeros(inputArray.shape)
        for i,z in enumerate(self.xVec):
            for j,y in enumerate(self.yVec):
                for k,partition in enumerate(partitions):
                    if (i,j) in partition:
                        partitionPlot[i,j] = k+1        
        return partitionPlot,partitions
    def findCenters(self):
        firstDerivAbs = abs(sndi.filters.sobel(self.effectivePotential,axis=0))+abs(sndi.filters.sobel(self.effectivePotential,axis=1))
        firstDerivAbs/= amax(firstDerivAbs)
        derivXX = sndi.filters.sobel(sndi.filters.sobel(self.effectivePotential,axis=0),axis=0)
        derivYY = sndi.filters.sobel(sndi.filters.sobel(self.effectivePotential,axis=1),axis=1)
        derivYX= sndi.filters.sobel(sndi.filters.sobel(self.effectivePotential,axis=1),axis=0)

        minMatrix = zeros(self.lateralDensity.shape)
        for i,x in enumerate(self.xVec):
            for j,y in enumerate(self.yVec):
                if firstDerivAbs[i,j]<=0.1 and derivXX[i,j]>0.0 and (derivXX[i,j]*derivYY[i,j]-derivYX[i,j]**2)>0.0:
                    minMatrix[i,j] = 1.0    
        partitionPlot,partitions = self.createPartitions(minMatrix)      
        centerList = []
        for partition in partitions:
            centerList += [average(array(partition),axis=0)]
        self.centerList = centerList
    def identifyCenters(self,guesses):
        labelList = []
        self.identificationSuccess = True
        for guess in guesses:
            guessLabel = guess[0]
            guessCoords = guess[1]
            distanceList = []
            for center in self.centerList:
                distanceList += [sum(abs(center-guessCoords))]
            coordIndex = distanceList.index(min(distanceList))
            if min(distanceList)>20.0:
                self.identificationSuccess = False
            labelList += [[guessLabel,self.centerList[coordIndex]]]
        self.labelList = labelList
######################    End of Center-finding stuff

    def watershedDoubleDot(self):
        watershedMap = zeros(self.effectivePotential.shape)
        watershedMap[0:75,81:144] = 2*ones((75,63))
        watershedMap[75:150,81:144] = 3*ones((75,63))
        watershedMap[62:75,70:81] = 4*ones((13,11))
        watershedMap[75:88,70:81] = 5*ones((13,11))
        self.watershedMap = watershedMap
        #plt.clf()
        #im1 = plt.imshow(transpose(self.watershedMap),origin='lower')
        #plt.colorbar(im1)
        #plt.savefig('watershedMap.png')
        #plt.cla()
        #plt.clf()
        
    def watershedQuadDot(self):
        watershedMap = zeros(self.effectivePotential.shape)
        watershedMap[43:60,43:70] = 2*ones((17,27)) #LQ1
        watershedMap[66:87,43:70] = 3*ones((21,27)) #RQ1
        watershedMap[66:87,80:105] = 4*ones((21,25)) #LQ2
        watershedMap[92:108,80:105] = 5*ones((16,25)) #RQ2

        self.watershedMap = watershedMap
        #plt.clf()
        #im1 = plt.imshow(transpose(self.watershedMap),origin='lower',aspect=self.Ly/self.Lx)
        #plt.colorbar(im1)
        #plt.savefig('watershedMap.png')
        #plt.cla()
        #plt.clf()

    
    def determineWatershedRegions(self):
        self.labelListRegions = []
        for i,label in enumerate(self.labelList):
            self.labelListRegions += [self.watershedMap[label[1][0],label[1][1]]]                    
    
    def computeDensities(self):
        self.densities = []
        self.detunings = []
        for region in self.labelListRegions:
            densityMap = self.lateralDensity*map(lambda x:x==region,self.watershedMap)
            potentialMap = self.effectivePotential*map(lambda x:x==region,self.watershedMap)
            self.detunings += [amin(potentialMap)]
            areaMap = map(lambda x:x!=0.0,densityMap)
            #print trapz(trapz(areaMap,self.yVec),self.xVec)
            self.densities += [trapz(trapz(densityMap,self.yVec),self.xVec)]            

#######  Computes the tunneling coefficient.
    def computeTunnelings(self,tunnelingPairs):
        self.Ts = []
        for tunnelPair in tunnelingPairs:
            potentialFunc = RectBivariateSpline(range(len(self.xVec)),range(len(self.yVec)),self.effectivePotential)
            potentialFunc = vectorize(potentialFunc)
            Lx = self.xVec[1]-self.xVec[0]
            Ly = self.yVec[1]-self.yVec[0]
            def bigT(points):
                startY = points[1]; startX = points[0]
                endY = points[3]; endX = points[2]
                ySlice = linspace(startY,endY,100)
                xSlice = linspace(startX,endX,100)
                rSize = sqrt(((ySlice[0]-ySlice[1])*Ly)**2+((xSlice[0]-xSlice[1])*Lx)**2)
                potentialSlice = potentialFunc(xSlice,ySlice)
                def zeroThreshold(x):
                    if x<0:
                        return 0.0
                    else:
                        return x
                potentialSlice2 = map(zeroThreshold,potentialSlice)
                if max(potentialSlice)<=0.0:
                    return 2.0
                else:
                    # Note: I am abs'ing the trapz since the values of rSlice might run backwards
                    expFactor = exp(-2*abs(trapz(sqrt(2*self.mt/self.hbar**2*array(potentialSlice2)*self.q),dx=rSize*1e-9)))
                    return expFactor

            startRegion = self.labelListRegions[tunnelPair[0]]
            endRegion = self.labelListRegions[tunnelPair[1]]

            def objectiveFunc(points):
                startY = points[1]; startX = points[0]
                endY = points[3]; endX = points[2]
                startBooleanMap = map(lambda x: (x==startRegion).astype(float), self.watershedMap)
                endBooleanMap = map(lambda x: (x==endRegion).astype(float), self.watershedMap)
                if startX >=len(self.xVec) or startX>=len(self.xVec) or endY>=len(self.yVec) or endY>=len(self.yVec):
                    return nan
                if (1.0-startBooleanMap[int(startX)][int(startY)])>7e-2 or (1.0-endBooleanMap[int(endX)][int(endY)])>7e-2 or potentialFunc(endX,endY)>0.0 or potentialFunc(endX,endY)>0.0 :
                    return 2.0
                else:
                    TVal = bigT(points)
                    if TVal == 2.0:
                        return 0.0
                    else:
                        return 1.0-TVal
            point0 = self.labelList[tunnelPair[0]][1]
            point1 = self.labelList[tunnelPair[1]][1]
            points = [point0[0],point0[1],point1[0],point1[1]]
            minObj = minimize(objectiveFunc,points,method='Nelder-Mead')
            [x0,y0,x1,y1] = minObj.x
            if minObj.fun==2:
                self.Ts += [0.0]
            else:
                self.Ts += [1-minObj.fun]

if __name__ == '__main__':
    filename = str(sys.argv[1])
    #
    # 2 qubit, 4 dot barrier device
    #
    #print "Tunnel"
    #guesses = [['leftDotQ1',[50,55]],['rightDotQ1',[75,55]],['leftDotQ2',[75,100]],['rightDotQ2',[100,100]]]
    guesses = [['leftDot',[65,75]],['rightDot',[80,75]],['leftLead',[50,110]],['rightLead',[90,110]]]
    
    tunnelingPairs = [[0,1],[0,2],[1,3]]
    
    startTime = time.time() 
    
    procDat = processData(filename)
    
    procDat.processPotential(guesses,tunnelingPairs)
    densities = procDat.densities
    Ts = procDat.Ts
    distinguishedRegions = len(set(procDat.labelListRegions))
    densityString = ' '.join(map(lambda x:str(x),densities))
    Tstring = ' '.join(map(lambda x:str(x),Ts))
    writeString = densityString+' '+Tstring+'\n'
    
    ##outputFile = open('Output2.txt','r')
    ##outputLines  = outputFile.readlines()
    ##outputFile.close()
    ##outputLines += [writeString]
    ##outputFile = open('Output2.txt','w')
    ##outputFile.writelines(outputLines)
    ##outputFile.close()
    #print 'Qubit 1: Detuning: '+str(procDat.detunings[0]-procDat.detunings[1])
    #print 'Qubit 2: Detuning: '+str(procDat.detunings[2]-procDat.detunings[3])
    sys.stdout.write(str(densities[0])+','+str(densities[1])+','+str(densities[2])+','+str(densities[3])+','+str(Ts[0])+','+str(Ts[1])+','+str(Ts[2]))
    '''sys.stdout.write('Qubit 1: Left Dot: ' + str(densities[0]))
    sys.stdout.write('Qubit 1: Right Dot: ' + str(densities[1]))
    sys.stdout.write('Qubit 2: Left Lead: ' + str(densities[2]))
    sys.stdout.write('Qubit 2: Right Lead: ' + str(densities[3]))
    sys.stdout.write('Between Dots: ' + '%e'%(Ts[0]))
    sys.stdout.write('Left Lead/Dot: ' + '%e'%(Ts[1]))
    sys.stdout.write('Right Lead/Dot: ' + '%e'%(Ts[2]))    '''           
    #print i
    #print densityString + Tstring