from numpy import *
import pymatlab
import os
import time

class model:
    def __init__(self,limitCoresTo = None):
        if limitCoresTo is None:
            os.system('comsol server &')
        else:
            os.system('comsol server -np '+str(int(limitCoresTo))+' &')
        #time.sleep(1)
        self.session = pymatlab.session_factory()
        self.session.run('mphstart(2036)')
        self.session.run('import com.comsol.model.*')
        self.session.run('import com.comsol.model.util.*')
    
    def loadModelFile(self,filename):
        self.session.putvalue('MSCRIPT',"""model = mphload('{FILENAME}');""".format(FILENAME=filename))
        self.session.run('eval(MSCRIPT)')
    
    def changeHeight(self,feature,height,geom='geom1',mesh='mesh1'):
        args = {'FEATURE' : feature, 'HEIGHT': height, 'GEOM': geom, 'MESH': mesh}
        self.session.run("""height = {HEIGHT}""".format(**args))
        mscript = """model.geom('{GEOM}').feature('{FEATURE}').set('quickz',height)
while true
    try
        model.mesh('{MESH}').run()
        break
    catch e
        height = height+0.00001
        model.geom('{GEOM}').feature('{FEATURE}').set('quickz',height)
    end
end""".format(**args)
        self.session.putvalue('MSCRIPT',mscript)
        self.session.run('eval(MSCRIPT)')
    
    def getVoltage(self, feature):
        self.session.run("""V = str2double(model.physics('es').feature('{FEATURE}').getString('V0'))""".format(FEATURE=feature))
        return self.session.getvalue('V')
    
    def setVoltage(self, feature, value):
        self.session.run("""model.physics('es').feature('{FEATURE}').set('V0',{VALUE})""".format(FEATURE=feature,VALUE=str(value)))
    
    def getParam(self, feature):
        self.session.run("""V = str2double(model.param.get('{FEATURE}'))""".format(FEATURE=feature))
        return self.session.getvalue('V')
    
    def setParam(self, feature, value):
        self.session.run("""model.param.set('{FEATURE}',{VALUE})""".format(FEATURE=feature,VALUE=str(value)))
        
    def runSimulation(self,study='std1'):
        self.session.run("""model.study('{STUDY}').run()""".format(STUDY=study))
    
    def exportData(self,filename,data='data1',solnum=None):
        args = {'FILENAME' : filename, 'DATA' : data, 'SOLNUM' : solnum}
        if solnum is None:
            mscript = """model.result.export('{DATA}').set('filename','{FILENAME}');
model.result.export('{DATA}').run;""".format(**args)
        else:
            mscript = """model.result.export('{DATA}').set('solnum',{SOLNUM});
model.result.export('{DATA}').set('filename','{FILENAME}');
model.result.export('{DATA}').run;""".format(**args)
        self.session.putvalue('MSCRIPT',mscript)
        self.session.run('eval(MSCRIPT)')
    
    def initializeSweep(self,gateNames,paramMatrix,solution='sol4',solver='s1',sweep='p1'):
        gateNames = "gates = "+repr(gateNames).replace('[','{').replace(']','}').replace(',','')
        self.session.putvalue('params',paramMatrix)
        args = {'SOLUTION' : solution, 'SOLVER' : solver, 'SWEEP' : sweep}
        mscript = gateNames + """
model.sol('{SOLUTION}').feature('{SOLVER}').feature('{SWEEP}').set('pname',gates);
model.sol('{SOLUTION}').feature('{SOLVER}').feature('{SWEEP}').set('plistarr',params);""".format(**args)
        self.session.putvalue('MSCRIPT',mscript)
        self.session.run('eval(MSCRIPT)')
