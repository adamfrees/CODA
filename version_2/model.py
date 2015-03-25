from numpy import *
import pymatlab
import os
import time

class model:
    def __init__(self):
        os.system('comsol server &')
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
    
    def runSimulation(self,study='std1'):
        self.session.run("""model.study('{STUDY}').run()""".format(STUDY=study))
    
    def exportData(self,filename,data='data1'):
        args = {'FILENAME' : filename, 'DATA' : data}
        mscript = """model.result.export('{DATA}').set('filename','{FILENAME}');
model.result.export('{DATA}').run;""".format(**args)
        self.session.putvalue('MSCRIPT',mscript)
        self.session.run('eval(MSCRIPT)')
