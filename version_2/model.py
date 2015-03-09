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
        self.session.putvalue('MSCRIPT',"""model = mphload('"""+filename+"""');""")
        self.session.run('eval(MSCRIPT)')
    
    def changeHeight(self,feature,height,geom='geom1',mesh='mesh1'):
        self.session.run("""height = """+str(height))
        mscript = """model.geom('"""+geom+"""').feature('"""+feature+"""').set('quickz',height)
while true
    try
        model.mesh('"""+mesh+"""').run()
        break
    catch e
        height = height+0.00001
        model.geom('"""+geom+"""').feature('"""+feature+"""').set('quickz',height)
    end
end"""
        self.session.putvalue('MSCRIPT',mscript)
        self.session.run('eval(MSCRIPT)')
    
    def getVoltage(self, feature):
        self.session.run("""V = str2double(model.physics('es').feature('"""+feature+"""').getString('V0'))""")
        return self.session.getvalue('V')
    
    def setVoltage(self, feature, value):
        self.session.run("""model.physics('es').feature('"""+feature+"""').set('V0',"""+str(value)+""")""")
    
    def runSimulation(self,study='std1'):
        self.session.run("""model.study('"""+study+"""').run()""")
    
    def exportData(self,filename,data='data1'):
        mscript = """model.result.export('"""+data+"""').set('filename','"""+filename+"""');
model.result.export('"""+data+"""').run;"""
        self.session.putvalue('MSCRIPT',mscript)
        self.session.run('eval(MSCRIPT)')
