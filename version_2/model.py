from numpy import *
import pymatlab

class model:
    def __init__(self):
        self.session = pymatlab.session_factory()
        self.session.run('mphstart(2036)')
        self.session.run('import com.comsol.model.*')
        self.session.run('import com.comsol.model.util.*')
    
    def loadModelFile(self,filename):
        self.session.run('model = mphload('+filename+');')
    
    def changeHeight(self,feature,height):
        self.session.putvalue('height',height)
        mscript = """model.geom('geom1').feature('"""+feature+"""').set('quickz',height)
        while true
            try
                model.mesh('mesh1').run()
                break
            catch e
                height = height+0.00001
                model.geom('geom1').feature('wp4').set('quickz',height)
            end
        end"""
        self.session.putvalue('MSCRIPT',mscript)
        self.session.run('eval(MSCRIPT)')