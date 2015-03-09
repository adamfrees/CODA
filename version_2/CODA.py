from numpy import *
import pymatlab
from model import *

order = [3,1,5,14,15,6,7,8,9,19,10,11,16,12,17]

dDot = model()
dDot.loadModelFile('/home/frees/Dropbox/_UW/Scale-up/Comsol_test_1/var_height.mph')
dDot.changeHeight('wp4',0.08)
#for i in order:
#    print dDot.getVoltage('pot'+str(i))
i=3
v = dDot.getVoltage('pot'+str(i))
print v
dDot.setVoltage('pot'+str(i),v+0.005)
print dDot.getVoltage('pot'+str(i))
