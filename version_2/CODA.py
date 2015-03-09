from numpy import *
import pymatlab
from model import *

dDot = model()
dDot.loadModelFile('/Users/adamfrees/Dropbox/_UW/CODA/Var_Height/Data')
dDot.changeHeight('wp4',0.08)