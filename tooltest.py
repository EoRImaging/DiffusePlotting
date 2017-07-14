from gtools import *
import resource
from datetime import datetime
#_, hard = resource.getrlimit(resource.RLIMIT_DATA)
#resource.setrlimit(resource.RLIMIT_DATA, (2**32,hard))
#print('Max memory usage restricted to 4 GB')

startTime = datetime.now()

names = ['./catalogs/1130788624_source_array.sav']#, './catalogs/1130784064_source_array.sav', './catalogs/1130781304_source_array.sav']
nside = 4096
variances = [0.1]#, 0.1, 0.03, 0.01]
caps = [5.0, 4.0, 3.0, 2.0, 1.0]
for fname in names:
    for var in variances:
        for cap in caps:
            makeGaussPlot(fname, nside, var, cap, loadData = True, colorbarLevels=10)
print('tooltest.py executed in '+str(datetime.now()-startTime)+' seconds.')
