from gtools import *

names = ['./catalogs/1130788624_source_array.sav']#, './catalogs/1130784064_source_array.sav', './catalogs/1130781304_source_array.sav']
nside = 1024
variances = [0.1]#[0.3, 0.1, 0.03, 0.01]
caps = [100.0, 50.0, 20.0]
for fname in names:
    for var in variances:
        for cap in caps:
            makeGaussPlot(fname, nside, var, cap)
