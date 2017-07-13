from gtools import *

#fname = './catalogs/1130788624_source_array.sav'
fname = './catalogs/1130784064_source_array.sav'
#fname = './catalogs/1130781304_source_array.sav'
nside = 1024
var = 30000.0
cap = 100.0
makeGaussPlot(fname, nside, var, cap)
