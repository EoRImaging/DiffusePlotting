# radcos_fio.py

import sys, os
import numpy as np 
filename = '1130781304_run1_catalog.txt';

data = np.zeros((2,3))
with open(filename) as f:
	next(f)
	for line in f:
		linedata = f.read()
		linedata = linedata.split(',');
		# print(linedata)
		extended_flag = linedata[3]
		if extended_flag == 0:
			continue;
			# TBI
		num_elements = linedata[4]
		linedata = np.asarray(linedata[5:])
		linedata = linedata.reshape(-1,3)
		print(linedata.shape)
		print(data.shape)
		data = np.concatenate((data,linedata),axis=0)
data = data[2:,:]
print(data)
data.dtype
