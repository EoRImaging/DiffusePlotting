# radcos_fio.py

import sys, os
import numpy as np 
filename = '1130781304_run1_catalog.txt';
debugging = 0

def dprint(str):
	if debugging==1:
		print(str)
		

data = np.zeros((2,3)) # Hacky concatenation target
print(data.dtype)
with open(filename) as f:
	next(f) # Skip header
	i=1;
	# Iterate over lines in file; each line corresponds to a source
	for line in f:
		i = i+1;
		dprint('Line: ' + str(i))
		linedata = line[:-2].split(','); # Separate by commas
		# print(linedata)
		#if i>=2665:
		#	print(linedata)
		
		extended_flag = linedata[3] # Check if extended source
		if extended_flag == 0: # If not, drop single line onto final array and move on to the next
			linedata = np.asarray(linedata[0:3],dtype=np.float32)
			data = np.concatenate((data,linedata),axis=0)
			continue;
		# If the source is extended, take all data following the
		# header and append onto final array
		linedata = np.asarray(linedata[5:],dtype=np.float32)
		linedata = linedata.reshape(-1,3)
		dprint(linedata.shape)
		dprint(data.shape)
		data = np.concatenate((data,linedata),axis=0)
data = data[2:,:] # Trim off the hacky original array
print(data)
print(data.dtype)
