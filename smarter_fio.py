# smarter_fio.py

import sys, os
import numpy as np 
filename = '1130781304_run1_catalog.txt';
debugging = 0

def dprint(str):
	if debugging==1:
		print(str)
		
class Source(object):
	"""
	
	Document!!!
	
	"""

	def __init__(self, RA, dec, flux, extended, elements, elementDataArray):
		self.RA = RA
		self.dec = dec
		self.flux = flux
		self.ext = extended
		if extended:
			self.elements = elements
			self.elRA = elementDataArray[:,0]
			self.elDec = elementDataArray[:,1]
			self.elFlux = elementDataArray[:,2]
	
	def Gaussian(self, raArray, decArray, cap):
		
	
	