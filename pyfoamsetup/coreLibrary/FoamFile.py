import numpy as np
from collections import OrderedDict
import os

class Dict():
	def __init__(self):
		self.libFolder = os.path.dirname(os.path.realpath(__file__))
		f = open(self.libFolder + '/fileHeader.txt')
		
		self.header = f.read()

		f.close()

	def findMaxLength(self, dictVar):
		maxLen = 0
		for k, v in dictVar.items():
			if len(k) > maxLen:
				maxLen = len(k)

		return maxLen

	def writeDictVar(self, f, dictVar, nrTabs=1):
		maxLen = self.findMaxLength(dictVar)

		for k, v in dictVar.items():
			if type(v) == type(dictVar):
				f.write(nrTabs*'\t' + k + '\n' + nrTabs*'\t' + '{\n')
				self.writeDictVar(f, v, nrTabs=nrTabs+1)
				f.write(nrTabs*'\t' + '}' + '\n')
			else:
				if k =='':
					f.write(nrTabs*'\t' + '{0};\n'.format(v))
				else:
					nSpaces = maxLen - len(k) + 4
					f.write(nrTabs*'\t' + '{0}'.format(k) + nSpaces*' ' + '{0};\n'.format(v))

	def writeDict(self, f, dictVar, dictName, wrapper=True, writeDictName=True, nrTabs=1):
		if wrapper:
			if writeDictName:
				f.write(dictName)
			f.write('\n{\n')
			
			tabs = '\t'
		else:
			tabs = ''

		self.writeDictVar(f, dictVar, nrTabs=nrTabs)

		if wrapper:
			f.write('}\n\n')

