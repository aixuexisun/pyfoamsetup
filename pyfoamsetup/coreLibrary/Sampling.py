import numpy as np
from collections import OrderedDict
import os
import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self):
		super().__init__()

		self.nSamples = 0
		self.sampleList = []
		self.sampleNameList = []
		self.fileName = 'samples'

	def samplePressureOnPatch(self, name, patchName, writeInterval = 10, pName='p'):
		self.nSamples += 1
		self.sampleNameList.append(name)

		sample = OrderedDict()

		sample['type']                = 'surfaces'
		sample['functionObjectLibs']  = '("libsampling.so")'
		sample['writeControl']        = 'timeStep'
		sample['writeInterval']       = writeInterval
		sample['surfaceFormat']       = 'raw'
		sample['fields']              = '\n\t(\n\t\t{0}\n\t)'.format(pName)
		sample['interpolationScheme'] = 'cellPoint'
		sample['surfaces']            = '\n\t(\n\t\twalls\n\t\t{\n\t\t\ttype        patch;\n\t\t\tpatches     ('+patchName+');\n\t\t\ttriangulate false;\n\t\t}\n\t)'

		self.sampleList.append(sample)

	def write(self, folder):
		if self.nSamples > 0:
			filePath = folder + self.fileName

			f = open(filePath, 'w')

			f.write(self.header)
			for i in range(self.nSamples):
				self.writeDict(f, self.sampleList[i], self.sampleNameList[i])

			f.close()

