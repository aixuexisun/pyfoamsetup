from collections import OrderedDict
import os

import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self, nCPUs, method):
		super().__init__()

		self.FoamFile = OrderedDict()
		self.FoamFile['version'] = 2.0
		self.FoamFile['format']  ='ascii'
		self.FoamFile['class']   ='dictionary'
		self.FoamFile['object']  ='snappyHexMeshDict'

		self.nCPUs = nCPUs
		self.method = method

		if self.method == 'simple':
			if self.nCPUs == 16:
				self.simpleCoeffs = [4, 4, 1]
			elif self.nCPUs == 12:
				self.simpleCoeffs = [4, 3, 1]
			elif self.nCPUs == 8:
				self.simpleCoeffs = [4, 2, 1]
			elif self.nCPUs == 4:
				self.simpleCoeffs = [2, 2, 1]
			elif self.nCPUs == 2:
				self.simpleCoeffs = [2, 1, 1]
			else:
				self.simpleCoeffs = [self.nCPUs, 1, 1]

	def write(self, folder, ending=''):
		filePath = folder + 'decomposeParDict' + ending

		f = open(filePath, 'w')

		# Write FoamFile header
		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		f.write('numberOfSubdomains {:.0f};\n\n'.format(self.nCPUs))
		f.write('method             {};\n\n'.format(self.method))

		if self.method == 'scotch':
			f.write('scotchCoeffs\n{\n}\n')
		elif self.method == 'simple':
			f.write('simpleCoeffs\n{\n')
			f.write('\tn     ({:.0f} {:.0f} {:.0f});\n'.format(self.simpleCoeffs[0], self.simpleCoeffs[1], self.simpleCoeffs[2]))
			f.write('\tdelta 0.001;\n')
			f.write('}\n')

		f.close()

