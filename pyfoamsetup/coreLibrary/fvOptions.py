import numpy as np
from collections import OrderedDict

class actuationDisk():
	def __init__(self, name, Ct, Cp, diskArea, diskDir, upstreamPoint, cellSet):
		self.settings = OrderedDict()

		self.name = name

		self.settings['selectionMode'] = 'cellSet'
		self.settings['cellSet']       = cellSet
		self.settings['fieldNames']    = '(U)'
		self.settings['Ct']            = Ct
		self.settings['Cp']            = Cp
		self.settings['diskArea']      = diskArea
		self.settings['diskDir']       = diskDir
		self.settings['upstreamPoint'] = upstreamPoint

	def generateString(self):
		string = self.name +'\n{\n'
		string += '\ttype      actuationDiskSource;\n\tactive    on;\n\n'
		string += '\tactuationDiskSourceCoeffs\n\t{\n'

		for k, v in self.settings.items():
			string += '\t\t{0} {1};\n'.format(k, v)

		string += '\t}\n}\n'

		return string

class goldsteinActuationDisk():
	def __init__(self, name, Thrust, Torque, p1, p2, Rp, Rh, cellSet):
		self.settings = OrderedDict()

		self.name = name

		self.settings['selectionMode'] = 'cellSet'
		self.settings['cellSet']       = cellSet
		self.settings['fieldNames']    = '(U)'
		self.settings['Thrust']        = Thrust
		self.settings['Torque']        = Torque
		self.settings['p1']            = '({:.6f} {:.6f} {:.6f})'.format(p1[0], p1[1], p1[2])
		self.settings['p2']            = '({:.6f} {:.6f} {:.6f})'.format(p2[0], p2[1], p2[2])
		self.settings['Rp']            = Rp
		self.settings['Rh']            = Rh

	def generateString(self):
		string = self.name +'\n{\n'
		string += '\ttype      goldsteinActuationDiskSource;\n\tactive    on;\n\n'
		string += '\tgoldsteinActuationDiskSourceCoeffs\n\t{\n'

		for k, v in self.settings.items():
			string += '\t\t{0} {1};\n'.format(k, v)

		string += '\t}\n}\n'

		return string

def header():
	return 'FoamFile\n{\n\tversion    2.0;\n\tformat    ascii;\n\tclass    dictionary;\n\tlocation    "constant";\n\tobject    fvOptions;\n}\n\n'

