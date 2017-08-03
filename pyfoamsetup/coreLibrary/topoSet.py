import numpy as np
from collections import OrderedDict

class boxToCell():
	def __init__(self, name, pMin, pMax):
		self.pMin = pMin
		self.pMax = pMax

		self.settings = OrderedDict()

		self.settings['name']   = name
		self.settings['type']   = 'cellSet'
		self.settings['action'] = 'new'
		self.settings['source'] = 'boxToCell'

	def generateString(self):
		string = '\t{\n'

		for k, v in self.settings.items():
			string += '\t\t{0} {1};\n'.format(k, v)

		numberString = '({:.6f} {:.6f} {:.6f}) ({:.6f} {:.6f} {:.6f})'.format(self.pMin[0], self.pMin[1], self.pMin[2], self.pMax[0], self.pMax[1], self.pMax[2])


		string += '\t\tsourceInfo\n\t\t{\n\t\t\tbox '+numberString+';\n\t\t}\n'
		string += '\t}\n\n'

		return string

class cylinderToCell():
	def __init__(self, name, p1, p2, radius):
		self.p1     = p1
		self.p2     = p2
		self.radius = radius

		self.settings = OrderedDict()

		self.settings['name']   = name
		self.settings['type']   = 'cellSet'
		self.settings['action'] = 'new'
		self.settings['source'] = 'cylinderToCell'

	def generateString(self):
		string = '\t{\n'

		for k, v in self.settings.items():
			string += '\t\t{0} {1};\n'.format(k, v)

		string += '\t\tsourceInfo\n\t\t{\n'

		string += '\t\t\tp1 ({:.6f} {:.6f} {:.6f});\n'.format(self.p1[0], self.p1[1], self.p1[2])
		string += '\t\t\tp2 ({:.6f} {:.6f} {:.6f});\n'.format(self.p2[0], self.p2[1], self.p2[2])
		string += '\t\t\tradius {:.6f};\n'.format(self.radius)

		string += '\t\t}\n'
		
		string += '\t}\n\n'

		return string


def header():
	return 'FoamFile\n{\n\tversion    2.0;\n\tformat    ascii;\n\tclass    dictionary;\n\tobject    topoSetDict;\n}\n\n'