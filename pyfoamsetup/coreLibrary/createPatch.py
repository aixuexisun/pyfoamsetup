import numpy as np
from collections import OrderedDict
import os

import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self):
		super().__init__()

		self.FoamFile = OrderedDict()
		self.FoamFile['version']  = 2.0
		self.FoamFile['format']   ='ascii'
		self.FoamFile['class']    ='dictionary'
		self.FoamFile['location'] ='"system"'
		self.FoamFile['object']   ='createPatchDict'

		self.patchList = []

	def addPatch(self, name, patches, neighbourPatch):
		patch = OrderedDict()

		patch['name'] = name

		patch['patchInfo'] = OrderedDict()
		patch['patchInfo']['type'] = 'cyclicAMI'
		patch['patchInfo']['matchTolerance'] = 0.001
		patch['patchInfo']['neighbourPatch'] = neighbourPatch
		patch['patchInfo']['transform']      = 'noOrdering'

		patch['constructFrom'] = 'patches'
		patch['patches'] = '({})'.format(patches)

		self.patchList.append(patch)

	def addTranslationalPatch(self, name, patches, neighbourPatch, translationVector):
		patch = OrderedDict()

		patch['name'] = name

		patch['patchInfo'] = OrderedDict()
		patch['patchInfo']['type'] = 'cyclicAMI'
		patch['patchInfo']['matchTolerance'] = 0.001
		patch['patchInfo']['neighbourPatch'] = neighbourPatch
		patch['patchInfo']['transform']      = 'translational'
		patch['patchInfo']['separationVector'] = '({0} {1} {2})'.format(translationVector[0], translationVector[1], translationVector[2])

		patch['constructFrom'] = 'patches'
		patch['patches'] = '({})'.format(patches)

		self.patchList.append(patch)

	def write(self, folder, ending=''):
		f = open(folder + 'createPatchDict' + ending, 'w')

		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		f.write('pointSync false;\n\n')

		f.write('patches\n(\n')

		for i in range(len(self.patchList)):
			self.writeDict(f, self.patchList[i], 'patch', wrapper=True, writeDictName=False)


		f.write(');\n')

		f.close()
