import numpy as np
from collections import OrderedDict

import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self, solver, localTimeStep=False):
		super().__init__()

		self.solver = solver

		self.FoamFile = OrderedDict()
		self.FoamFile['version'] = 2.0
		self.FoamFile['format']  ='ascii'
		self.FoamFile['class']   ='dictionary'
		self.FoamFile['object']  ='controlDict'

		self.Settings = OrderedDict()
		self.Settings['application']       = self.solver
		self.Settings['startFrom']         = 'latestTime'
		self.Settings['startTime']         = 0
		self.Settings['stopAt']            = 'endTime'
		self.Settings['endTime']           = 100
		self.Settings['deltaT']            = 1

		if (self.solver == 'pimpleFoam' or self.solver == 'interFoam' or self.solver == 'interDyMFoam' or self.solver=='pimpleDyMFoam') and not(localTimeStep):
			self.Settings['writeControl']  = 'adjustableRunTime'
		else:
			self.Settings['writeControl']  = 'timeStep'

		self.Settings['writeInterval']     = 100
		self.Settings['purgeWrite']        = 2
		self.Settings['writeFormat']       = 'ascii'
		self.Settings['writePrecision']    = 9
		self.Settings['writeCompression']  = 'uncompressed'
		self.Settings['timeFormat']        = 'general'
		self.Settings['timePrecision']     = 9
		self.Settings['runTimeModifiable'] = 'yes' 

		if (self.solver == 'pimpleFoam' or self.solver == 'interFoam' or self.solver == 'interDyMFoam' or self.solver == 'pimpleDyMFoam') and not(localTimeStep):
			self.Settings['adjustTimeStep'] = 'no'
			self.Settings['maxCo']          = 1
			self.Settings['maxDeltaT']      = 0.1

		if (self.solver == 'interFoam' or self.solver == 'interDyMFoam') and not(localTimeStep):
			self.Settings['maxAlphaCo'] = 1

		self.FunctionList = []

	def write(self, folder):
		filePath = folder + 'controlDict'

		f = open(filePath, 'w')

		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')
		self.writeDict(f, self.Settings, 'settings', wrapper=False, nrTabs=0)

		f.write('\n')

		f.write('functions\n{\n')
		if len(self.FunctionList) > 0:
			for item in self.FunctionList:
				f.write(item)
				f.write('\n')
		f.write('}\n')

		f.close()