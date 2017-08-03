import numpy as np
import os
from collections import OrderedDict
import pyfoamsetup.coreLibrary.FileHandling as FileHandling
import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self, U, L, nu, model, modelType):
		super().__init__()

		self.U         = U
		self.L         = L
		self.nu        = nu
		self.model     = model
		self.modelType = modelType  
		self.I         = 0.01
		self.viscRatio = 10
		self.beta      = 0.075

		self.FoamFile = OrderedDict()
		self.FoamFile['version']  = 2.0
		self.FoamFile['format']   ='ascii'
		self.FoamFile['class']    ='dictionary'
		self.FoamFile['location'] ='"constant"'
		self.FoamFile['object']   ='turbulenceProperties'

		if self.modelType == 'RAS':
			self.modelTypeDict = OrderedDict()
			self.modelTypeDict['RASModel']    = self.model
			self.modelTypeDict['turbulence']  = 'on'
			self.modelTypeDict['printCoeffs'] = 'on'
		elif self.modelType == 'LES':
			self.modelTypeDict = OrderedDict()
			self.modelTypeDict['LESModel']    = self.model
			self.modelTypeDict['turbulence']  = 'on'
			self.modelTypeDict['delta']       = 'cubeRootVol'
			
			self.modelTypeDict['printCoeffs'] = 'on'

			self.modelTypeDict['cubeRootVolCoeffs'] = OrderedDict()
			self.modelTypeDict['cubeRootVolCoeffs']['deltaCoeff'] = 1

			self.modelTypeDict['PrandtlCoeffs'] = OrderedDict()
			self.modelTypeDict['PrandtlCoeffs']['delta'] = 'cubeRootVol'
			self.modelTypeDict['PrandtlCoeffs']['cubeRootVolCoeffs'] = OrderedDict()
			self.modelTypeDict['PrandtlCoeffs']['cubeRootVolCoeffs']['deltaCoeff'] = 1
			self.modelTypeDict['PrandtlCoeffs']['smoothCoeffs'] = OrderedDict()
			self.modelTypeDict['PrandtlCoeffs']['smoothCoeffs']['delta'] = 'cubeRootVol'
			self.modelTypeDict['PrandtlCoeffs']['smoothCoeffs']['cubeRootVolCoeffs'] = OrderedDict()
			self.modelTypeDict['PrandtlCoeffs']['smoothCoeffs']['cubeRootVolCoeffs']['deltaCoeff'] = 1
			self.modelTypeDict['PrandtlCoeffs']['smoothCoeffs']['maxDeltaRatio'] = 1.1
			self.modelTypeDict['PrandtlCoeffs']['Cdelta'] = 0.158

			self.modelTypeDict['vanDriestCoeffs'] = OrderedDict()
			self.modelTypeDict['vanDriestCoeffs']['delta'] = 'cubeRootVol'
			self.modelTypeDict['vanDriestCoeffs']['cubeRootVolCoeffs'] = OrderedDict()
			self.modelTypeDict['vanDriestCoeffs']['cubeRootVolCoeffs']['deltaCoeff'] = 1
			self.modelTypeDict['vanDriestCoeffs']['smoothCoeffs'] = OrderedDict()
			self.modelTypeDict['vanDriestCoeffs']['smoothCoeffs']['delta'] = 'cubeRootVol'
			self.modelTypeDict['vanDriestCoeffs']['smoothCoeffs']['cubeRootVolCoeffs'] = OrderedDict()
			self.modelTypeDict['vanDriestCoeffs']['smoothCoeffs']['cubeRootVolCoeffs']['deltaCoeff'] = 1
			self.modelTypeDict['vanDriestCoeffs']['smoothCoeffs']['maxDeltaRatio'] = 1.1
			self.modelTypeDict['vanDriestCoeffs']['Adelta'] = 26
			self.modelTypeDict['vanDriestCoeffs']['Cdelta'] = 0.158

			self.modelTypeDict['smoothCoeffs'] = OrderedDict()
			self.modelTypeDict['smoothCoeffs']['delta'] = 'cubeRootVol'
			self.modelTypeDict['smoothCoeffs']['cubeRootVolCoeffs'] = OrderedDict()
			self.modelTypeDict['smoothCoeffs']['cubeRootVolCoeffs']['deltaCoeff'] = 1
			self.modelTypeDict['smoothCoeffs']['maxDeltaRatio'] = 1.1

	def calculateInitialValues(self):
		self.k       = (3/2)*(self.U*self.I)**2
		self.omega   = 0.09*self.k/(self.nu*self.viscRatio)
		self.epsilon = 0.09*self.k**2/(self.nu*self.viscRatio)
		self.nut     = self.viscRatio*self.nu 
		self.Re      = self.L*self.U/self.nu
		self.ybl     = 0.382/self.Re**(1/5)

	def writeInitialFile(self, caseFolder):
		includeFolder = caseFolder + '/0/include'
		os.makedirs(includeFolder)

		f = open(includeFolder+'/initialConditions', 'w')
		f.write('flowVelocity    ({:.6f} 0 0);\n'.format(self.U))
		f.write('Umean            {:.6f};\n'.format(self.U))
		f.write('turbulentKE      {:.6e};\n'.format(self.k))
		f.write('turbulentOmega   {:.6e};\n'.format(self.omega))
		f.write('turbulentEpsilon {:.6e};\n'.format(self.epsilon))
		f.write('turbulentNut     {:.6e};\n'.format(self.nut))
		f.write('#inputMode       merge;')
		f.close()

	def write(self, folder, ending=''):
		f = open(folder + 'turbulenceProperties' + ending, 'w')

		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		# turbulence type
		f.write('simulationType {};\n\n'.format(self.modelType))

		self.writeDict(f, self.modelTypeDict, self.modelType)

		f.close()