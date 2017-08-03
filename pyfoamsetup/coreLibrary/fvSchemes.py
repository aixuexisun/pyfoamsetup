import numpy as np
from collections import OrderedDict
import os
import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self, solver, turbulenceModel, schemeType='robustAndAccurate', timeStepScheme='Euler', turbulenceType='RAS'):
		super().__init__()

		self.FoamFile = OrderedDict()
		self.FoamFile['version']  = 2.0
		self.FoamFile['format']   ='ascii'
		self.FoamFile['class']    ='dictionary'
		self.FoamFile['location'] ='"system"'
		self.FoamFile['object']   ='fvSchemes'

		self.kOmegaList             = ['kOmega', 'kOmegaSST']
		self.kOmegaLESList          = ['kOmegaSSTDES', 'kOmegaSSTDDES', 'kOmegaSSTIDDES']
		self.kEpsilonList           = ['kEpsilon', 'realizableKE']
		self.spalartAllmarasList    = ['SpalartAllmaras']
		self.spalartAllmarasLESList = ['SpalartAllmarasDES', 'SpalartAllmarasDDES', 'SpalartAllmarasIDDES']
		self.LESList                = ['Smagorinsky', 'WALE']

		# ------------------------ Time integration ----------------------------------------------------
		self.ddtSchemes = OrderedDict()
		if solver == 'simpleFoam':
			self.ddtSchemes['default'] = 'steadyState'
		else:
			if timeStepScheme == 'CrankNicolson':
				self.ddtSchemes['default'] = 'CrankNicolson 0.75'
			elif timeStepScheme == 'CrankNicolson 0.9':
				self.ddtSchemes['default'] = 'CrankNicolson 0.9'
			elif timeStepScheme == 'backward':
				self.ddtSchemes['default'] = 'backward'
			elif timeStepScheme == 'Euler':
				self.ddtSchemes['default'] = 'Euler'
			elif timeStepScheme == 'localEuler':
				self.ddtSchemes['default'] = 'localEuler'
			else:
				self.ddtSchemes['default'] = 'Euler'

		# ------------------------ Grad schemes --------------------------------------------------------
		self.gradSchemes = OrderedDict()
		self.gradSchemes['default'] = 'Gauss linear'

		if schemeType == 'robustAndAccurate':
			self.gradSchemes['limitedGrad'] = 'cellLimited Gauss linear 1'

		# ------------------------ Div schemes --------------------------------------------------------
		self.divSchemes = OrderedDict()

		if solver == 'interFoam' or solver == 'interDyMFoam':		
			divPhiU_name = 'div(rhoPhi,U)'
		else:
			divPhiU_name = 'div(phi,U)'

		if solver == 'simpleFoam':
			preString = 'bounded '
		else:
			preString = ''

		if schemeType == 'accurate':
			self.divSchemes[divPhiU_name] = preString + 'Gauss LUST grad(U)'
		elif schemeType == 'robustAndAccurate':
			if solver == 'interFoam' or solver == 'interDyMFoam':
				self.divSchemes[divPhiU_name] = preString + 'Gauss linearUpwindV grad(U)'
			else:
				self.divSchemes[divPhiU_name] = preString + 'Gauss linearUpwindV grad(U)'
		else:
			self.divSchemes[divPhiU_name] = preString + 'Gauss upwind'

		if solver == 'interFoam' or solver == 'interDyMFoam':
			self.divSchemes['div(phi,alpha)']   = 'Gauss vanLeer'
			self.divSchemes['div(phirb,alpha)'] = 'Gauss linear'

		turbulenceVariables = []
		if turbulenceModel in self.kOmegaList:
			turbulenceVariables = ['k', 'omega']
		elif turbulenceModel in self.kOmegaLESList:
			turbulenceVariables = ['k', 'omega']
		elif turbulenceModel in self.kEpsilonList:
			turbulenceVariables = ['k', 'epsilon']
		elif turbulenceModel in self.spalartAllmarasList:
			turbulenceVariables = ['nuTilda']
		elif turbulenceModel in self.spalartAllmarasLESList:
			turbulenceVariables = ['nuTilda']
		elif turbulenceModel in self.LESList:
			turbulenceVariables = ['k', 'B', 'nuTilda']
		else:
			raise ValueError('Unknown turbulence model in the fvScheme file!')

		if schemeType == 'accurate':
			for i in range(len(turbulenceVariables)):
				self.divSchemes['div(phi,{})'.format(turbulenceVariables[i])] = preString + 'Gauss linearUpwind  default'
		elif schemeType == 'robustAndAccurate':
			if solver == 'interFoam' or solver == 'interDyMFoam':
				for i in range(len(turbulenceVariables)):
					self.divSchemes['div(phi,{})'.format(turbulenceVariables[i])] = preString + 'Gauss linearUpwind limitedGrad'
			else:
				for i in range(len(turbulenceVariables)):
					self.divSchemes['div(phi,{})'.format(turbulenceVariables[i])] = preString + 'Gauss linearUpwind limitedGrad'
		else:
			for i in range(len(turbulenceVariables)):
				self.divSchemes['div(phi,{})'.format(turbulenceVariables[i])] = preString + 'Gauss upwind'

		if solver == 'interFoam' or solver == 'interDyMFoam':
			self.divSchemes['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'
		else:
			self.divSchemes['div((nuEff*dev2(T(grad(U)))))'] = 'Gauss linear'

		# ------------------------ Laplacian schemes --------------------------------------------------------
		self.laplacianSchemes = OrderedDict()
		self.laplacianSchemes['default'] = 'Gauss linear corrected'
		#self.laplacianSchemes['laplacian(diffusivity,cellDisplacement)'] = 'Gauss linear corrected'

		# ------------------------ Interpolation schemes --------------------------------------------------------
		self.interpolationSchemes = OrderedDict()
		self.interpolationSchemes['default'] = 'linear'

		# ------------------------ snGrad schemes --------------------------------------------------------
		self.snGradSchemes = OrderedDict()
		self.snGradSchemes['default'] = 'corrected'

		self.wallDist = OrderedDict()
		self.wallDist['method']       = 'meshWave'
		self.wallDist['correctWalls'] = 'true'

	def write(self, folder, ending=''):
		f = open(folder + 'fvSchemes' + ending, 'w')

		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		self.writeDict(f, self.ddtSchemes, 'ddtSchemes')
		self.writeDict(f, self.gradSchemes, 'gradSchemes')
		self.writeDict(f, self.divSchemes, 'divSchemes')
		self.writeDict(f, self.laplacianSchemes, 'laplacianSchemes')
		self.writeDict(f, self.interpolationSchemes, 'interpolationSchemes')
		self.writeDict(f, self.snGradSchemes, 'snGradSchemes')
		self.writeDict(f, self.wallDist, 'wallDist')

		f.close()
