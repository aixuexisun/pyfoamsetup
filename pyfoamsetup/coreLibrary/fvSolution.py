import numpy as np
from collections import OrderedDict
import os
import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self, solver, turbulenceModel, version=3, localTimeStep=False):
		super().__init__()

		self.solver        = solver
		self.version       = version
		self.localTimeStep = localTimeStep

		self.kOmegaList             = ['kOmega', 'kOmegaSST']
		self.kOmegaLESList          = ['kOmegaSSTDES', 'kOmegaSSTDDES', 'kOmegaSSTIDDES']
		self.kEpsilonList           = ['kEpsilon', 'realizableKE']
		self.spalartAllmarasList    = ['SpalartAllmaras']
		self.spalartAllmarasLESList = ['SpalartAllmarasDES', 'SpalartAllmarasDDES', 'SpalartAllmarasIDDES']
		self.LESList                = ['Smagorinsky', 'WALE']

		if self.solver == 'pimpleFoam' or self.solver == 'pimpleDyMFoam' or self.solver == 'interFoam' or self.solver == 'interDyMFoam':
			self.solverName = 'PIMPLE'
		elif self.solver == 'simpleFoam':
			self.solverName = 'SIMPLE'
		else:
			raise ValueError('Unknown solver in the fvSolution file!')

		self.turbulenceVariables = []
		if turbulenceModel in self.kOmegaList:
			self.turbulenceVariables = ['k', 'omega']
		elif turbulenceModel in self.kOmegaLESList:
			self.turbulenceVariables = ['k', 'omega', 'B']
		elif turbulenceModel in self.kEpsilonList:
			self.turbulenceVariables = ['k', 'epsilon']
		elif turbulenceModel in self.spalartAllmarasList:
			self.turbulenceVariables = ['nuTilda']
		elif turbulenceModel in self.spalartAllmarasLESList:
			self.turbulenceVariables = ['nuTilda', 'B']
		elif turbulenceModel in self.LESList:
			turbulenceVariables = ['k', 'B', 'nuTilda']
		else:
			raise ValueError('Unknown turbulence model in the fvSolution file!')

		self.velocitySolverKey = '"(U'
		for i in range(len(self.turbulenceVariables)):
			self.velocitySolverKey += '|{}'.format(self.turbulenceVariables[i])
		self.velocitySolverKey += ')"'

		if self.solverName == 'PIMPLE':
			self.nOuterCorrectors         = 3
			self.nCorrectors              = 2
			self.nNonOrthogonalCorrectors = 0
		elif self.solverName == 'SIMPLE':
			self.nNonOrthogonalCorrectors = 1
			self.nOuterCorrectors         = 0
			self.nCorrectors              = 0

		if self.solver == 'interFoam' or self.solver == 'interDyMFoam':
			self.pressureRelaxation = 1
			self.velocityRelaxation = 1

			self.pressureTolerance = 1e-7
			self.velocityTolerance = 1e-8

			self.pressureRelTol      = 0.001
			self.velocityRelTol      = 0.001
			self.pressureFinalRelTol = 1e-4
			self.velocityFinalRelTol = 0.0
		elif self.solver == 'simpleFoam':
			self.pressureRelaxation = 0.9
			self.velocityRelaxation = 0.9

			self.pressureTolerance = 1e-8
			self.velocityTolerance = 1e-9

			self.pressureRelTol = 0.01
			self.velocityRelTol = 0.0
		else:
			self.pressureRelaxation = 0.3
			self.velocityRelaxation = 0.7

			self.pressureTolerance = 1e-6
			self.velocityTolerance = 1e-7

			self.pressureRelTol      = 0.001
			self.velocityRelTol      = 0.001
			self.pressureFinalRelTol = 0.0001
			self.velocityFinalRelTol = 0.0

		self.pressureResidual   = 1e-5
		self.velocityResidual   = 1e-6
		self.turbulenceResidual = 1e-5

	def setup(self):
		self.FoamFile = OrderedDict()
		self.FoamFile['version']  = 2.0
		self.FoamFile['format']   ='ascii'
		self.FoamFile['class']    ='dictionary'
		self.FoamFile['location'] ='"system"'
		self.FoamFile['object']   ='fvSolution'

		self.solverDict = OrderedDict()

		if self.solver == 'interFoam' or self.solver == 'interDyMFoam':
			self.solverDict['"alpha.water.*"'] = OrderedDict()
			self.solverDict['"alpha.water.*"']['nAlphaCorr']      = 3 
			self.solverDict['"alpha.water.*"']['nAlphaSubCycles'] = 1
			self.solverDict['"alpha.water.*"']['cAlpha']          = 1
			self.solverDict['"alpha.water.*"']['icAlpha']         = 0

			self.solverDict['"alpha.water.*"']['MULESCorr']          = 'yes'
			self.solverDict['"alpha.water.*"']['nLimiterIter']       = 10
			self.solverDict['"alpha.water.*"']['alphaApplyPrevCorr'] = 'yes'

			self.solverDict['"alpha.water.*"']['solver']             = 'smoothSolver'
			self.solverDict['"alpha.water.*"']['smoother']           = 'symGaussSeidel'
			self.solverDict['"alpha.water.*"']['tolerance']          = 1e-10
			self.solverDict['"alpha.water.*"']['relTol']             = 0
			self.solverDict['"alpha.water.*"']['minIter']            = 1

			self.solverDict['"pcorr.*"'] = OrderedDict()
			self.solverDict['"pcorr.*"']['solver']                = 'GAMG'
			self.solverDict['"pcorr.*"']['smoother']              = 'DIC'
			self.solverDict['"pcorr.*"']['tolerance']             = 0.1
			self.solverDict['"pcorr.*"']['relTol']                = 0
			if self.version < 4:
				self.solverDict['"pcorr.*"']['agglomerator']          = 'faceAreaPair'
				self.solverDict['"pcorr.*"']['mergeLevels']           = 1
				self.solverDict['"pcorr.*"']['nCellsInCoarsestLevel'] = 10
				self.solverDict['"pcorr.*"']['cacheAgglomeration']    = 'true'
			
			pName = 'p_rgh'
		else:
			pName = 'p'

		if self.solver == 'pimpleDyMFoam':
			self.solverDict['"pcorr.*"'] = OrderedDict()
			self.solverDict['"pcorr.*"']['solver']                = 'GAMG'
			self.solverDict['"pcorr.*"']['smoother']              = 'GaussSeidel'
			self.solverDict['"pcorr.*"']['tolerance']             = 0.1
			self.solverDict['"pcorr.*"']['relTol']                = 0
			if self.version < 4:
				self.solverDict['"pcorr.*"']['agglomerator']          = 'faceAreaPair'
				self.solverDict['"pcorr.*"']['mergeLevels']           = 1
				self.solverDict['"pcorr.*"']['nCellsInCoarsestLevel'] = 10
				self.solverDict['"pcorr.*"']['cacheAgglomeration']    = 'true'
			

		self.solverDict[pName] = OrderedDict()
		self.solverDict[pName]['solver']                = 'GAMG'
		self.solverDict[pName]['smoother']              = 'GaussSeidel'
		self.solverDict[pName]['tolerance']             = self.pressureTolerance
		self.solverDict[pName]['relTol']                = self.pressureRelTol
		
		if self.version < 4:
			self.solverDict[pName]['nPreSweeps']            = 0
			self.solverDict[pName]['nPostSweeps']           = 2
			self.solverDict[pName]['cacheAgglomeration']    = 'true'
			self.solverDict[pName]['agglomerator']          = 'faceAreaPair'
			self.solverDict[pName]['nCellsInCoarsestLevel'] = 10
			self.solverDict[pName]['mergeLevels']           = 1

		if self.solver != 'simpleFoam':
			self.solverDict[pName + 'Final'] = OrderedDict()
			self.solverDict[pName + 'Final']['']       = '$' + pName
			self.solverDict[pName + 'Final']['relTol'] = self.pressureFinalRelTol

		self.solverDict[self.velocitySolverKey] = OrderedDict()
		self.solverDict[self.velocitySolverKey]['solver']    = 'smoothSolver'
		self.solverDict[self.velocitySolverKey]['smoother']  = 'GaussSeidel'
		self.solverDict[self.velocitySolverKey]['tolerance'] = self.velocityTolerance
		self.solverDict[self.velocitySolverKey]['relTol']    = self.velocityRelTol
		if self.solver == 'interFoam' or self.solver == 'interDyMFoam':
			self.solverDict[self.velocitySolverKey]['minIter'] = 1
		if self.version < 4:
			self.solverDict[self.velocitySolverKey]['nSweeps']   = 2

		if self.solver != 'simpleFoam':
			self.solverDict[self.velocitySolverKey[0:-1] + 'Final"'] = OrderedDict()
			self.solverDict[self.velocitySolverKey[0:-1] + 'Final"'][''] = '$U'
			self.solverDict[self.velocitySolverKey[0:-1] + 'Final"']['relTol'] = self.velocityFinalRelTol



		self.solverDict['"cellDisplacement.*"'] = OrderedDict()
		self.solverDict['"cellDisplacement.*"']['solver']                = 'GAMG'
		self.solverDict['"cellDisplacement.*"']['smoother']              = 'GaussSeidel'
		self.solverDict['"cellDisplacement.*"']['tolerance']             =  1e-9
		self.solverDict['"cellDisplacement.*"']['relTol']                = 0
		if self.version < 4:
			self.solverDict['"cellDisplacement.*"']['agglomerator']          = 'faceAreaPair'
			self.solverDict['"cellDisplacement.*"']['mergeLevels']           = 1
			self.solverDict['"cellDisplacement.*"']['nCellsInCoarsestLevel'] = 10
			self.solverDict['"cellDisplacement.*"']['cacheAgglomeration']    = 'true'
		
		self.solverSettings = OrderedDict()

		if self.solverName == 'PIMPLE':
			if self.solver == 'interFoam' or self.solver == 'interDyMFoam':
				self.solverSettings['momentumPredictor'] = 'yes'
			else:
				self.solverSettings['momentumPredictor'] = 'yes'

			self.solverSettings['nOuterCorrectors']         = self.nOuterCorrectors
			self.solverSettings['nCorrectors']              = self.nCorrectors
			self.solverSettings['nNonOrthogonalCorrectors'] = self.nNonOrthogonalCorrectors

			if self.localTimeStep:
				self.solverSettings['maxCo']      = 10
				self.solverSettings['maxAlphaCo'] = 5

				self.solverSettings['rDeltaTSmoothingCoeff'] = 0.05
				self.solverSettings['rDeltaTDampingCoeff']   = 0.5
				self.solverSettings['nAlphaSpreadIter']      = int(0)
				self.solverSettings['nAlphaSweepIter']       = int(0)
				self.solverSettings['maxDeltaT']             = int(1)

			if self.solver == 'pimpleDyMFoam' or self.solver == 'interDyMFoam':
				self.solverSettings['correctPhi']      = '        yes'
				self.solverSettings['moveMeshOuterCorrectors'] = 'yes'
				self.solverSettings['turbOnFinalIterOnly'] = '    yes'
		elif self.solverName == 'SIMPLE':
			self.solverSettings['nNonOrthogonalCorrectors'] = self.nNonOrthogonalCorrectors
			self.solverSettings['consistent'] = 'yes'
			self.solverSettings['residualControl'] = OrderedDict()
			self.solverSettings['residualControl']['p'] = self.pressureResidual
			
			self.solverSettings['residualControl']['U'] = self.velocityResidual

		self.relaxationFactors = OrderedDict()
		self.relaxationFactors['fields']    = OrderedDict()
		self.relaxationFactors['fields']['p'] = self.pressureRelaxation
		if self.solverName == 'PIMPLE':
			self.relaxationFactors['fields']['pFinal'] = 1

		self.relaxationFactors['equations'] = OrderedDict()
		self.relaxationFactors['equations'][self.velocitySolverKey] = self.velocityRelaxation
		if self.solverName == 'PIMPLE':
			self.relaxationFactors['equations'][self.velocitySolverKey[0:-1] + 'Final"'] = 1

		self.cache = OrderedDict()
		self.cache[''] = 'grad(U)'


	def write(self, folder, ending=''):
		self.setup()

		f = open(folder + 'fvSolution' + ending, 'w')

		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		self.writeDict(f, self.solverDict, 'solvers')

		self.writeDict(f, self.solverSettings, self.solverName)
		self.writeDict(f, self.relaxationFactors, 'relaxationFactors')

		self.writeDict(f, self.cache, 'cache')

		f.close()