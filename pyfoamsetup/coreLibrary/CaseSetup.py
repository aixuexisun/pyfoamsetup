import numpy as np
import scipy.interpolate as interpolate
import os
import subprocess
import shutil
from collections import OrderedDict

from pyfoamsetup.coreLibrary import *

import multiprocessing
import sys

class motion():
	def __init__(self, t, x, y, omega, origin):
		self.t      = t
		self.x      = x
		self.y      = y
		self.omega  = omega
		self.origin = origin

		self.n = len(t)

		self.x_spl     = interpolate.splrep(self.t, self.x)
		self.y_spl     = interpolate.splrep(self.t, self.y)
		self.omega_spl = interpolate.splrep(self.t, self.omega)

	def write(self, filePath):
		f = open(filePath, 'w')

		f.write('(\n')
		for i in range(self.n):
			f.write('\t({:.12f} ({:.12f} {:.12f} {:.12f}))\n'.format(self.t[i], self.x[i], self.y[i], self.omega[i]))

		f.write(');')

		f.close()

class CaseSetup:
	def __init__(self, runName, patchList, L, U, A, nu, rho, applicationName):
		self.runName = runName
		self.L       = L
		self.U       = U
		self.A       = A
		self.nu      = nu
		self.rho     = rho
		self.Re      = self.L*self.U/self.nu

		self.applicationName = applicationName

		self.superComputer = False # Switch that turns on settings for simulations on a cluster

		self.wallFunction = True

		self.patchList = patchList
		
		self.snappyDict = SnappyHexMesh.Dict()

		self.forceOrigo = np.zeros(3)

		self.turbulenceModel = 'kOmegaSST' # Available models: kOmega, kOmegaSST, kEpsilon, realizableKE, SpalartAllmaras, SpalartAllmarasDDES, 
		self.turbulenceType  = 'RAS'

		self.schemeType     = 'robustAndAccurate'
		self.timeStepScheme = 'Euler'
		self.localTimeStep  = False
		self.adjustTimeStep = True

		self.nOuterCorrectors         = 3
		self.nCorrectors              = 2
		self.nNonOrthogonalCorrectors = 0

		self.fvOptionsList   = []
		self.topoSetList     = []
		self.fvOptionsNames  = []

		self.decompositionMethod = 'scotch' # Options: scotch, simple
		self.setComputerSettings()
		self.setMeshSettings()
		self.setDefaultCellLengths()
		self.setDefaultTimeStep()
		self.calculateBoundaryLayerHeight()
		self.setFolderPaths()

		self.endTime = 20*self.L/self.U

		# Force dict options
		self.splitPatches = False

		# OpenFoam version
		if sys.platform == 'darwin':
			self.version = 4
		else:
			self.version = 4

		self.samples = Sampling.Dict()

		self.purgeWrite = 2

	def setFolderPaths(self):
		self.foamPath       = os.environ['FOAM_RUN']
		self.mainRunFolder  = self.foamPath + '/' + self.applicationName
		self.baseFolder     = self.homePath + '/BaseFolder'

		self.runFolder      = self.mainRunFolder + '/' + self.runName + '/'
		self.systemFolder   = self.runFolder + 'system/'
		self.boundaryFolder = self.runFolder + '0/'
		self.constantFolder = self.runFolder + 'constant/'
		self.geometryFolder = self.runFolder + 'constant/triSurface/'

		# Check if main path exists
		if not os.path.exists(self.mainRunFolder):
			os.makedirs(self.mainRunFolder)

	def createMissingFolders(self):
		pathList = ['constant', 'constant/triSurface', 'system']

		for path in pathList:
			if not(os.path.exists(self.runFolder+path)):
				os.mkdir(self.runFolder+path)

	def setComputerSettings(self, nCPU=None):
		if sys.platform == 'darwin':
			cpu_count = multiprocessing.cpu_count()
		else:
			cpu_count = multiprocessing.cpu_count()/2

		if cpu_count > 16:
			self.nCPUs_mesh = 16
		else:
			self.nCPUs_mesh = cpu_count

		if not(nCPU==None):
			self.nCPUs = nCPU
		else:
			self.nCPUs = cpu_count
	
	def setMeshSettings(self):
		self.nCellsBetweenLevels = 5
		self.maxRefinementLevel  = 10
		self.minRefinementLevel  = 2
		self.finalLayerFactor    = 0.5

		if self.wallFunction:
			self.layerExpansion    = 1.3
			self.minLayerExpansion = 1.1
			self.nrLayers          = 3
			self.maxNrLayers       = 8
			self.yPlusMin          = 30
		else:
			self.nrLayers          = 5
			self.maxNrLayers       = 20
			self.yPlusMin          = 0
			self.layerExpansion    = 1.25
			self.minLayerExpansion = 1.1

	def setDefaultCellLengths(self):
		self.maxSmallestSize = 0.001*self.L
		self.maxBaseSize     = 1.0*self.L
		self.viscousLength   = 0.01*self.L

		if self.wallFunction:
			self.yPlus = 60
		else:
			self.yPlus = 1.0

	def setDefaultTimeStep(self):
		self.deltaT     = 0.01*self.L/self.U
		self.maxDeltaT  = 0.01*self.L/self.U
		self.maxCo      = 10
		self.maxAlphaCo = 5

	def changeTimeStep(self, timeStepFactor):
		self.setDefaultTimeStep()

		self.deltaT     *= timeStepFactor
		self.maxDeltaT  *= timeStepFactor
		self.maxCo      *= timeStepFactor
		self.maxAlphaCo *= timeStepFactor

	def changeCellLengths(self, cellLengthFactor):
		self.setDefaultCellLengths()

		self.maxBaseSize     *= cellLengthFactor
		self.maxSmallestSize *= cellLengthFactor
		self.viscousLength   *= cellLengthFactor

		self.yPlus *= cellLengthFactor

	def setSolver(self, solver, endTime=None, localTimeStep=False):
		self.solver = solver

		if self.solver == 'interFoam' or self.solver == 'interDyMFoam':
			self.localTimeStep = localTimeStep

			self.nOuterCorrectors         = 1
			self.nCorrectors              = 2
			self.nNonOrthogonalCorrectors = 0

			if localTimeStep:
				self.deltaT  = 1
		elif self.solver == 'simpleFoam':
			self.deltaT = 1
		elif self.solver == 'pimpleFoam' or self.solver == 'pimpleDyMFoam':
			self.nOuterCorrectors         = 3
			self.nCorrectors              = 2
			self.nNonOrthogonalCorrectors = 0
		else:
			print('Unknown solver!')

		if endTime != None:
			self.endTime       = endTime
			self.writeInterval = endTime/10

		self.fvSolutionDict = fvSolution.Dict(self.solver, self.turbulenceModel, version=self.version, localTimeStep=localTimeStep)

	def setLocalTimeStep(self):
		self.localTimeStep = True
		self.setSolver('interFoam', localTimeStep=True)
		#self.fvSolutionDict.solverSettings['maxCo']      = self.maxCo
		#self.fvSolutionDict.solverSettings['maxAlphaCo'] = self.maxAlphaCo
		self.timeStepScheme = 'localEuler'
		self.deltaT         = 1
		self.maxDeltaT      = 1
		self.adjustTimeStep = False

	def setTurbulenceModel(self, turbulenceModel, turbulenceType='RAS'):
		self.turbulenceModel = turbulenceModel
		self.turbulenceType  = turbulenceType

		self.fvSolutionDict = fvSolution.Dict(self.solver, self.turbulenceModel)

	def Cf(self, Re):
		# Numerical friction line based on k-omega SST
		a1 =  0.089
		a2 = -0.283
		a3 =  4.73e-3
		a4 = -2.43e-5

		Cf = a1*Re**(a2+a3*np.log(Re) + a4*np.log(Re)**2) 

		return Cf

	def calculateBoundaryLayerHeight(self):
		self.boundaryLayerHeight = 0.37*self.L/self.Re**(1/5)

	def calculateFirstLayerThickness(self):
		Cf = self.Cf(self.Re)
		U  = self.Re*self.nu/self.L # Calculate representative velocity from the Reynolds number
		
		Tau = Cf*0.5*U**2

		u_t = np.sqrt(Tau)

		self.t_first = 2*self.yPlus*self.nu/u_t

	def calculateyPlus(self):
		Cf = self.Cf(self.Re)
		U = self.Re*self.nu/self.L

		Tau = Cf*0.5*U**2

		u_t = np.sqrt(Tau)

		self.yPlus = self.t_first*u_t/(2*self.nu)

	def adjustNumberOfLayers(self):
		self.calculateBoundaryLayerHeight()

		delta = self.boundaryLayerHeight
		a     = self.t_first
		r     = self.layerExpansion

		self.nrLayers = int(np.round(np.log(1 - (delta/a)*(1-r))/np.log(r)))

	def calculateBaseSize(self):
		self.calculateFirstLayerThickness()

		if self.nrLayers > 0:
			self.t_final      = self.t_first*self.layerExpansion**(self.nrLayers - 1)
			self.smallestSize = self.t_final/self.finalLayerFactor

			while self.smallestSize < self.maxSmallestSize and self.nrLayers < self.maxNrLayers:
				self.nrLayers     += 1
				self.t_final      = self.t_first*self.layerExpansion**(self.nrLayers - 1)
				self.smallestSize = self.t_final/self.finalLayerFactor

			if self.smallestSize > self.maxSmallestSize:
				'''Try to reduce layer expansion, or if that dont work, reduce yPlus value'''
				necessaryExpansion = (self.maxSmallestSize*self.finalLayerFactor/self.t_first)**(1/(self.nrLayers-1))

				if necessaryExpansion > self.minLayerExpansion:
					self.layerExpansion = necessaryExpansion
					self.t_final        = self.t_first*self.layerExpansion**(self.nrLayers - 1)
					self.smallestSize   = self.t_final/self.finalLayerFactor
				else:	
					self.layerExpansion = self.minLayerExpansion
					self.smallestSize = self.maxSmallestSize
					self.t_final      = self.smallestSize*self.finalLayerFactor
					self.t_first      = self.t_final/(self.layerExpansion**(self.nrLayers - 1))
					self.calculateyPlus()

					if self.wallFunction and self.yPlus < self.yPlusMin:
						print('Warning! y+ below {:.3f}!'.format(self.yPlusMin))

		else:
			self.smallestSize = self.t_first
			self.t_final = self.t_first
			self.layerExpansion = 1

			if self.smallestSize > self.maxSmallestSize:
				self.smallestSize = self.maxSmallestSize

				self.calculateyPlus()

				if self.wallFunction and self.yPlus < self.yPlusMin:
					print('Warning! y+ below {:.3f}!'.format(self.yPlusMin))
	
		self.baseSize = self.smallestSize*2**(self.maxRefinementLevel)

		while self.baseSize > self.maxBaseSize and self.maxRefinementLevel > self.minRefinementLevel:
			self.maxRefinementLevel -= 1
			self.baseSize = self.smallestSize*2**(self.maxRefinementLevel)

	def writeCaseFiles(self):
		# Delete old folder with runName inside towing tank folder
		if os.path.exists(self.runFolder):
			shutil.rmtree(self.runFolder)
			
		# Copy over base folder to run folder
		shutil.copytree(self.baseFolder, self.runFolder)
		self.createMissingFolders()

		# Copy over all additional system files
		FileHandling.adjustSystemFolder(self.runFolder, self.solver)

		# Write mesh files
		self.writeMesh()

		# Set right controlDict settings
		self.setControlDict()
		self.controlDict.write(self.systemFolder)
		# Set up numerical schemes
		fvSchemesDict = fvSchemes.Dict(self.solver, self.turbulenceModel, schemeType=self.schemeType, timeStepScheme=self.timeStepScheme)
		fvSchemesDict.write(self.systemFolder)

		# Write fvSolver file
		if self.fvSolutionDict.solverName == 'PIMPLE':
			self.fvSolutionDict.nOuterCorrectors         = self.nOuterCorrectors
			self.fvSolutionDict.nCorrectors              = self.nCorrectors
			self.fvSolutionDict.nNonOrthogonalCorrectors = self.nNonOrthogonalCorrectors
		self.fvSolutionDict.write(self.systemFolder)

		# Set right turbulence model
		turbulencePropertiesDict = turbulenceProperties.Dict(self.U, self.L, self.nu, self.turbulenceModel, self.turbulenceType)
		turbulencePropertiesDict.calculateInitialValues()
		turbulencePropertiesDict.writeInitialFile(self.runFolder)
		turbulencePropertiesDict.write(self.constantFolder)
		# Adjust turbulence boundaries
		if not(self.wallFunction):
			for i in range(len(self.patchList)):
				FileHandling.changeBoundary(self.boundaryFolder+'epsilon', self.patchList[i], '\t\ttype    epsilonLowReWallFunction;\n\t\tvalue    $internalField;\n')
				FileHandling.changeBoundary(self.boundaryFolder+'k', self.patchList[i], '\t\ttype    kLowReWallFunction;\n\t\tvalue    $internalField;\n')
				FileHandling.changeBoundary(self.boundaryFolder +'nut', self.patchList[i], '\t\ttype     nutLowReWallFunction;\n\t\tvalue    $internalField;\n')

		# Write transport properties
		if self.solver == 'interFoam' or self.solver == 'interDyMFoam':
			transportPropertiesDict = transportProperties.Dict([self.nu, self.nu_air], rho=[self.rho, self.rho_air], freeSurface=True)
		else:
			transportPropertiesDict = transportProperties.Dict(self.nu, rho=self.rho, freeSurface=False)

		transportPropertiesDict.write(self.constantFolder)

		# Write decomposeParDict
		decomposeParDict = DecomposePar.Dict(self.nCPUs, self.decompositionMethod)
		decomposeParDict.write(self.systemFolder)

		decomposeParDictMesh = DecomposePar.Dict(self.nCPUs_mesh, 'simple')
		decomposeParDictMesh.write(self.systemFolder, ending='.mesh')

		# Write forceDict
		if self.solver == 'interFoam' or self.solver == 'interDyMFoam':
			interFoam = True
		else:
			interFoam = False

		if len(self.patchList) > 0:
			forces = Forces.Dict(self.patchList, self.rho, self.U, self.L, self.A, splitPatches=self.splitPatches, interFoam=interFoam, version=self.version)

			forces.CofR = self.forceOrigo
			forces.write(self.systemFolder)

		self.samples.write(self.systemFolder)

		# Write fvOptions file
		if len(self.fvOptionsList) > 0:
			f = open(self.constantFolder+'fvOptions', 'w')
			f.write(fvOptions.header())

			for option in self.fvOptionsList:
				f.write(option.generateString())

			f.close()

			f = open(self.systemFolder+'controlDict', 'a')
			f.write('\nlibs\n(\n\t"libfvOptions.so"\n\t"libmyfvOptions.so"\n);\n')
			f.close()

		# Write topoSet file
		if len(self.topoSetList) > 0:
			f = open(self.systemFolder+'topoSetDict', 'w')
			f.write(topoSet.header())

			f.write('actions\n(\n')

			for action in self.topoSetList:
				f.write(action.generateString())

			f.write(');\n')
			f.close()

			#self.createPatchDict.write(self.systemFolder)

		# Write sim info
		self.writeSimInfo()

	def setControlDict(self):
		self.controlDict = ControlDict.Dict(self.solver, localTimeStep=self.localTimeStep)

		if not(self.adjustTimeStep):
			self.deltaT    = min(self.deltaT, self.maxCo*self.smallestSize/self.U)
			self.maxDeltaT = min(self.deltaT, self.maxCo**self.smallestSize/self.U)

		if self.solver == 'simpleFoam' or self.localTimeStep:
			self.deltaT    = int(1)
			self.maxDeltaT = int(1)

			self.endTime       = int(self.endTime)
			self.writeInterval = int(self.writeInterval)
		else:
			self.deltaT        = roundSignificant(self.deltaT,        significantDigits=4)
			self.maxDeltaT     = roundSignificant(self.maxDeltaT,     significantDigits=4)
			self.endTime       = roundSignificant(self.endTime,       significantDigits=4)
			self.writeInterval = roundSignificant(self.writeInterval, significantDigits=4)
			self.maxCo         = roundSignificant(self.maxCo,         significantDigits=4)
			self.maxAlphaCo    = roundSignificant(self.maxAlphaCo,    significantDigits=4)
		
		self.controlDict.Settings['deltaT']        = self.deltaT
		self.controlDict.Settings['endTime']       = self.endTime
		self.controlDict.Settings['writeInterval'] = self.writeInterval

		self.controlDict.Settings['purgeWrite'] = self.purgeWrite

		if self.solver != 'simpleFoam' and not(self.localTimeStep):
			self.controlDict.Settings['maxCo']     = self.maxCo
			self.controlDict.Settings['maxDeltaT'] = self.maxDeltaT

			if self.adjustTimeStep:
				self.controlDict.Settings['adjustTimeStep'] = 'yes'
			else:
				self.controlDict.Settings['adjustTimeStep'] = 'no'
				if not(self.localTimeStep):
					self.controlDict.Settings['writeControl']   = 'runTime'

		if (self.solver == 'interFoam' or self.solver == 'interDyMFoam') and not(self.localTimeStep):
			self.controlDict.Settings['maxAlphaCo'] = self.maxAlphaCo

		if len(self.patchList) > 0:
			self.controlDict.FunctionList = ['\t#include "forces";']

		if self.samples.nSamples > 0:
			self.controlDict.FunctionList.append('\t#include "' + self.samples.fileName + '";')

		if self.superComputer:
			self.controlDict.Settings['runTimeModifiable'] = 'no'

	def computeCellLevel(self, Length):
		return int(np.ceil(np.log(self.baseSize/Length)/np.log(2)))

	def writeSimInfo(self):
		f = open(self.runFolder+'simulationInfo.txt', 'w')
		f.write('Simulation information\n')

		f.write('Re {:.2f}\n'.format(self.Re))
		f.write('L  {:.6f}\n'.format(self.L))
		f.write('yPlus {:.2f}\n'.format(self.yPlus))
		f.write('nrLayers {:.0f}\n'.format(self.nrLayers))
		f.write('layerExpansion {:.6f}\n'.format(self.layerExpansion))
		f.write('smallestSize/L {:.6f}\n'.format(self.smallestSize/self.L))
		f.write('baseSize/L {:.6f}\n'.format(self.baseSize/self.L))
		f.write('boundaryLayerHeight/L {:.6f}\n'.format(self.boundaryLayerHeight/self.L))
		f.write('refinementLevels {:.0f}\n'.format(self.maxRefinementLevel))
		f.write('nCellsBetweenLevels {:.0f}\n'.format(self.nCellsBetweenLevels))
		f.write('solver {}\n'.format(self.solver))

		f.close()

	def addPropellerActuatorDisk(self, name, p1, r, diskDir, diskThickness, Thrust, Torque, cellSize, rh_factor=0.1):
		# Add cylinders
		cellSetName = name+'_cellSet'
		level = self.computeCellLevel(cellSize)

		rh = rh_factor*0.1

		diskThickness = max(2*cellSize, diskThickness)

		p2 = p1 + diskThickness*diskDir
		p2_snappy = p1 + r*diskDir

		point1String = '({:.6f} {:.6f} {:.6f})'.format(p1[0], p1[1], p1[2])
		point2String = '({:.6f} {:.6f} {:.6f})'.format(p2_snappy[0], p2_snappy[1], p2_snappy[2])
		radiusString = '{:.6f}'.format(r)

		self.snappyDict.addGeometry(name, 'searchableCylinder', {'point1':point1String, 'point2':point2String, 'radius':radiusString})
		self.snappyDict.addRefinementRegion(name, 'inside', np.array([1, level]))

		actuatorDisk = fvOptions.goldsteinActuationDisk(name, Thrust, Torque, p1, p2, r, rh, cellSetName)		
		cellSet      = topoSet.cylinderToCell(cellSetName, p1, p2, r)

		self.fvOptionsList.append(actuatorDisk)
		self.topoSetList.append(cellSet)

def writeRunScripts(caseNameList, sim, folderName=''):
	# Write run script
	if sim.superComputer:
		filePath = sim.mainRunFolder + folderName + 'viljeRun.sh'
		foalderNameLength = len(folderName)
		f = open(filePath, 'w')

		f.write('#!/bin/bash\n')
		f.write('#PBS -A nn4040k\n')
		f.write('#PBS -N driftAngles\n')
		f.write('#PBS -l walltime={:.0f}:00:00\n'.format(24*len(caseNameList)))
		f.write('#PBS -l select={:.0f}:ncpus=32:mpiprocs=16\n\n'.format(sim.nCPUs/16))
		f.write('module load gcc/4.9.1\n')
		f.write('module load mpt/2.13\n')
		f.write('module load openfoam/3.0+\n\n')
		f.write('cd $PBS_O_WORKDIR\n\n')

		for i in range(len(caseNameList)):
			f.write('cd {0}\n'.format(caseNameList[i][foalderNameLength:]))
			f.write('bash runSim.sh\n')
			f.write('cd $PBS_O_WORKDIR\n\n')
	else:
		f = open('run.sh', 'w')
		f.write('#!/bin/bash\n\n')
		f.write('cd {0}\n\n'.format(sim.mainRunFolder))

		for i in range(len(caseNameList)):
			f.write('cd {0}\n'.format(caseNameList[i]))
			f.write('bash mesh.sh\n')
			f.write('bash runSim.sh\n')
			f.write('cd {0}\n\n'.format(sim.mainRunFolder))

	f.close()

	if sim.superComputer:
		f = open('localMeshing.sh', 'w')
		f.write('#!/bin/bash\n\n')
		f.write('cd {0}\n\n'.format(sim.mainRunFolder))

		for i in range(len(caseNameList)):
			f.write('cd {0}\n'.format(caseNameList[i]))
			f.write('bash mesh.sh\n')
			f.write('cd {0}\n\n'.format(sim.mainRunFolder))

		f.close()

def roundSignificant(variable, significantDigits=4):
	decimals = -int(np.round(np.log10(variable))) + significantDigits

	return np.round(variable, decimals = decimals)