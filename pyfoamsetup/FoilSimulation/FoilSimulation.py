import numpy as np
import scipy.interpolate as interpolate
import os
import subprocess
import shutil
import multiprocessing
import sys
import copy

from pyfoamsetup.coreLibrary import *
import pyfoamsetup.coreLibrary.CaseSetup as CaseSetup

class FoilSimulation(CaseSetup.CaseSetup):
	def __init__(self, runName, U, c = 1, threeD = False, vilje = False, meshSetting='medium', fluid='air', timeStepSetting='medium'):
		if fluid == 'air':
			nu  = 1.45e-5
			rho = 1.226
		elif fluid == 'water':
			nu  = 1.19e-6
			rho = 999.1
		else:
			print('Unknown fluid! Assuming air...')
			nu  = 1.45e-5
			rho = 1.226

		A = c

		patchList = ['foil']
		self.nrFoils = 1

		self.wallFunction = True
		self.uniformWallRefinement = False

		# Call init from base class
		self.homePath = os.path.dirname(os.path.realpath(__file__))
		super().__init__(runName, patchList, c, U, A, nu, rho, 'FoilSimulation')
		self.splitPatches = True

		# Input variables
		self.threeD = threeD
		
		# Default mesh settings
		self.baseSize     = 1*self.L
		self.domainFront  = 5*self.L
		self.domainWake   = 10*self.L
		self.domainWidth  = 0.25*self.L
		
		# Default simulation settings
		self.setSolver('pimpleFoam')
		#self.timeStepScheme = 'CrankNicolson'
		self.endTime = np.round(50*self.L/self.U, decimals=6)
		self.writeInterval = self.endTime/10
	
		self.turbulenceModel = 'kOmegaSST' 
		self.turbulenceType  = 'RAS'

		self.schemeType = 'robustAndAccurate'

		self.writeWallPressure = False

		# List conating motion information
		self.motion     = [False]
		self.motionData = [False]

		self.createPatch = False

	def setDefaultCellLengths(self):
		super().setDefaultCellLengths()
		
		self.maxBaseSize     = 0.1*self.L   # Ensure fine enough background mesh
		self.maxSmallestSize = 0.001*self.L # Ensure that the geometry is captured
		self.viscousLength   = 0.01*self.L  # Resolution in the 'viscous wake'

	def setMeshSettings(self):
		super().setMeshSettings()
		if self.uniformWallRefinement:
			self.finalLayerFactor = 0.75
		else:
			self.finalLayerFactor = 1.0

	def addFoil(self, c = 1):
		if self.nrFoils == 1:
			self.patchList = ['foil1']

		self.nrFoils += 1
		self.patchList.append('foil{:.0f}'.format(self.nrFoils))
		self.motion.append(False)
		self.motionData.append(False)

	def addMotion(self, t, x, y, omega, origin, foilIndex = 0):
		if foilIndex > self.nrFoils:
			print('Too large foil index. Ignoring your request for motion!')

		if self.solver != 'pimpleDyMFoam':
			self.setSolver('pimpleDyMFoam')
			self.timeStepScheme   = 'Euler'
			self.adjustTimeStep   = True
			self.maxCo            = 10

			self.forceOrigo = np.array([origin[0], 0, 0])

			self.purgeWrite = 0
		
		self.motion[foilIndex] = True
		self.motionData[foilIndex] = motion(t, x, y, omega, origin)

		# Write more time steps
		self.writeInterval = np.round(1/48, decimals = 6)

	def writeBlockMesh(self):
		blockMesh = BlockMesh.Dict()

		# Calculate minimum values for domain size
		xBack  =  self.domainWake
		xFront = -self.domainFront

		yRight =  self.domainFront
		yLeft =  -self.domainFront

		if self.threeD:
			width     = np.round(self.domainWidth, decimals=6)
			z_nrCells = np.round(width/self.baseSize)
		else:
			width = np.round(self.baseSize, decimals=6)
			z_nrCells = 1

		# Calculate number of cells in each direction
		x_nrCells = np.ceil((xBack - xFront)/self.baseSize)
		y_nrCells = np.ceil((yRight - yLeft)/self.baseSize)

		# Readjust domain size to fit nr cells
		xLength = self.baseSize*x_nrCells
		yLength = self.baseSize*y_nrCells
		zLength = self.baseSize*z_nrCells

		wakeFraction  = (self.domainWake/(self.domainWake + self.domainFront))
		frontFraction = (self.domainFront/(self.domainWake + self.domainFront))
		xFront = np.round(-xLength*frontFraction, decimals=6)
		xBack  = np.round(xLength*wakeFraction, decimals=6)
	
		yRight = np.round(yLength/2, decimals=6)
		yLeft  = np.round(-yLength/2, decimals=6)

		# Add data to blockmesh and write
		if self.threeD:
			blockMesh.addVertex([xFront, yLeft,  -width/2])
			blockMesh.addVertex([xBack,  yLeft,  -width/2])
			blockMesh.addVertex([xBack,  yRight, -width/2])
			blockMesh.addVertex([xFront, yRight, -width/2])

			blockMesh.addVertex([xFront, yLeft,  width/2])
			blockMesh.addVertex([xBack,  yLeft,  width/2])
			blockMesh.addVertex([xBack,  yRight, width/2])
			blockMesh.addVertex([xFront, yRight, width/2])

			blockMesh.addBlock([x_nrCells, y_nrCells, z_nrCells])
		else:
			blockMesh.addBlock([x_nrCells, y_nrCells, 1])

			blockMesh.addVertex([xFront, yLeft,  -width])
			blockMesh.addVertex([xBack,  yLeft,  -width])
			blockMesh.addVertex([xBack,  yRight, -width])
			blockMesh.addVertex([xFront, yRight, -width])

			blockMesh.addVertex([xFront, yLeft,  0])
			blockMesh.addVertex([xBack,  yLeft,  0])
			blockMesh.addVertex([xBack,  yRight, 0])
			blockMesh.addVertex([xFront, yRight, 0])
		
		blockMesh.addBoundary('inlet',  'patch', [[0, 4, 7, 3], [3, 7, 6, 2], [0, 1, 5, 4]])
		blockMesh.addBoundary('outlet', 'patch', [[2, 6, 5, 1]])

		if self.threeD:
			blockMesh.addBoundary('leftMesh',  'patch', [[0, 3, 2, 1]])
			blockMesh.addBoundary('rightMesh', 'patch', [[4, 5, 6, 7]])
		else:
			blockMesh.addBoundary('left',  'empty', [[0, 3, 2, 1]])
			blockMesh.addBoundary('right', 'empty', [[4, 5, 6, 7]])

		blockMesh.write(self.systemFolder)

	def writeMesh(self):
		self.calculateBaseSize()
		self.writeBlockMesh()

		if self.uniformWallRefinement:
			minWallRefinement = self.maxRefinementLevel
		else:
			minWallRefinement = self.maxRefinementLevel-1

		# Initialize snappy hex mesh with foil
		if self.nrFoils == 1:
			self.snappyDict.addGeometry('foil.obj', 'triSurfaceMesh', {'name':'foil'})
			self.snappyDict.addRefinementSurface('foil', minWallRefinement, self.maxRefinementLevel, self.nrLayers)
			self.snappyDict.addFeature('foil.eMesh', self.maxRefinementLevel)
		else:
			for i in range(self.nrFoils):
				self.snappyDict.addGeometry('foil{:.0f}.obj'.format(i+1), 'triSurfaceMesh', {'name':'foil{:.0f}'.format(i+1)})
				self.snappyDict.addRefinementSurface('foil{:.0f}'.format(i+1), minWallRefinement, self.maxRefinementLevel, self.nrLayers)
				self.snappyDict.addFeature('foil{:.0f}.eMesh'.format(i+1), self.maxRefinementLevel)

		# Set up layer settings
		self.snappyDict.addLayersControls['relativeSizes']       = 'false'
		self.snappyDict.addLayersControls['finalLayerThickness'] = self.t_final
		self.snappyDict.addLayersControls['minThickness']        = 0.1*self.t_final
		self.snappyDict.addLayersControls['expansionRatio']      = self.layerExpansion

		self.snappyDict.castellatedMeshControls['nCellsBetweenLevels'] = int(self.nCellsBetweenLevels)

		# Set point in mesh
		if self.threeD:
			self.snappyDict.castellatedMeshControls['locationInMesh'] = '(-2.02 1.01 {:.6f})'.format(0.1*self.domainWidth)
		else:
			self.snappyDict.castellatedMeshControls['locationInMesh'] = '(-2.02 1.01 {:.6f})'.format(-0.1*self.baseSize)

			self.snappyDict.actionControl['snap'] = 'false'
			self.snappyDict.actionControl['addLayers'] = 'false'

			# Create second snappyHexMeshDict for snapping and layer stage
			snappyDict2 = copy.deepcopy(self.snappyDict)

			snappyDict2.actionControl['castellatedMesh'] = 'false'
			snappyDict2.actionControl['snap']            = 'true'
			snappyDict2.actionControl['addLayers']       = 'true'

			snappyDict2.meshQualityControls['minDeterminant'] = 1e-6

			snappyDict2.write(self.systemFolder, ending='.snap')

		self.snappyDict.write(self.systemFolder)

		if self.nrFoils == 1:
			self.snappyDict.writeSurfaceFeatureExtractDict(self.systemFolder, 'foil.obj')
		else:
			geoemtryFileList = []
			for i in range(self.nrFoils):
				geoemtryFileList.append('foil{:.0f}.obj'.format(i+1))
			self.snappyDict.writeSurfaceFeatureExtractDict(self.systemFolder, geoemtryFileList)

	def writeCaseFiles(self):
		# Update reference area
		if self.threeD:
			self.A = self.L*self.domainWidth

		super().writeCaseFiles()

		# Adjust boundary condition in case files
		if self.threeD:
			boundaryFiles = os.listdir(self.boundaryFolder)

			for i in range(len(boundaryFiles)):
				boundaryFile = boundaryFiles[i]
				filePath = self.boundaryFolder + boundaryFile

				if os.path.isfile(filePath):
					FileHandling.changeBoundary(filePath, 'left',  '\t\ttype    cyclicAMI;\n')
					FileHandling.changeBoundary(filePath, 'right', '\t\ttype    cyclicAMI;\n')

			createPatchDict = createPatch.Dict()
			createPatchDict.addTranslationalPatch('left', ('leftMesh'), 'right', [0, 0, self.domainWidth])
			createPatchDict.addTranslationalPatch('right', ('rightMesh'), 'left', [0, 0, -self.domainWidth])
			createPatchDict.write(self.systemFolder)

		# Right numerical schemes
		if self.turbulenceType == 'LES':
			os.remove(self.systemFolder+'fvSchemes')
			shutil.move(self.systemFolder+'fvSchemes.LES', self.systemFolder+'fvSchemes')
			os.remove(self.systemFolder+'fvSolution')
			shutil.move(self.systemFolder+'fvSolution.LES', self.systemFolder+'fvSolution')

		if self.nrFoils > 1:
			boundaryFiles = os.listdir(self.boundaryFolder)

			for i in range(len(boundaryFiles)):
				boundaryFile = boundaryFiles[i]
				filePath = self.boundaryFolder + boundaryFile

				if os.path.isfile(filePath):
					FileHandling.changeLine(filePath, 'foil',  '\tfoil1')

					for j in range(1, self.nrFoils):
						FileHandling.copyBoundary(filePath, 'foil1', 'foil{:.0f}'.format(j+1))

		if True in self.motion:
			f = open(self.systemFolder+'controlDict', 'a')
			f.write('\nlibs\n(\n\t"libmyfvMotionSolvers.so"\n);\n')
			f.close()

			dynamicMeshDict = dynamicMesh.displacement(self.patchList, farFieldPatch ='inlet')
			dynamicMeshDict.write(self.constantFolder)

			for i in range(self.nrFoils):
				if self.motion[i]:
					if self.nrFoils == 1:
						foilName = 'foil'
					else:
						foilName = 'foil{:.0f}'.format(i+1)

					self.motionData[i].write(self.runFolder+'motion_'+foilName)

					boundaryString = '\t\ttype\ttwoDimMotion;\n\t\torigin\t({:.6f} {:.6f} {:.6f});\n\t\ttimeDataFileName\t"$FOAM_CASE/motion_'.format(self.motionData[i].origin[0], self.motionData[i].origin[1], self.motionData[i].origin[2])+foilName+'";\n\t\tvalue\tuniform (0 0 0);\n'
					FileHandling.changeBoundary(self.boundaryFolder+'pointDisplacement', foilName, boundaryString)

					FileHandling.changeBoundary(self.boundaryFolder+'U', foilName, '\t\ttype    movingWallVelocity;\n\t\tvalue    uniform (0 0 0);\n')

		if self.createPatch:
			self.createPatchDict.write(self.systemFolder)

	def writeScripts(self, boundaryLayer=True):
		# Mesh
		f = open(self.runFolder+'/mesh.sh', 'w')

		f.write('#!/bin/bash\n\n')

		f.write('surfaceFeatureExtract\n')
		f.write('blockMesh\n')

		# Execute castellation in parallel
		f.write('mv system/decomposeParDict system/decomposeParDict.sim\n')
		f.write('mv system/decomposeParDict.mesh system/decomposeParDict\n')

		f.write('decomposePar\n')
			
		f.write('mpirun -np {:.0f} snappyHexMesh -overwrite -parallel\n'.format(self.nCPUs_mesh))

		f.write('reconstructParMesh -constant\n')
		f.write('rm -fr processor*\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.mesh\n')
		f.write('mv system/decomposeParDict.sim system/decomposeParDict\n')

		if self.threeD:
			f.write('createPatch -overwrite\n')
		else:
			f.write('extrudeMesh\n')

			f.write('mv system/snappyHexMeshDict system/snappyHexMeshDict.castellation\n')
			f.write('mv system/snappyHexMeshDict.snap system/snappyHexMeshDict\n')

			f.write('snappyHexMesh -overwrite\n')

		if len(self.topoSetList) > 0:
			f.write('topoSet\n')

		if self.createPatch:
			f.write('createPatch -overwrite\n')

		f.write('renumberMesh -overwrite\n')	
			
		f.close()

		# Simulation
		f = open(self.runFolder + '/runSim.sh', 'w')

		f.write('#!/bin/bash\n\n')

		if self.solver == 'simpleFoam':
			self.calculateBoundaryLayerHeight()
			f.write('applyBoundaryLayer -ybl {:.6f}\n\n'.format(self.boundaryLayerHeight))

		f.write('decomposePar\n')

		if self.vilje:
			f.write('mpiexec ' + self.solver + ' -parallel\n')
		else:
			if True in self.motion:
				f.write('mpirun -np {:.0f} '.format(self.nCPUs) + 'pimpleDyMFoam' + ' -parallel\n')
			else:
				f.write('mpirun -np {:.0f} '.format(self.nCPUs) + self.solver + ' -parallel\n')

		f.write('reconstructPar\n')

		f.write('rm -fr processor*\n')

		f.close()

	def addViscousWake(self, x0, y0, lengthFactor = 2, widthFactor = 0.5, expansion=1.5):
		nrLayers = copy.deepcopy(self.nrLayers)
		self.calculateBaseSize()
		self.nrLayers = nrLayers

		length0 = self.L*lengthFactor
		width0  = self.L*widthFactor

		maxLevel = self.computeCellLevel(self.viscousLength)

		if isinstance(x0, list):
			nrWakes = len(x0)
		else:
			nrWakes = 1

		for j in range(nrWakes):
			level = maxLevel
			for i in range(maxLevel):
				name = 'viscBox_wake{:.0f}_{:.0f}'.format(j+1, i+1)
				length = length0*expansion**(i)
				cellLength = self.baseSize/(2**level+1)

				if i == 0:
					width  = width0
				else:
					width  += 2*cellLength*self.nCellsBetweenLevels

				width  = width0*expansion**(i)

				distance = width/2

				if isinstance(x0, list):
					minString = '({:.6f} {:.6f} -999)'.format(x0[j],          y0[j] - 1e-6)
					maxString = '({:.6f} {:.6f} 999)'.format( x0[j] + length, y0[j] + 1e-6)
				else:
					minString = '({:.6f} {:.6f} -999)'.format(x0,          y0 - 1e-6)
					maxString = '({:.6f} {:.6f} 999)'.format( x0 + length, y0 + 1e-6)

				self.snappyDict.addGeometry(name, 'searchableBox', {'min':minString, 'max':maxString})
				self.snappyDict.addRefinementRegion(name, 'distance', np.array([distance, level]))

				level -= 1

		distance = width0/2
		if self.nrFoils == 1:
			self.snappyDict.addRefinementRegion('foil', 'distance', np.array([distance, maxLevel]))
		else:
			for i in range(self.nrFoils):
				self.snappyDict.addRefinementRegion('foil{:.0f}'.format(i+1), 'distance', np.array([distance, maxLevel]))

	def addRotation(self, omega):
		boundaryString = '\t\ttype rotatingWallVelocity;\n\t\torigin (0.5 0 0);\n\t\taxis (0 0 -1);\n\t\tomega {};\n'.format(omega)

		FileHandling.changeBoundary(self.boundaryFolder+'U', 'foil', boundaryString)


def writeRunScript(caseNameList, foilSim):
	# Write run script
	f = open('run.sh', 'w')
	f.write('#!/bin/bash\n\n')

	if sys.platform == 'darwin':
		f = open(foilSim.foamPath + '/run.sh', 'w')
	else:
		f = open('run.sh', 'w')

	f.write('#!/bin/bash\n\n')

	if sys.platform == 'darwin':
		f.write('cd $HOME/workingDir/OpenFOAM/run/FoilSimulation\n\n')
	else:
		f.write('cd $FOAM_RUN/FoilSimulation\n\n')

	for i in range(len(caseNameList)):
		f.write('cd {0}\n'.format(caseNameList[i]))
		f.write('sh mesh.sh\n')

		f.write('sh runSim.sh\n')

		if sys.platform == 'darwin':
			f.write('cd $HOME/workingDir/OpenFOAM/run/FoilSimulation\n\n')
		else:
			f.write('cd $FOAM_RUN/FoilSimulation\n\n')

	f.close()	