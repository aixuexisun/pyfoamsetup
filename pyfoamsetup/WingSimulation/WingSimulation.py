import numpy as np
import os
import subprocess
import shutil
import multiprocessing
import sys
import copy

from pyfoamsetup.coreLibrary import *
import pyfoamsetup.coreLibrary.CaseSetup as CaseSetup

class WingSimulation(CaseSetup.CaseSetup):
	def __init__(self, runName, U, L, S, A, vilje = False, meshSetting='medium', fluid='air', timeStepSetting='medium'):
		# Default environment settings
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

		patchList = ['wing']
		self.nrWings = 1

		self.uniformWallRefinement = False
		self.meshSetting           = meshSetting
		self.timeStepSetting       = timeStepSetting 

		# Call init from base class
		self.homePath = os.path.dirname(os.path.realpath(__file__))
		super().__init__(runName, patchList, L, U, A, nu, rho, 'WingSimulation')

		# Extra length dimension
		self.S = S  # Span of wing

		self.baseSize = 1

		self.domainWake_span   = 4
		self.domainFront_span  = 2
		self.domainWidth_span  = 2
		self.domainHeight_span = 2

		self.domainWake_chord   = 10
		self.domainFront_chord  = 5
		self.domainWidth_chord  = 5
		self.domainHeight_chord = 5

		self.noSlipGround = False
		self.viscousWake  = False
		self.WIGSim       = False
		self.heightAboveGround = 1*self.L
		
		self.setMeshSettings()
		
		self.setSolver('pimpleFoam')


		self.endTime       = np.round(15*self.L/self.U, decimals=6)
		self.writeInterval = self.endTime/10

		# List conating motion information
		self.motion     = [False]
		self.motionData = [False]

	def setDefaultCellLengths(self):
		super().setDefaultCellLengths()
		
		self.maxBaseSize     = 0.2*self.L   # Ensure fine enough background mesh
		self.maxSmallestSize = 0.002*self.L # Ensure that the geometry is captured
		self.viscousLength   = 0.01*self.L  # Resolution in the 'viscous wake'

	def setMeshSettings(self):
		super().setMeshSettings()
		self.nCellsBetweenLevels = 3

		if self.uniformWallRefinement:
			self.finalLayerFactor = 0.75
		else:
			self.finalLayerFactor = 1.0

	def writeBlockMesh(self):
		blockMesh = BlockMesh.Dict()

		# Calculate minimum values for domain size
		xBack  =  max(self.domainWake_span*self.S, self.domainWake_chord*self.L)
		xFront = -max(self.domainFront_span*self.S, self.domainFront_chord*self.L)
		yRight =  max(self.domainWidth_span*self.S, self.domainWidth_chord*self.L)

		if self.WIGSim:
			yLeft = -self.heightAboveGround
		else:
			yLeft  = -yRight

		zHeight = max(self.domainHeight_span*self.S, self.domainHeight_chord*self.L)
		zDepth  = 0

		# Calculate number of cells in each direction
		x_nrCells = np.ceil((xBack - xFront)/self.baseSize)
		y_nrCells = np.ceil((yRight - yLeft)/self.baseSize)
		z_nrCells = np.ceil((zHeight - zDepth)/self.baseSize)

		# Readjust domain size to fit nr cells
		xLength = self.baseSize*x_nrCells
		zLength = self.baseSize*z_nrCells
		yLength = self.baseSize*y_nrCells
		

		wakeFraction  = abs(xBack)/(abs(xBack) + abs(xFront))
		frontFraction = abs(xFront)/(abs(xBack) + abs(xFront))
		xFront = -xLength*frontFraction
		xBack  =  xLength*wakeFraction
	
		if self.WIGSim:
			yLeft  = -self.heightAboveGround
			yRight = yLeft + yLength
		else:
			yLeft = -yLength/2
			yRight =  yLength/2

		# Add data to blockmesh and write
		blockMesh.addVertex([xFront, yLeft,  zDepth])
		blockMesh.addVertex([xBack,  yLeft,  zDepth])
		blockMesh.addVertex([xBack,  yRight, zDepth])
		blockMesh.addVertex([xFront, yRight, zDepth])

		blockMesh.addVertex([xFront, yLeft,  zHeight])
		blockMesh.addVertex([xBack,  yLeft,  zHeight])
		blockMesh.addVertex([xBack,  yRight, zHeight])
		blockMesh.addVertex([xFront, yRight, zHeight])

		blockMesh.addBlock([x_nrCells, y_nrCells, z_nrCells])
		
		if self.WIGSim:
			blockMesh.addBoundary('inlet',  'patch', [[0, 4, 7, 3]])
			blockMesh.addBoundary('sides',  'patch', [[3, 2, 6, 7]])
			blockMesh.addBoundary('ground', 'patch', [[0, 1, 5, 4]])
		else:
			blockMesh.addBoundary('inlet', 'patch', [[0, 4, 7, 3]])
			blockMesh.addBoundary('sides',  'patch', [[3, 2, 6, 7], [0, 1, 5, 4]])

		if not(self.noSlipGround):
			blockMesh.addBoundary('bottom', 'patch', [[0, 3, 2, 1]])
		else:
			blockMesh.addBoundary('bottom', 'wall', [[0, 3, 2, 1]])

		blockMesh.addBoundary('top', 'patch', [[4, 5, 6, 7]])

		blockMesh.addBoundary('outlet', 'patch', [[2, 6, 5, 1]])

		blockMesh.write(self.systemFolder)

	def writeMesh(self):
		self.calculateBaseSize()
		self.writeBlockMesh()

		self.snappyDict.addGeometry('wing.obj', 'triSurfaceMesh', {'name':'wing'})
		if self.uniformWallRefinement:
			self.snappyDict.addRefinementSurface('wing', self.maxRefinementLevel, self.maxRefinementLevel, self.nrLayers)
		else:
			self.snappyDict.addRefinementSurface('wing', self.maxRefinementLevel-1, self.maxRefinementLevel, self.nrLayers)
		self.snappyDict.addFeature('wing.eMesh', self.maxRefinementLevel)

		# Set up layer settings
		self.snappyDict.addLayersControls['relativeSizes']       = 'false'
		self.snappyDict.addLayersControls['finalLayerThickness'] = self.t_final
		self.snappyDict.addLayersControls['minThickness']        = 0.5*self.t_final
		self.snappyDict.addLayersControls['expansionRatio']      = self.layerExpansion

		self.snappyDict.castellatedMeshControls['locationInMesh'] = '({:.3f} {:.3f} {:.3f})'.format(-0.5*self.L, 1.5*self.L, 1.5*self.S)

		self.snappyDict.castellatedMeshControls['nCellsBetweenLevels'] = int(self.nCellsBetweenLevels)

		if self.noSlipGround:
			level = self.maxRefinementLevel-4
			cellLength = self.baseSize/(2**level)
			minString = '({:.6f} {:.6f} {:.6f})'.format(-999, -999, 0)
			maxString = '({:.6f} {:.6f} {:.6f})'.format(999, 999, cellLength*self.nCellsBetweenLevels)
			name = 'bottomBox'
			self.snappyDict.addGeometry(name, 'searchableBox', {'min':minString, 'max':maxString})
			self.snappyDict.addRefinementRegion(name, 'inside', np.array([1, level]))

		if self.viscousWake:
			self.addViscousWake()

		self.snappyDict.write(self.systemFolder)
		self.snappyDict.writeSurfaceFeatureExtractDict(self.systemFolder, 'wing.obj')

	def writeCaseFiles(self):
		super().writeCaseFiles()

		if self.turbulenceType == 'LES':
			os.remove(self.systemFolder+'fvSchemes')
			shutil.move(self.systemFolder+'fvSchemes.LES', self.systemFolder+'fvSchemes')

		if not(self.noSlipGround):
			boundaryFiles = os.listdir(self.boundaryFolder)

			# Change bottom boundary condition
			for f in boundaryFiles:
				if os.path.isfile(self.boundaryFolder + f) and f != 'pointDisplacement':
					FileHandling.changeBoundary(self.boundaryFolder + f, 'bottom', '\t\ttype    slip;\n')


		if self.WIGSim:
			boundaryFiles = os.listdir(self.boundaryFolder)

			# Create new boundary for the ground
			for f in boundaryFiles:
				if os.path.isfile(self.boundaryFolder + f):
					FileHandling.copyBoundary(self.boundaryFolder + f, 'inlet', 'ground')
					FileHandling.changeBoundary(self.boundaryFolder + f, 'ground', '\t\ttype    slip;\n')

		if True in self.motion:
			f = open(self.systemFolder+'controlDict', 'a')
			f.write('\nlibs\n(\n\t"libmyfvMotionSolvers.so"\n);\n')
			f.close()

			dynamicMeshDict = dynamicMesh.displacement(self.patchList, farFieldPatch ='inlet')
			dynamicMeshDict.write(self.constantFolder)

			for i in range(self.nrWings):
				if self.motion[i]:
					if self.nrWings == 1:
						wingName = 'wing'
					else:
						wingName = 'wing{:.0f}'.format(i+1)

					self.motionData[i].write(self.runFolder+'motion_'+wingName)

					boundaryString = '\t\ttype\ttwoDimMotion;\n\t\torigin\t({:.6f} {:.6f} {:.6f});\n\t\ttimeDataFileName\t"$FOAM_CASE/motion_'.format(self.motionData[i].origin[0], self.motionData[i].origin[1], self.motionData[i].origin[2])+wingName+'";\n\t\tvalue\tuniform (0 0 0);\n'
					FileHandling.changeBoundary(self.boundaryFolder+'pointDisplacement', wingName, boundaryString)

					FileHandling.changeBoundary(self.boundaryFolder+'U', wingName, '\t\ttype    movingWallVelocity;\n\t\tvalue    uniform (0 0 0);\n')

					if self.noSlipGround:
						FileHandling.changeBoundary(self.boundaryFolder+'U', 'bottom', '\t\ttype    movingWallVelocity;\n\t\tvalue    uniform (0 0 0);\n')

	def writeScripts(self):
		# ------ Mesh --------------------
		f = open(self.runFolder+'/mesh.sh', 'w')

		f.write('#!/bin/bash\n\n')

		if self.snappyDict.snapControls['explicitFeatureSnap'] == 'true':
			f.write('surfaceFeatureExtract\n')

		f.write('blockMesh\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.sim\n')
		f.write('mv system/decomposeParDict.mesh system/decomposeParDict\n')
		f.write('decomposePar\n')

		f.write('mpirun -np {:.0f} snappyHexMesh -overwrite -parallel\n'.format(self.nCPUs_mesh))

		f.write('reconstructParMesh -constant\n')
		f.write('rm -fr processor*\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.mesh\n')
		f.write('mv system/decomposeParDict.sim system/decomposeParDict\n')

		if len(self.topoSetList) > 0:
			f.write('topoSet\n')

		f.write('renumberMesh -overwrite\n')

		f.close()

		# ------- Simulation ---------------------
		f = open(self.runFolder + '/runSim.sh', 'w')

		f.write('#!/bin/bash\n\n')

		#if self.solver == 'simpleFoam':
		#	self.calculateBoundaryLayerHeight()
		#	f.write('applyBoundaryLayer -ybl {:.6f}\n\n'.format(self.boundaryLayerHeight))
		
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

	def setWIGSim(self, heightAboveGround):
		self.WIGSim            = True
		self.heightAboveGround = heightAboveGround

	def addMotion(self, t, x, y, omega, origin, wingIndex = 0):
		if wingIndex > self.nrWings:
			print('Too large foil index. Ignoring your request for motion!')

		if self.solver != 'pimpleDyMFoam':
			self.setSolver('pimpleDyMFoam')
			self.timeStepScheme   = 'Euler'
			self.adjustTimeStep   = True
			self.maxCo            = 10
			self.timeStepFactor   = 1

			if self.timeStepSetting == 'large':
				self.timeStepFactor *= np.sqrt(2.0)
				self.deltaT         *= np.sqrt(2.0)
				self.maxDeltaT      *= np.sqrt(2.0)
				self.maxCo          *= np.sqrt(2.0)
			elif self.timeStepSetting == 'veryLarge':
				self.timeStepFactor *= 2.0
				self.deltaT         *= 2.0
				self.maxDeltaT      *= 2.0
				self.maxCo          *= 2.0
			elif self.timeStepSetting == 'small':
				self.timeStepFactor /= np.sqrt(2.0)
				self.deltaT         /= np.sqrt(2.0)
				self.maxDeltaT      /= np.sqrt(2.0)
				self.maxCo          /= np.sqrt(2.0)
			elif self.timeStepSetting == 'verySmall':
				self.timeStepFactor /= 2.0
				self.deltaT         /= 2.0
				self.maxDeltaT      /= 2.0
				self.maxCo          /= 2.0

			self.forceOrigo = np.array([origin[0], 0, 0])

			self.purgeWrite = 0
		
		self.motion[wingIndex] = True
		self.motionData[wingIndex] = motion(t, x, y, omega, origin)

		# Write more time steps
		self.writeInterval = np.round(1/4, decimals = 6)

	def setViscousWake(self, x0, y0, z0):
		self.viscousWake = True

		self.visc_x0 = x0
		self.visc_y0 = y0
		self.visc_z0 = z0

		self.viscLengthFactor = 2
		self.viscRadiusFactor = 0.25
		self.viscExpansion    = 2

	def addViscousWake(self):
		maxLevel = self.computeCellLevel(self.viscousLength)

		radius0 = self.viscRadiusFactor*self.L
		length0 = self.viscLengthFactor*self.L
		x0 = self.visc_x0
		y0 = self.visc_y0
		z0 = self.visc_z0

		level = maxLevel
		for i in range(maxLevel):
			cellLength = self.baseSize/(2**level+1)

			name = 'viscBox{:.0f}'.format(i+1)
			length = length0*self.viscExpansion**(i)
			radius = radius0*self.viscExpansion**(i)

			point1String = '({:.6f} {:.6f} {:.6f})'.format(x0, y0, z0)
			point2String = '({:.6f} {:.6f} {:.6f})'.format(x0+length, y0, z0)
			radiusString = '{:.6f}'.format(radius)

			self.snappyDict.addGeometry(name, 'searchableCylinder', {'point1':point1String, 'point2':point2String, 'radius':radiusString})
			self.snappyDict.addRefinementRegion(name, 'inside', np.array([1, level]))

			level -= 1

		distance = radius0
		self.snappyDict.addRefinementRegion('wing', 'distance', np.array([distance, maxLevel]))

	def addRotation(self, omega):
		boundaryString = '\t\ttype rotatingWallVelocity;\n\t\torigin (0.5 0 0);\n\t\taxis (0 0 -1);\n\t\tomega {};\n'.format(omega)

		FileHandling.changeBoundary(self.boundaryFolder+'U', 'wing', boundaryString)