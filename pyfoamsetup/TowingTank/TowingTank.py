import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import os
import sys
import copy
import shutil
import multiprocessing
from collections import OrderedDict

import polymesh.mesh as Mesh

from pyfoamsetup.coreLibrary import *
import pyfoamsetup.coreLibrary.CaseSetup as CaseSetup

class TowingTank(CaseSetup.CaseSetup):
	def __init__(self, runName, L, T, B, A, U, freeSurface=True, symmetry=False, rudder=False, defaultSettingsType='drift'):
		# Default environment settings
		rho = 999.1
		nu  = 1.14e-6
		
		self.rho_air = 1.226
		self.nu_air  = 1.45e-5

		self.rho_fullScale = 1025

		if rudder:
			patchList = ['ship', 'rudder']
		else:
			patchList = ['ship']

		self.defaultSettingsType = defaultSettingsType # waveResistance, drift
		self.setDefaultSettings()

		self.twoPartShip          = False
		self.slipWalls            = False
		self.uniformWallRefinment = self.defaultSettings['uniformWallRefinment']
		self.rudderSim            = False # Switch to generate proper mesh around the rudder

		# Extra input variables
		self.freeSurface = freeSurface
		self.symmetry    = symmetry
		self.rudder      = rudder

		# Call init from base class
		self.homePath = os.path.dirname(os.path.realpath(__file__))
		super().__init__(runName, patchList, L, U, A, nu, rho, 'TowingTank')

		self.T = T
		self.B = B

		if self.symmetry:
			self.A  = A/2
		else:
			self.A = A

		self.g  = 9.81
		self.Fr = self.U/np.sqrt(self.L*self.g)

		# Dimensions of the domain
		self.baseSize         = 1
		self.domainWake       = 6
		self.domainFront      = 2
		self.domainWidth      = 2
		self.domainHeight     = 1
		self.domainHeight_air = 2
		self.domainDepth      = 2
		
		self.freeSurfaceWidth     = 0.01*self.L
		self.freeSurfaceExpansion = 2
		self.freeSurfaceRefinementMethod = 'snapAfter' # 'snapFirst' or 'snapAfter' 
	
		self.schemeType     = 'robustAndAccurate' #robustAndAccurate, accurate
		self.timeStepScheme = 'Euler' # Euler, CrankNicolson, backwards
		self.wallFunction   = True

		self.adjustTimeStep = True

		if freeSurface:
			self.setSolver('interFoam')
		else:
			self.setSolver('pimpleFoam')	

		self.setTurbulenceModel('kOmegaSST', turbulenceType='RAS') # kOmegaSST, realizableKE

		self.nNonOrthogonalCorrectors = 0

		if not(freeSurface):
			self.endTime       = 5.0*self.L/self.U
		else:
			self.endTime       = 15*self.L/self.U
		
		self.writeInterval = self.endTime/10
		
		# Options
		self.freeHeaveAndPitch = False
		self.viscousWake       = False
		self.kelvinWake        = False
		self.actuatorDisk      = False
		self.airResistanceSim  = False

	def setMeshSettings(self, yPlus=None):
		super().setMeshSettings()
		if self.wallFunction:
			self.nrLayers = 5
			self.maxNrLayers = 8
			self.layerExpansion = 1.25
		else:
			self.maxNrLayers = 15

		if self.uniformWallRefinment:
			self.finalLayerFactor = 0.75
		else:
			self.finalLayerFactor = 1.0
		
		self.nrRudderLayers = 0
		self.nCellsBetweenLevels = 5

		if yPlus != None:
			self.yPlus = yPlus

	def setDefaultCellLengths(self):
		super().setDefaultCellLengths()

		self.maxSmallestSize   = 0.001*self.L
		self.maxBaseSize       = 0.1*self.L
		self.maxRudderSize     = 1.1*0.001*self.L
		self.kelvinWakeLength  = 1.1*0.008*self.L
		self.viscousLength     = 1.1*0.004*self.L
		self.freeSurfaceLength = 0.004*self.L

	def changeCellLengths(self, cellLengthFactor):
		super().changeCellLengths(cellLengthFactor)

		self.maxRudderSize     *= cellLengthFactor
		self.kelvinWakeLength  *= cellLengthFactor
		self.freeSurfaceLength *= cellLengthFactor

	def setDefaultTimeStep(self):
		self.deltaT    = 0.005*self.L/self.U 
		self.maxDeltaT = 0.005*self.L/self.U

		if self.freeSurface:
			self.maxCo = 10
		else:
			self.maxCo = 5

		self.maxAlphaCo = 5

	def setTwoPartShip(self):
		self.twoPartShip = True
		if self.rudder:
			self.patchList = ['deck', 'hull', 'rudder']
		else:
			self.patchList = ['deck', 'hull']

	def setFreeHeaveAndPitch(self, volume, x_b, z_b=0.0):
		if self.solver == 'interFoam':
			self.freeHeaveAndPitch = True
			self.solver = 'interDyMFoam'

			self.volume = volume
			self.mass   = self.volume*self.rho

			# Correct mass and volume if half the model is used
			if self.symmetry:
				self.mass   /= 2
				self.volume /= 2

			self.centreOfMass = np.array([x_b, 0, z_b])

			# Estimate moment of intertia by assuming box shape
			self.I = (self.mass/12)*np.array([self.B**2+self.T**2, self.L**2+self.T**2, self.L**2+self.B**2])

			self.centreOfRotation = None

		else:
			print('Solver is not interFoam. Why do you want free heave and pitch? Ignoring your request...')

	def setAirResistanceSim(self):
		self.airResistanceSim = True

		if not(self.freeSurface):
			self.rho = self.rho_air
			self.nu  = self.nu_air

			self.Re  = self.L*self.U/self.nu

		self.setMeshSettings()

		if not(self.freeSurface):
			self.nCellsBetweenLevels = 3

	def writeMesh(self):
		if self.rudderSim:
			self.maxSmallestSize = self.maxRudderSize

		self.calculateBaseSize()
		self.writeBlockMesh()

		if self.twoPartShip:
			# Initialize snappy hex mesh with ship
			self.snappyDict.addGeometry('hull.obj', 'triSurfaceMesh', {'name':'hull'})

			if self.uniformWallRefinment:
				self.snappyDict.addRefinementSurface('hull', self.maxRefinementLevel, self.maxRefinementLevel, self.nrLayers)
			else:
				self.snappyDict.addRefinementSurface('hull', self.maxRefinementLevel-1, self.maxRefinementLevel, self.nrLayers)

			# Initialize snappy hex mesh with ship
			self.snappyDict.addGeometry('deck.obj', 'triSurfaceMesh', {'name':'deck'})
			self.snappyDict.addRefinementSurface('deck', self.maxRefinementLevel-2, self.maxRefinementLevel-2, 0)
		else:
			# Initialize snappy hex mesh with ship
			self.snappyDict.addGeometry('ship.obj', 'triSurfaceMesh', {'name':'ship'})

			if self.uniformWallRefinment:
				self.snappyDict.addRefinementSurface('ship', self.maxRefinementLevel, self.maxRefinementLevel, self.nrLayers)
			else:
				self.snappyDict.addRefinementSurface('ship', self.maxRefinementLevel-1, self.maxRefinementLevel, self.nrLayers)

		if self.rudder:
			minRudderLevel = self.computeCellLevel(self.maxRudderSize)
			rudderLevel = max(minRudderLevel, self.maxRefinementLevel)
			self.snappyDict.addGeometry('rudder.obj', 'triSurfaceMesh', {'name':'rudder'})

			self.snappyDict.addRefinementSurface('rudder', rudderLevel, rudderLevel, int(self.nrRudderLayers))
			
		self.snappyDict.addFeature('ship.eMesh', self.maxRefinementLevel)

		# Set up layer settings
		self.snappyDict.addLayersControls['relativeSizes']       = 'false'
		self.snappyDict.addLayersControls['finalLayerThickness'] = self.t_final
		self.snappyDict.addLayersControls['minThickness']        = 0.5*self.t_final
		self.snappyDict.addLayersControls['expansionRatio']      = self.layerExpansion

		if self.airResistanceSim and not(self.freeSurface):
			self.snappyDict.castellatedMeshControls['locationInMesh']      = '({:.3f} {:.3f} {:.3f})'.format(-0.5*self.L, 1.5*self.B, 1*self.T)
		else:
			self.snappyDict.castellatedMeshControls['locationInMesh']      = '({:.3f} {:.3f} {:.3f})'.format(-0.5*self.L, 1.5*self.B, -2*self.T)
		self.snappyDict.castellatedMeshControls['nCellsBetweenLevels'] = int(self.nCellsBetweenLevels)

		if self.viscousWake:
			self.addViscousWake()

		if self.actuatorDisk:
			self.addActuatorDisk()

		if self.freeSurface:
			self.writeFreeSurfaceRefinement()

			if self.freeSurfaceRefinementMethod == 'snapAfter':
				self.snappyDict.actionControl['snap']      = 'false'
				self.snappyDict.actionControl['addLayers'] = 'false'

				snappySnap = copy.deepcopy(self.snappyDict)
				snappySnap.actionControl['castellatedMesh'] = 'false'
				snappySnap.actionControl['snap']            = 'true'
				snappySnap.actionControl['addLayers']       = 'true'
			
				snappySnap.write(self.systemFolder, ending='.snap')
		
		self.snappyDict.write(self.systemFolder)
		self.snappyDict.writeSurfaceFeatureExtractDict(self.systemFolder, 'ship.obj')

		if self.kelvinWake:
			self.writeKelvinWake()

	def writeBlockMesh(self):
		blockMesh = BlockMesh.Dict()

		# Calculate minimum values for domain size
		xBack  =  self.domainWake*self.L
		xFront = -self.domainFront*self.L

		yRight =  self.domainWidth*self.L
		if self.symmetry:
			yLeft = 0
		else:
			yLeft = -self.domainWidth*self.L

		if self.airResistanceSim and not(self.freeSurface):
			zDepth = 0
		else:
			zDepth = -self.domainDepth*self.L
		
		if self.freeSurface:
			zHeight = self.domainHeight*self.L
		elif self.airResistanceSim:
			zHeight = self.domainHeight_air*self.L
		else:
			zHeight = 0 

		# Calculate number of cells in each direction
		x_nrCells = np.ceil((xBack - xFront)/self.baseSize)
		y_nrCells = np.ceil((yRight - yLeft)/self.baseSize)

		if self.freeSurface or self.airResistanceSim:
			z_nrCells_top = np.ceil(zHeight/self.baseSize)
		else:
			z_nrCells_top = 0

		z_nrCells_bot = np.ceil((-zDepth)/self.baseSize)

		z_nrCells = z_nrCells_bot + z_nrCells_top

		# Readjust domain size to fit nr cells
		xLength = self.baseSize*x_nrCells
		yLength = self.baseSize*y_nrCells

		wakeFraction  = (self.domainWake/(self.domainWake + self.domainFront))
		frontFraction = (self.domainFront/(self.domainWake + self.domainFront))
		xFront = -xLength*frontFraction
		xBack  =  xLength*wakeFraction
		if self.symmetry:
			yRight = yLength
			yLeft  = 0
		else:
			yRight =  yLength/2
			yLeft  = -yLength/2

		zDepth = -self.baseSize*z_nrCells_bot

		if self.freeSurface:
			zHeight = self.baseSize*z_nrCells_top

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

		if self.airResistanceSim and not(self.freeSurface):
			if self.symmetry:
				blockMesh.addBoundary('inlet', 'patch', [[0, 4, 7, 3],[3, 2, 6, 7], [4, 5, 6, 7]])
				blockMesh.addBoundary('midPlane', 'symmetryPlane', [[0, 1, 5, 4]])
			else:
				blockMesh.addBoundary('inlet', 'patch', [[0, 4, 7, 3],[3, 2, 6, 7], [4, 5, 6, 7], [0, 1, 5, 4]])
		else:
			if self.symmetry:
				blockMesh.addBoundary('inlet', 'patch', [[0, 4, 7, 3],[3, 2, 6, 7], [0, 3, 2, 1]])
				blockMesh.addBoundary('midPlane', 'symmetryPlane', [[0, 1, 5, 4]])
			else:
				blockMesh.addBoundary('inlet', 'patch', [[0, 4, 7, 3],[3, 2, 6, 7], [0, 3, 2, 1], [0, 1, 5, 4]])

		blockMesh.addBoundary('outlet', 'patch', [[2, 6, 5, 1]])

		if self.freeSurface:
			blockMesh.addBoundary('top', 'patch', [[4, 5, 6, 7]])
		elif self.airResistanceSim and not(self.freeSurface):
			blockMesh.addBoundary('top', 'symmetryPlane', [[0, 3, 2, 1]])
		else:
			blockMesh.addBoundary('top', 'symmetryPlane', [[4, 5, 6, 7]])

		blockMesh.write(self.systemFolder)

	def writeFreeSurfaceRefinement(self):
		self.freeSurfaceRefinement = self.computeCellLevel(self.freeSurfaceLength)

		maxBoxWidth = self.freeSurfaceWidth*self.freeSurfaceExpansion**(self.freeSurfaceRefinement)

		for i in range(self.freeSurfaceRefinement):
			cellLength     = self.baseSize/(2**i)
			cellVolume     = cellLength**3

			boxWidth   = max(maxBoxWidth/self.freeSurfaceExpansion**i, 4*cellLength)

			# Write files
			f = open(self.systemFolder+'topoSetDict.{:.0f}'.format(i), 'w')

			# Write header
			f.write('FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\n\tclass\tdictionary;\n\tlocation\t"system";\n\tobject\ttopoSetDict;\n}\n')

			f.write('\nactions\n(\n')
	
			f.write('\t{\n')
			f.write('\t\tname   c0;\n')
			f.write('\t\ttype   cellSet;\n')
			f.write('\t\taction new;\n')
			f.write('\t\tsource boxToCell;\n')
			f.write('\t\tsourceInfo\n')
			f.write('\t\t{\n')
			f.write('\t\t\tbox (-999 -999 {:.6f}) (999 999 {:.8f});\n'.format(-boxWidth/2, boxWidth/2))
			f.write('\t\t}\n\t}\n')

			f.write('\t{\n')
			f.write('\t\tname   c0;\n')
			f.write('\t\ttype   cellSet;\n')
			f.write('\t\taction delete;\n')
			f.write('\t\tsource fieldToCell;\n')
			f.write('\t\tsourceInfo\n')
			f.write('\t\t{\n')
			if self.version == 3:
				f.write('\t\t\tfieldName cellVolume;\n')
			else:
				f.write('\t\t\tfield V;\n')
			f.write('\t\t\tmin       0;\n')
			f.write('\t\t\tmax       {:.8f};\n'.format(cellVolume*0.75))
			f.write('\t\t}\n\t}\n')

			f.write(');\n')

			f.close()

	def writeCaseFiles(self):
		if self.airResistanceSim:
			if self.twoPartShip:
				self.splitPatches = True

		super().writeCaseFiles()

		# Initial files
		if self.freeSurface:
			os.remove(self.boundaryFolder+'p')

		else:
			boundaryFiles = os.listdir(self.boundaryFolder)

			# Change top boundary condition
			for f in boundaryFiles:
				if os.path.isfile(self.boundaryFolder + f):
					FileHandling.changeBoundary(self.boundaryFolder + f, 'top', '\t\ttype    symmetryPlane;\n')

			# Change velocity outlet
			FileHandling.changeBoundary(self.boundaryFolder+'U', 'outlet', '\t\ttype    zeroGradient;\n')

			# Remove free surface transport properties
			os.remove(self.runFolder+'constant/g')
			os.remove(self.boundaryFolder+'alpha.water')
			os.remove(self.boundaryFolder+'p_rgh')

		# Remove midPlane patch if no symmetry
		if not(self.symmetry):
			boundaryFiles = os.listdir(self.boundaryFolder)
			for f in boundaryFiles:
				if os.path.isfile(self.boundaryFolder + f):
					FileHandling.deleteBoundary(self.boundaryFolder + f, 'midPlane')

		# Adjust U boundary conditions if the ship is free in heave and pitch
		if self.freeHeaveAndPitch:
			dynamicMeshDict = dynamicMesh.sixDoFRigidBodyMotion(self.patchList, self.mass, self.centreOfMass, self.I, self.L, centreOfRotation=self.centreOfRotation)
			dynamicMeshDict.write(self.constantFolder)
			FileHandling.changeBoundary(self.boundaryFolder+'U', 'ship', '\t\ttype    movingWallVelocity;\n\t\tvalue    uniform (0 0 0);\n')

		if self.rudder:
			boundaryFiles = os.listdir(self.boundaryFolder)

			for i in range(len(boundaryFiles)):
				boundaryFile = boundaryFiles[i]
				filePath = self.boundaryFolder + boundaryFile

				if os.path.isfile(filePath):
					FileHandling.copyBoundary(filePath, 'ship', 'rudder')

		# Adjust initial conditions if the ship is split in two
		if self.twoPartShip:
			boundaryFiles = os.listdir(self.boundaryFolder)

			for i in range(len(boundaryFiles)):
				boundaryFile = boundaryFiles[i]
				filePath = self.boundaryFolder + boundaryFile

				if os.path.isfile(filePath):
					FileHandling.changeWord(filePath, 'ship', 'hull')
					FileHandling.copyBoundary(filePath, 'hull', 'deck')

		if self.actuatorDisk:
			createPatchDict = createPatch.Dict()
			createPatchDict.addPatch('AMI1', 'actuatorDisk', 'AMI2')
			createPatchDict.addPatch('AMI2', 'actuatorDisk_slave', 'AMI1')
			createPatchDict.write(self.systemFolder)

		self.writeShipInfo()

	def writeScripts(self):
		# Mesh
		f = open(self.runFolder+'/mesh.sh', 'w')

		f.write('#!/bin/bash\n\n')

		if self.snappyDict.snapControls['explicitFeatureSnap'] == 'true' and not(self.twoPartShip):
			f.write('surfaceFeatureExtract\n')

		f.write('blockMesh\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.sim\n')
		f.write('mv system/decomposeParDict.mesh system/decomposeParDict\n')
		f.write('decomposePar\n')

		f.write('mpirun -np {:.0f} snappyHexMesh -overwrite -parallel\n'.format(self.nCPUs_mesh))

		f.write('reconstructParMesh -constant\n')
		f.write('rm -fr processor*\n')

		# Add free surface refinement if necessary
		if self.freeSurface:
			f.write('\nCOUNTER=0\n')
			f.write('while [ $COUNTER -lt {:.0f} ]; do\n'.format(self.freeSurfaceRefinement))
			if self.version == 3:
				f.write('\tcellVolume\n')
			else:
				f.write('\tpostProcess -func writeCellVolumes\n')

			f.write('\ttopoSet -dict system/topoSetDict.$COUNTER\n')
			f.write('\trefineMesh -overwrite\n')
			if self.version == 3:
				f.write('rm 0/cellVolume\n')
			else:
				f.write('rm 0/V\n')

			f.write('\tlet COUNTER=COUNTER+1\n')
			f.write('done\n\n')

			if self.freeSurfaceRefinementMethod == 'snapAfter':
				# Run snapping and layer generation
				f.write('mv system/snappyHexMeshDict system/snappyHexMeshDict.0\n')
				f.write('mv system/snappyHexMeshDict.snap system/snappyHexMeshDict\n')

				f.write('decomposePar\n')
				f.write('mpirun -np {:.0f} snappyHexMesh -overwrite -parallel\n'.format(self.nCPUs_mesh))

				f.write('reconstructParMesh -constant\n')
				f.write('rm -fr processor*\n')

		if self.actuatorDisk:
			f.write('createPatch -overwrite\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.mesh\n')
		f.write('mv system/decomposeParDict.sim system/decomposeParDict\n')

		f.write('renumberMesh -overwrite\n')

		if len(self.topoSetList) > 0:
			f.write('topoSet\n')

		f.close()

		# Simulation
		f = open(self.runFolder + '/runSim.sh', 'w')

		f.write('#!/bin/bash\n\n')

		if self.freeSurface:
			f.write('setFields\n\n')
		elif self.solver == 'simpleFoam':
			self.calculateBoundaryLayerHeight()
			f.write('applyBoundaryLayer -ybl {:.6f}\n\n'.format(self.boundaryLayerHeight))
		
		f.write('decomposePar\n')

		if self.superComputer:
			if self.freeHeaveAndPitch:
				f.write('mpiexec interDyMFoam -parallel\n')
			else:
				f.write('mpiexec ' + self.solver + ' -parallel\n')
		else:
			if self.freeHeaveAndPitch:
				f.write('mpirun -np {:.0f} interDyMFoam -parallel\n'.format(self.nCPUs))
			else:
				f.write('mpirun -np {:.0f} '.format(self.nCPUs) + self.solver + ' -parallel\n')

		f.write('reconstructPar\n')

		f.write('rm -fr processor*\n')

		f.close()

	def setViscousWake(self, x0, y0):
		self.viscousWake = True

		self.visc_x0 = x0
		self.visc_y0 = y0

		self.viscLengthFactor = 0.25
		self.viscWidthFactor  = 0.1
		self.viscDepthFactor  = 1.0
		self.viscExpansion    = 1.5

	def addViscousWake(self):
		self.viscMaxLevel = self.computeCellLevel(self.viscousLength)
		cellLength = self.baseSize/(2**self.viscMaxLevel)

		length0 = self.L*self.viscLengthFactor
		width0  = self.L*self.viscWidthFactor
		depth0  = self.T*self.viscDepthFactor
		if self.airResistanceSim:
			height0 = 0.25*self.L
		else:
			height0 = self.freeSurfaceWidth

		level = self.viscMaxLevel
		for i in range(self.viscMaxLevel):
			cellLength = self.baseSize/(2**level+1)

			name = 'viscBox{:.0f}'.format(i+1)
			length = length0*self.viscExpansion**(i)

			if i == 0:
				width  = width0
				depth  = depth0
				height = height0 
			else:
				width  += 2*cellLength*self.nCellsBetweenLevels

			distance = width/2
				
			minString = '({:.6f} {:.6f} {:.6f})'.format(self.visc_x0,          self.visc_y0 - 1e-6, -depth)
			maxString = '({:.6f} {:.6f} {:.6f})'.format(self.visc_x0 + length, self.visc_y0 + 1e-6,  height)

			self.snappyDict.addGeometry(name, 'searchableBox', {'min':minString, 'max':maxString})
			self.snappyDict.addRefinementRegion(name, 'distance', np.array([distance, level]))

			level -= 1

		# Outside refinement
		distance = width0/2

		if self.twoPartShip:
			self.snappyDict.addRefinementRegion('hull', 'distance', np.array([distance, self.viscMaxLevel]))
		else:
			self.snappyDict.addRefinementRegion('ship', 'distance', np.array([distance, self.viscMaxLevel]))

	def setActuatorDisk(self, p1, diskDir, radius, Thrust, Torque):
		self.actuatorDisk = True

		self.prop_p1 = p1
		self.prop_diskDir = diskDir
		self.prop_radius = radius
		self.prop_Thrust = Thrust
		self.prop_Torque = Torque

	def addActuatorDisk(self):
		# Add cylinders
		name = 'actuatorDisk'
		level = self.computeCellLevel(self.smallestSize)

		p1 = self.prop_p1
		diskLength = max(2*self.smallestSize, 0.1*self.prop_radius)
		p2 = self.prop_p1 + diskLength*self.prop_diskDir

		point1String = '({:.6f} {:.6f} {:.6f})'.format(p1[0], p1[1], p1[2])
		point2String = '({:.6f} {:.6f} {:.6f})'.format(p2[0], p2[1], p2[2])
		radiusString = '{:.6f}'.format(self.prop_radius)

		extraArgumentDict = OrderedDict()
		extraArgumentDict['faceType'] = 'boundary'
		extraArgumentDict['cellZone'] = name
		extraArgumentDict['faceZone'] = name
		extraArgumentDict['cellSet']  = name
		extraArgumentDict['cellZoneInside'] = 'inside'

		self.snappyDict.addGeometry(name, 'searchableCylinder', {'point1':point1String, 'point2':point2String, 'radius':radiusString})
		self.snappyDict.addRefinementSurface(name, level, level, 0, extraArgumentDict=extraArgumentDict)
		self.snappyDict.addRefinementRegion(name, 'inside', np.array([1, level]))

		actuatorDisk = fvOptions.goldsteinActuationDisk('propeller', self.prop_Thrust, self.prop_Torque, p1, p2, self.prop_radius, self.prop_radius*0.2, 'actuatorDisk')		
		cellSet      = topoSet.cylinderToCell('actuatorDisk', p1, p2, self.prop_radius)

		self.fvOptionsList.append(actuatorDisk)
		self.topoSetList.append(cellSet)

	def setKelvinWake(self):
		self.kelvinWake = True

		self.kelvinLengthFactor    = 0.5
		self.kelvinThicknessFactor = 0.25
		self.kelvinDepthFactor     = 1.0
		self.kelvinExpansion       = 1.5

	def writeKelvinWake(self):
		self.kelvinMaxLevel = self.computeCellLevel(self.kelvinWakeLength)

		thickness0 = self.L*self.kelvinThicknessFactor
		length0    = self.L+self.L*self.kelvinLengthFactor
		depth0     = self.kelvinDepthFactor*self.T

		cellLength = self.baseSize/(2**self.kelvinMaxLevel)
		height0    = self.freeSurfaceWidth + self.nCellsBetweenLevels*cellLength
		
		angle = 19.47*np.pi/180

		face_verts = np.array([0, 1, 2, 3, 
				               0, 4, 5, 1,
				               1, 5, 6, 2,
				               6, 7, 3, 2,
				               3, 7, 4, 0,
				               4, 7, 6, 5])

		face_nrVerts = np.array([4, 4, 4, 4, 4, 4])

		level = self.kelvinMaxLevel

		length    = np.zeros(level)
		thickness = np.zeros(level)
		yWake     = np.zeros(level)

		for l in range(self.kelvinMaxLevel):
			cellLength = self.baseSize/(2**level+1)

			name = 'kelvinWake{:.0f}'.format(l+1)

			if l == 0:
				depth        = depth0
				height       = height0
				thickness[l] = thickness0
			else:
				depth        += 1.1*cellLength*self.nCellsBetweenLevels
				height       += 1.1*cellLength*self.nCellsBetweenLevels
				thickness[l]  = thickness[l-1] + cellLength*self.nCellsBetweenLevels

			length[l] = length0*self.kelvinExpansion**(l)
			yWake[l]  = thickness[l]/2 + (length[l]+thickness[l])*np.tan(angle)

			if l != 0:
				cellExpansionsLength = cellLength*self.nCellsBetweenLevels
				geometryLength       = (length[l-1]+thickness[l])*np.tan(angle) - (length[l-1]+thickness[l-1])*np.tan(angle)

				if cellExpansionsLength > geometryLength:
					yWake[l] += cellExpansionsLength-geometryLength 

			verts = np.array([[-thickness[l],  thickness[l]/2, height],
							  [-thickness[l], -thickness[l]/2, height],
							  [length[l],     -yWake[l],       height],
							  [length[l],      yWake[l],       height],
							  [-thickness[l],  thickness[l]/2, -depth],
							  [-thickness[l], -thickness[l]/2, -depth],
							  [length[l],     -yWake[l],       -depth],
							  [length[l],      yWake[l],       -depth]])

			mesh = Mesh.Mesh(verts, face_verts, face_nrVerts, simple=True)
			mesh.exportObj(self.geometryFolder+name+'.obj')

			f = open(self.systemFolder+'snappyHexMeshDict', 'r')
			lines = f.readlines()
			f.close()

			f = open(self.systemFolder+'snappyHexMeshDict', 'w')

			i = 0
			while i < len(lines):
				line = lines[i].strip().split()
				addGeometry = False
				addRegion = False

				if len(line) > 0:
					if line[0] == 'geometry':
						addGeometry = True
					elif line[0] =='refinementRegions':
						addRegion = True

				if addGeometry:
					f.write(lines[i])
					f.write(lines[i+1])
					f.write('\t' + name + '.obj' + '\n\t{\n\t\ttype triSurfaceMesh;\n')
					f.write('\t\tname '+ name + ';\n')
					f.write('\t}\n\n')

					i += 2
				elif addRegion:
					f.write(lines[i])
					f.write(lines[i+1])
					f.write('\t\t'+name+'\n\t\t{\n')
					f.write('\t\t\tmode inside;\n')
					f.write('\t\t\tlevels ((1.0 {:.0f}));\n'.format(level))
					f.write('\t\t}\n\n')

					i += 2
				else:
					f.write(lines[i])
					i += 1

			level -= 1 

	def writeShipInfo(self):
		f = open(self.runFolder + 'shipInfo.txt', 'w')

		f.write('rho {:.6f}\n'.format(self.rho))
		f.write('nu  {:.6f}\n'.format(self.nu))
		f.write('U   {:.6f}\n'.format(self.U))
		f.write('Re  {:.6f}\n'.format(self.Re))
		f.write('Fr  {:.6f}\n'.format(self.Fr))

		if self.actuatorDisk:
			f.write('Thrust  {:.6f}\n'.format(self.prop_Thrust))

		f.close()

	def setDefaultSettings(self):
		self.defaultSettings = {}

		if self.defaultSettingsType == 'waveResistance':
			self.defaultSettings['uniformWallRefinment'] = False
		else:
			self.defaultSettings['uniformWallRefinment'] = True