import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import os
import sys
import copy
import subprocess as sp
import shutil
import multiprocessing
from collections import OrderedDict

import polymesh.mesh as Mesh

from pyfoamsetup.coreLibrary import *

from pyfoamsetup.TowingTank.TowingTank import *

class HydrofoilTowingTank(TowingTank):
	def __init__(self, runName, c, T, B, A, U, freeSurface=True, symmetry=False):
		super().__init__(runName, c, T, B, A, U, freeSurface=freeSurface, symmetry=symmetry)

		# Change patch list
		self.patchList = ['wing']

		# Mesh settings
		self.domainWake_span    = 2
		self.domainWake_chord   = 40
		self.domainFront_span   = 1
		self.domainFront_chord  = 20
		self.domainWidth_chord  = 10
		self.domainWidth_span   = 2
		self.domainHeight_span  = 0.5
		self.domainHeight_chord = 5
		self.domainDepth_span   = 0.5
		self.domainDepth_chord  = 10
		
		self.freeSurfaceWidth     = 0.25*self.L
		self.freeSurfaceExpansion = 2

		self.uniformWallRefinment = True

		self.setMeshSettings()

		if not(freeSurface):
			self.endTime = 15*self.L/self.U
		else:
			self.endTime = 30*self.L/self.U

	def setDefaultCellLengths(self):
		self.maxSmallestSize   = 0.005*self.L
		self.viscousLength     = 1.1*0.04*self.L
		self.maxBaseSize       = 1.1*0.25*self.L
		self.freeSurfaceLength = 1.1*0.08*self.L

		self.cellLengthVariables = [self.maxSmallestSize, self.maxBaseSize, self.viscousLength, self.freeSurfaceLength]

	def setMeshSettings(self):
		super().setMeshSettings()
		if self.wallFunction:
			self.nrLayers = 5
			self.maxNrLayers = 10
			self.layerExpansion = 1.25
		else:
			self.maxNrLayers = 20

		if self.uniformWallRefinment:
			self.finalLayerFactor = 0.5
		else:
			self.finalLayerFactor = 1.0

		self.nCellsBetweenLevels = 3

	def writeBlockMesh(self):
		blockMesh = BlockMesh.Dict()

		# Calculate minimum values for domain size
		xBack  =  max(self.domainWake_span*self.B, self.domainWake_chord*self.L)
		xFront = -max(self.domainFront_span*self.B, self.domainFront_chord*self.L)

		yRight =  self.B/2 + max(self.domainWidth_span*self.B/2, self.domainWidth_chord*self.L)
		if self.symmetry:
			yLeft = 0
		else:
			yLeft = -yRight

		zDepth = -self.T-max(self.domainDepth_span*self.B, self.domainDepth_chord*self.L)
		
		if self.freeSurface:
			zHeight = max(self.domainHeight_span*self.B, self.domainHeight_chord*self.L)
		else:
			zHeight = max(self.domainHeight_span*self.B, self.domainHeight_chord*self.L) 

		# Calculate number of cells in each direction
		x_nrCells     = np.ceil((xBack  - xFront)/self.baseSize)
		y_nrCells     = np.ceil((yRight - yLeft)/self.baseSize)
		z_nrCells_top = np.ceil(zHeight/self.baseSize)
		z_nrCells_bot = np.ceil((-zDepth)/self.baseSize)

		z_nrCells = z_nrCells_bot + z_nrCells_top

		# Readjust domain size to fit nr cells
		xLength = self.baseSize*x_nrCells
		yLength = self.baseSize*y_nrCells

		wakeFraction  = (xBack/(-xFront + xBack))
		frontFraction = (-xFront/(-xFront + xBack))
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

	def writeMesh(self):
		self.calculateBaseSize()
		self.writeBlockMesh()

		self.snappyDict.addGeometry('wing.obj', 'triSurfaceMesh', {'name':'wing'})
		if self.uniformWallRefinment:
			self.snappyDict.addRefinementSurface('wing', self.maxRefinementLevel, self.maxRefinementLevel, self.nrLayers)
		else:
			self.snappyDict.addRefinementSurface('wing', self.maxRefinementLevel-1, self.maxRefinementLevel, self.nrLayers)
		
		self.snappyDict.addFeature('wing.eMesh', self.maxRefinementLevel)

		# Set up layer settings
		self.snappyDict.addLayersControls['relativeSizes']       = 'false'
		self.snappyDict.addLayersControls['finalLayerThickness'] = self.t_final
		self.snappyDict.addLayersControls['minThickness']        = 0.5*self.t_final
		self.snappyDict.addLayersControls['expansionRatio']      = self.layerExpansion

		self.snappyDict.castellatedMeshControls['locationInMesh']      = '({:.3f} {:.3f} {:.3f})'.format(-1.0*self.L, 0.5*self.B, -1.1*self.T)
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
		self.snappyDict.writeSurfaceFeatureExtractDict(self.systemFolder, 'wing.obj')

	def writeCaseFiles(self):
		super().writeCaseFiles()

		boundaryFiles = os.listdir(self.boundaryFolder)

		for i in range(len(boundaryFiles)):
			boundaryFile = boundaryFiles[i]
			filePath = self.boundaryFolder + boundaryFile

			if os.path.isfile(filePath):
				FileHandling.changeWord(filePath, 'ship', 'wing')

		# Change forces file
		FileHandling.changeLine(self.systemFolder+'forces', 'liftDir', '\tliftDir    (0 0 1);')
		FileHandling.changeLine(self.systemFolder+'forces', 'pitchAxis', '\tpitchAxis    (0 1 0);')

	def setViscousWake(self, x0):
		self.viscousWake = True

		self.visc_x0 = x0

		self.viscLengthFactor = 5
		self.viscWidthFactor  = 0.75
		self.viscExpansion    = 1.5

	def addViscousWake(self):
		self.viscMaxLevel = self.computeCellLevel(self.viscousLength)

		length0 = self.L*self.viscLengthFactor
		width0  = self.L*self.viscWidthFactor

		height0 = self.freeSurfaceWidth

		level = self.viscMaxLevel
		for i in range(self.viscMaxLevel):
			cellLength = self.baseSize/(2**level+1)

			name      = 'viscBox{:.0f}'.format(i+1)

			length = length0*self.viscExpansion**(i)
			width  = width0*self.viscExpansion**(i)
				
			minString = '({:.6f} {:.6f} {:.6f})'.format(self.visc_x0,          -(self.B/2 + width), -self.T - width)
			maxString = '({:.6f} {:.6f} {:.6f})'.format(self.visc_x0 + length,  (self.B/2 + width), -self.T + width)

			self.snappyDict.addGeometry(name, 'searchableBox', {'min':minString, 'max':maxString})
			self.snappyDict.addRefinementRegion(name, 'inside', np.array([level, level]))

			level -= 1

		# Outside refinement
		distance = width0

		self.snappyDict.addRefinementRegion('wing', 'distance', np.array([distance, self.viscMaxLevel]))
