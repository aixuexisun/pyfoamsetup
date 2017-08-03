import numpy as np
import os
import subprocess
import shutil
import multiprocessing
import sys
import copy
from collections import OrderedDict

from pyfoamsetup.coreLibrary import *
import pyfoamsetup.coreLibrary.CaseSetup as CaseSetup

class PropellerSimulation(CaseSetup.CaseSetup):
	def __init__(self, runName, c, D, Re, J, fluid='air', rotationAxis='right', pushPull='push'):
		# Default environment settings
		if fluid == 'air':
			rho = 1.226
			nu  = 1.45e-05
		elif fluid == 'water':
			rho = 1000
			nu  = 1.19e-6

		self.D     = D
		self.r     = D/2
		self.J     = J
		self.n     = Re*nu/(np.sqrt((J*D)**2 + (0.7*self.r*2*np.pi)**2)*c)
		self.omega = 2*np.pi*self.n
		self.U_r   = 0.7*(D/2)*self.omega
		
		U = self.J*self.n*self.D
		A = np.pi*(D/2)**2

		self.rotationAxis = rotationAxis

		patchList = ['propeller']

		# Call init from base class
		self.homePath = os.path.dirname(os.path.realpath(__file__))
		super().__init__(runName, patchList, c, U, A, nu, rho, 'PropellerSimulation')

		# Reset reynolds number from input
		self.Re = Re

		# Time step settings
		self.maxDegreesPrTimeStep = 2
		self.numberOfRevolutions = 4

		self.baseSize     = 1
		self.domainWake   = 6
		self.domainFront  = 4
		self.domainWidth  = 4

		self.rotatingCylinderRadius = 0.75
		self.rotatingCylinderLength = 1

		self.setMeshSettings()
		self.nrLayers = 0

		self.setSolver('pimpleDyMFoam')
		self.adjustTimeStep = False

	def setDefaultCellLengths(self):
		super().setDefaultCellLengths()
		
		self.maxBaseSize     = 0.1*self.D # Ensure that the geometry is captured!
		self.maxSmallestSize = 0.01*self.L
		self.viscousLength   = 0.02*self.D
		
	def writeBlockMesh(self):
		blockMesh = BlockMesh.Dict()

		# Calculate minimum values for domain size
		xBack  =  self.domainWake*self.D
		xFront = -self.domainFront*self.D

		yRight =  self.domainWidth*self.D
		yLeft  = -self.domainWidth*self.D

		zHeight =  self.domainWidth*self.D
		zDepth  = -self.domainWidth*self.D

		# Calculate number of cells in each direction
		x_nrCells = np.ceil((xBack - xFront)/self.baseSize)
		y_nrCells = np.ceil((yRight - yLeft)/self.baseSize)
		z_nrCells = np.ceil((zHeight - zDepth)/self.baseSize)

		# Readjust domain size to fit nr cells
		xLength = self.baseSize*x_nrCells
		yLength = self.baseSize*y_nrCells
		zLength = self.baseSize*z_nrCells

		wakeFraction  = (self.domainWake/(self.domainWake + self.domainFront))
		frontFraction = (self.domainFront/(self.domainWake + self.domainFront))
		xFront = -xLength*frontFraction
		xBack  =  xLength*wakeFraction
	
		yRight =  yLength/2
		yLeft  = -yLength/2

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
	
		blockMesh.addBoundary('inlet', 'patch', [[0, 4, 7, 3],[3, 2, 6, 7], [4, 5, 6, 7], [0, 1, 5, 4], [0, 3, 2, 1]])

		blockMesh.addBoundary('outlet', 'patch', [[2, 6, 5, 1]])

		blockMesh.write(self.systemFolder)

	def writeMesh(self):
		self.calculateBaseSize()
		self.writeBlockMesh()

		# Add geometry
		self.snappyDict.addGeometry('propeller.obj', 'triSurfaceMesh', {'name':'propeller'})
		self.snappyDict.addRefinementSurface('propeller', self.maxRefinementLevel-1, self.maxRefinementLevel, self.nrLayers)
		self.snappyDict.addFeature('propeller.eMesh', self.maxRefinementLevel)

		self.snappyDict.addGeometry('propellerStem.obj', 'triSurfaceMesh', {'name':'propellerStem'})
		self.snappyDict.addRefinementSurface('propellerStem', self.maxRefinementLevel-3, self.maxRefinementLevel-3, 0)

		# Add cylinders
		name = 'rotatingCylinder'
		length = self.rotatingCylinderLength*self.D
		radius = self.rotatingCylinderRadius*self.D

		x0 = 0
		level = self.maxRefinementLevel-2

		point1String = '({:.6f} {:.6f} {:.6f})'.format(x0, 0, 0)
		point2String = '({:.6f} {:.6f} {:.6f})'.format(x0+length, 0, 0)
		radiusString = '{:.6f}'.format(radius)

		extraArgumentDict = OrderedDict()
		extraArgumentDict['faceType'] = 'boundary'
		extraArgumentDict['cellZone'] = name
		extraArgumentDict['faceZone'] = name
		extraArgumentDict['cellZoneInside'] = 'inside'

		self.snappyDict.addGeometry(name, 'searchableCylinder', {'point1':point1String, 'point2':point2String, 'radius':radiusString})
		self.snappyDict.addRefinementSurface(name, level, level, 0, extraArgumentDict=extraArgumentDict)
		self.snappyDict.addRefinementRegion(name, 'inside', np.array([1, level]))

		# Set up layer settings
		self.snappyDict.addLayersControls['relativeSizes']       = 'false'
		self.snappyDict.addLayersControls['finalLayerThickness'] = self.t_final
		self.snappyDict.addLayersControls['minThickness']        = 0.5*self.t_final
		self.snappyDict.addLayersControls['expansionRatio']      = self.layerExpansion

		self.snappyDict.castellatedMeshControls['locationInMesh']           = '({:.3f} {:.3f} {:.3f})'.format(-1.03*self.D, 1.04*self.D, 1.3*self.D)

		self.snappyDict.castellatedMeshControls['nCellsBetweenLevels'] = int(self.nCellsBetweenLevels)

		self.snappyDict.write(self.systemFolder)
		self.snappyDict.writeSurfaceFeatureExtractDict(self.systemFolder, 'propeller.obj')

	def writeCaseFiles(self):
		# Recalculate time stepping
		self.deltaT        = np.round(self.maxDegreesPrTimeStep/(self.n*360), decimals=8)
		self.maxDeltaT     = np.round(self.maxDegreesPrTimeStep/(self.n*360), decimals=8)
		self.endTime       = np.round(self.numberOfRevolutions/self.n, decimals=8)
		self.writeInterval = np.round(self.endTime/10, decimals = 8)

		super().writeCaseFiles()

		FileHandling.changeLine(self.constantFolder+'dynamicMeshDict', 'omega', '\t\tomega       {:.6f};'.format(self.omega))

		if self.rotationAxis == 'left':
			FileHandling.changeLine(self.constantFolder+'dynamicMeshDict', 'axis', '\t\taxis       (-1 0 0);')

		self.writePropInfo()

		createPatchDict = createPatch.Dict()
		createPatchDict.addPatch('AMI1', 'rotatingCylinder', 'AMI2')
		createPatchDict.addPatch('AMI2', 'rotatingCylinder_slave', 'AMI1')
		createPatchDict.write(self.systemFolder)

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

		f.write('createPatch -overwrite\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.mesh\n')
		f.write('mv system/decomposeParDict.sim system/decomposeParDict\n')

		f.write('renumberMesh -overwrite\n')

		if len(self.topoSetList) > 0:
			f.write('topoSet\n')

		f.close()

		# ------- Simulation ---------------------
		f = open(self.runFolder + '/runSim.sh', 'w')

		f.write('#!/bin/bash\n\n')
		
		f.write('decomposePar\n')

		if self.vilje:
			f.write('mpiexec ' + self.solver + ' -parallel\n')
		else:
			f.write('mpirun -np {:.0f} '.format(self.nCPUs) + self.solver + ' -parallel\n')

		f.write('reconstructPar\n')

		f.write('rm -fr processor*\n')

		f.close()

	def addViscousWake(self, x0, y0, z0, lengthFactor = 4, radiusFactor = 1.0, expansion = 2):
		# Ensure that mesh size is calculated
		self.turbulence = Turbulence.Properties(self.U, self.L, self.nu, self.turbulenceModel, self.turbulenceType)
		self.turbulence.calculateInitialValues()
		self.calculateBaseSize()

		maxLevel = int(np.floor(np.log(self.baseSize/self.viscousLength)/np.log(2)))
		print(maxLevel)

		radius0 = radiusFactor*self.D
		length0 = lengthFactor*self.D 

		level = maxLevel
		for i in range(maxLevel):
			cellLength = self.baseSize/(2**level+1)

			name = 'viscWake{:.0f}'.format(i+1)
			length = length0*expansion**(i)
			radius = radius0*expansion**(i)

			point1String = '({:.6f} {:.6f} {:.6f})'.format(x0, y0, z0)
			point2String = '({:.6f} {:.6f} {:.6f})'.format(x0+length, y0, z0)
			radiusString = '{:.6f}'.format(radius)

			self.snappyDict.addGeometry(name, 'searchableCylinder', {'point1':point1String, 'point2':point2String, 'radius':radiusString})
			self.snappyDict.addRefinementRegion(name, 'inside', np.array([1, level]))

			level -= 1

	def writePropInfo(self):
		f = open(self.runFolder + 'propInfo.txt', 'w')

		f.write('D     {:.6f}\n'.format(self.D))
		f.write('c     {:.6f}\n'.format(self.L))
		f.write('Re    {:.6f}\n'.format(self.Re))
		f.write('J     {:.6f}\n'.format(self.J))
		f.write('n     {:.6f}\n'.format(self.n))
		f.write('omega {:.6f}\n'.format(self.omega))
		f.write('rho   {:.6f}\n'.format(self.rho))
		f.write('U     {:.6f}\n'.format(self.U))
		f.write('U_R   {:.6f}\n'.format(self.U_r))

		f.close()

class ActuatorDisk(CaseSetup):
	def __init__(self, runName, U, D, CT, CQ=0.0, rh_factor=0.1, alpha=0, fluid='air', meshSetting='medium', vilje=False):
		# Default environment settings
		if fluid == 'air':
			rho = 1.226
			nu  = 1.45e-05
		elif fluid == 'water':
			rho = 1000
			nu  = 1.19e-6

		self.D     = D
		self.r     = D/2
		self.r_h   = rh_factor*self.r
		self.CT    = CT
		self.CQ    = CQ
		self.alpha = alpha
		
		A = np.pi*self.r**2

		patchList = []

		# Call init from base class
		super().__init__(runName, patchList, 0.5*self.r, U, A, nu, rho, vilje)

		self.Thrust = 0.5*self.A*self.CT*self.U**2
		self.Torque = 0.5*self.A*self.CQ*self.U**2*self.D

		# Essential folder paths
		self.foamPath       = os.environ['FOAM_RUN']
		self.mainRunFolder  = self.foamPath + '/PropellerSimulation'
		self.homePath       = os.path.dirname(os.path.realpath(__file__))
		self.setFolderPaths()

		self.maxBaseSize     = 0.5*self.D
		self.maxSmallestSize = 0.01*self.D 
		self.actuatorDiskLength = 0.05*self.D
		self.viscousLength   = 0.02*self.D

		# Default mesh settings
		if meshSetting == 'fine':
			self.maxBaseSize     /= np.sqrt(2)
			self.maxSmallestSize /= np.sqrt(2)
			self.viscousLength   /= np.sqrt(2)
		elif meshSetting == 'veryFine':
			self.maxBaseSize     /= 2
			self.maxSmallestSize /= 2
			self.viscousLength   /= 2
		elif meshSetting == 'coarse':
			self.maxBaseSize     *= np.sqrt(2)
			self.maxSmallestSize *= np.sqrt(2)
			self.viscousLength   *= np.sqrt(2)
		elif meshSetting == 'veryCoarse':
			self.maxBaseSize     *= 2
			self.maxSmallestSize *= 2
			self.viscousLength   *= 2

		self.baseSize     = 1
		self.domainWake   = 10
		self.domainFront  = 5
		self.domainWidth  = 5

		self.setMeshSettings()
		self.setSolver('simpleFoam')

	def calculateBaseSize(self):
		self.smallestSize = self.maxSmallestSize

		self.baseSize = self.smallestSize*2**(self.maxRefinementLevel)
		while self.baseSize > self.maxBaseSize and self.maxRefinementLevel > self.minRefinementLevel:
			self.maxRefinementLevel -= 1
			self.baseSize = self.smallestSize*2**(self.maxRefinementLevel)

	def writeBlockMesh(self):
		blockMesh = BlockMesh.Dict()

		# Calculate minimum values for domain size
		xBack  =  self.domainWake*self.D
		xFront = -self.domainFront*self.D

		yRight =  self.domainWidth*self.D
		yLeft  = -self.domainWidth*self.D

		zHeight =  self.domainWidth*self.D
		zDepth  = -self.domainWidth*self.D

		# Calculate number of cells in each direction
		x_nrCells = np.ceil((xBack - xFront)/self.baseSize)
		y_nrCells = np.ceil((yRight - yLeft)/self.baseSize)
		z_nrCells = np.ceil((zHeight - zDepth)/self.baseSize)

		# Readjust domain size to fit nr cells
		xLength = self.baseSize*x_nrCells
		yLength = self.baseSize*y_nrCells
		zLength = self.baseSize*z_nrCells

		wakeFraction  = (self.domainWake/(self.domainWake + self.domainFront))
		frontFraction = (self.domainFront/(self.domainWake + self.domainFront))
		xFront = -xLength*frontFraction
		xBack  =  xLength*wakeFraction
	
		yRight =  yLength/2
		yLeft  = -yLength/2

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
	
		blockMesh.addBoundary('inlet', 'patch', [[0, 4, 7, 3],[3, 2, 6, 7], [4, 5, 6, 7], [0, 1, 5, 4], [0, 3, 2, 1]])

		blockMesh.addBoundary('outlet', 'patch', [[2, 6, 5, 1]])

		blockMesh.write(self.systemFolder)

	def writeMesh(self):
		self.calculateBaseSize()
		self.writeBlockMesh()

		name = 'actuatorDisk'
		p1   = np.array([0, 0, 0])
		r    = self.r
		diskDir = np.array([np.cos(self.alpha), -np.sin(self.alpha), 0])
		diskThickness = self.actuatorDiskLength

		self.addPropellerActuatorDisk(name, p1, r, diskDir, diskThickness, self.Thrust, self.Torque, self.smallestSize)

		self.snappyDict.write(self.systemFolder)

	def writeCaseFiles(self):
		super().writeCaseFiles()

	def addViscousWake(self, x0, y0, z0, lengthFactor = 4, radiusFactor = 1.0, expansion = 2):
		# Ensure that mesh size is calculated
		self.calculateBaseSize()

		maxLevel = int(np.floor(np.log(self.baseSize/self.viscousLength)/np.log(2)))

		radius0 = radiusFactor*self.D
		length0 = lengthFactor*self.D 

		level = maxLevel
		for i in range(maxLevel):
			cellLength = self.baseSize/(2**level+1)

			name = 'viscWake{:.0f}'.format(i+1)
			length = length0*expansion**(i)
			radius = radius0*expansion**(i)

			point1String = '({:.6f} {:.6f} {:.6f})'.format(x0, y0, z0)
			point2String = '({:.6f} {:.6f} {:.6f})'.format(x0+length, y0, z0)
			radiusString = '{:.6f}'.format(radius)

			self.snappyDict.addGeometry(name, 'searchableCylinder', {'point1':point1String, 'point2':point2String, 'radius':radiusString})
			self.snappyDict.addRefinementRegion(name, 'inside', np.array([1, level]))

			level -= 1

	def writeScripts(self):
		# ------ Mesh --------------------
		f = open(self.runFolder+'/mesh.sh', 'w')

		f.write('#!/bin/bash\n\n')

		f.write('blockMesh\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.sim\n')
		f.write('mv system/decomposeParDict.mesh system/decomposeParDict\n')
		f.write('decomposePar\n')

		f.write('mpirun -np {:.0f} snappyHexMesh -overwrite -parallel\n'.format(self.nCPUs_mesh))

		f.write('reconstructParMesh -constant\n')
		f.write('rm -fr processor*\n')

		f.write('mv system/decomposeParDict system/decomposeParDict.mesh\n')
		f.write('mv system/decomposeParDict.sim system/decomposeParDict\n')

		f.write('createPatch -overwrite\n')

		if len(self.topoSetList) > 0:
			f.write('topoSet\n')

		f.close()

		# ------- Simulation ---------------------
		f = open(self.runFolder + '/runSim.sh', 'w')

		f.write('#!/bin/bash\n\n')
		
		f.write('decomposePar\n')

		if self.vilje:
			f.write('mpiexec ' + self.solver + ' -parallel\n')
		else:
			f.write('mpirun -np {:.0f} '.format(self.nCPUs) + self.solver + ' -parallel\n')

		f.write('reconstructPar\n')

		f.write('rm -fr processor*\n')

		f.close()

def writeRunScripts(caseNameList, propSim, folderName=''):
	# Write run script
	f = open('run.sh', 'w')
	f.write('#!/bin/bash\n\n')
	f.write('cd $FOAM_RUN/PropellerSimulation\n\n')

	for i in range(len(caseNameList)):
		f.write('cd {0}\n'.format(caseNameList[i]))
		f.write('bash mesh.sh\n')
		f.write('bash runSim.sh\n')
		f.write('cd $FOAM_RUN/TowingTank\n\n')

	f.close()
