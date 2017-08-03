import numpy as np
from collections import OrderedDict
import os

import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self):
		super().__init__()

		self.FoamFile = OrderedDict()

		self.FoamFile['version']  = 2.0
		self.FoamFile['format']   ='ascii'
		self.FoamFile['class']    ='dictionary'
		self.FoamFile['location'] ='"constant"'
		self.FoamFile['object']   ='dynamicMeshDict'

		self.mainSetup = OrderedDict()
		self.mainSetup['dynamicFvMesh'] = 'dynamicMotionSolverFvMesh'

		self.solverSetup = OrderedDict()

	def write(self, folder, ending=''):
		f = open(folder + 'dynamicMeshDict' + ending, 'w')

		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		self.writeDict(f, self.mainSetup, 'mainSetup', wrapper=False, nrTabs=0)

		f.write('\n')

		self.writeDict(f, self.solverSetup, self.mainSetup['solver']+'Coeffs')

		f.close()

class displacement(Dict):
	def __init__(self, patchList, solverType='SBRStress', diffusivityType='quadratic inverseDistance', diffusivityFactor=1, farFieldPatch = None):
		super().__init__()

		self.mainSetup['motionSolverLibs'] = '("libfvMotionSolvers.so")'
		self.mainSetup['solver']           = 'displacement'+solverType

		patchListString = '('
		for i in range(len(patchList)):
			patchListString += patchList[i]

			if i != len(patchList)-1:
				patchListString += ' '
		patchListString += ')'

		self.solverSetup['diffusivity'] = diffusivityType + ' ' + patchListString

		if farFieldPatch != None:
			self.solverSetup['interpolation'] = 'patchCorrected (' + patchListString + ' (' + farFieldPatch + '))'

class rigidBodyMotion(Dict):
	def __init__(self, body, patchList, mass, centreOfMass):
		super().__init__()

		self.mainSetup['motionSolverLibs'] = '("librigidBodyMeshMotion.so")'
		self.mainSetup['solver']           = 'rigidBodyMotion'

		self.solverSetup['report'] = 'on'
		self.solverSetup['solver'] = OrderedDict()
		self.solverSetup['solver']['type']   = 'Newmark'

		self.solverSetup['accelerationRelaxation'] = 0.4

		self.solverSetup['bodies'] = OrderedDict()

		self.solverSetup['bodies'][body] = OrderedDict()
		self.solverSetup['bodies'][body]['type']         = 'rigidBody'
		self.solverSetup['bodies'][body]['parent']       = 'root'
		self.solverSetup['bodies'][body]['centreOfMass'] = '({:.9f} {:.9f} {:.9f})'.format(centreOfMass[0], centreOfMass[1], centreOfMass[2])
		self.solverSetup['bodies'][body]['mass']         = mass

class sixDoFRigidBodyMotion(Dict):
	def __init__(self, patchList, mass, centreOfMass, momentOfInertia, L, centreOfRotation=None, accelerationRelaxation=0.4, accelerationDamping=0.9):
		super().__init__()

		self.mainSetup['motionSolverLibs'] = '("libsixDoFRigidBodyMotion.so")'
		self.mainSetup['solver']           = 'sixDoFRigidBodyMotion'

		patchListString = '('
		for i in range(len(patchList)):
			patchListString += patchList[i]

			if i != len(patchList)-1:
				patchListString += ' '
		patchListString += ')'

		self.mainSetup['diffusivity'] = 'quadratic inverseDistance 1' + patchListString

		self.solverSetup['patches'] = patchListString

		self.solverSetup['innerDistance'] = 0.1*L
		self.solverSetup['outerDistance'] = 0.4*L

		self.solverSetup['centreOfMass']    = '({:.9f} {:.9f} {:.9f})'.format(centreOfMass[0], centreOfMass[1], centreOfMass[2])
		self.solverSetup['mass']            = mass
		self.solverSetup['momentOfInertia'] = '({:.9f} {:.9f} {:.9f})'.format(momentOfInertia[0], momentOfInertia[1], momentOfInertia[2])
		self.solverSetup['rhoInf']          = 1
		self.solverSetup['report']          = 'on'
		self.solverSetup['reportToFile']    = 'on'

		self.solverSetup['value']                  = 'uniform (0 0 0)'
		self.solverSetup['accelerationRelaxation'] = accelerationRelaxation
		self.solverSetup['accelerationDamping']    = accelerationDamping

		self.solverSetup['solver'] = OrderedDict()
		self.solverSetup['solver']['type']   = 'Newmark'

		self.solverSetup['constraints'] = OrderedDict()

		self.solverSetup['constraints']['zAxis'] = OrderedDict()
		self.solverSetup['constraints']['zAxis']['sixDoFRigidBodyMotionConstraint'] = 'line'
		self.solverSetup['constraints']['zAxis']['direction']                       = '(0 0 1)'
		if type(centreOfRotation) == np.ndarray:
			self.solverSetup['constraints']['zAxis']['centreOfRotation'] = '({:.9f} {:.9f} {:.9f})'.format(centreOfRotation[0], centreOfRotation[1], centreOfRotation[2])

		self.solverSetup['constraints']['yPlane'] = OrderedDict()
		self.solverSetup['constraints']['yPlane']['sixDoFRigidBodyMotionConstraint'] = 'axis'
		self.solverSetup['constraints']['yPlane']['axis']                            = '(0 1 0)'

		'''self.solverSetup['restraints'] = OrderedDict()
		self.solverSetup['restraints']['translationDamper'] = OrderedDict()
		self.solverSetup['restraints']['translationDamper']['sixDoFRigidBodyMotionRestraint'] = 'linearDamper' 
		self.solverSetup['restraints']['translationDamper']['coeff']                          = 8596

		self.solverSetup['restraints']['rotationDamper'] = OrderedDict()
		self.solverSetup['restraints']['rotationDamper']['sixDoFRigidBodyMotionRestraint'] = 'sphericalAngularDamper'
		self.solverSetup['restraints']['rotationDamper']['coeff'] = 11586'''
		