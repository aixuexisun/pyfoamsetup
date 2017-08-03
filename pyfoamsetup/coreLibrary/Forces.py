import numpy as np
from collections import OrderedDict
import os
import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self, patchList, rho, U, L, A, splitPatches=False, interFoam=False, version=3):
		super().__init__()

		self.patchList = patchList
		self.nrPatches = len(patchList)
		self.rho       = rho
		self.U         = U
		self.L         = L
		self.A         = A

		self.liftDir   = [0, 1, 0]
		self.dragDir   = [1, 0, 0]
		self.pitchAxis = [0, 0, 1]
		self.CofR      = [0, 0, 0]

		self.splitPatches = splitPatches

		self.forces = OrderedDict()
		self.forces['type']               = 'forces'
		self.forces['functionObjectLibs'] = '( "libforces.so" )'

		if version == 4:
			rho = 'rho'
			U   = 'U'
		else:
			rho = 'rhoName'
			U   = 'UName'

		if interFoam:
			self.forces[rho] = 'rho'
		else:
			self.forces[rho] = 'rhoInf'

		self.forces['log']            = 'yes'

		if version == 4:
			self.forces['writeControl']  = 'timeStep'
			self.forces['writeInterval'] = 1
		else:
			self.forces['outputControl']  = 'timeStep'
			self.forces['outputInterval'] = 1

		self.forceCoeffs = OrderedDict()
		self.forceCoeffs['type']               = 'forceCoeffs'
		self.forceCoeffs['functionObjectLibs'] = '( "libforces.so" )'

		self.forceCoeffs['log'] = 'yes'
		if version == 4:
			self.forceCoeffs['writeControl']  = 'timeStep'
			self.forceCoeffs['writeInterval'] = 1
		else:
			self.forceCoeffs['outputControl']  = 'timeStep'
			self.forceCoeffs['outputInterval'] = 1
		
		if interFoam:
			self.forceCoeffs[rho] = 'rho'
		else:
			self.forceCoeffs[rho] = 'rhoInf'

	def write(self, folder):
		filePath = folder + 'forces'

		f = open(filePath, 'w')

		f.write(self.header)

		# ----- forces -------------------------------------
		if self.splitPatches and self.nrPatches > 1:
			for i in range(self.nrPatches):
				f.write('\n\nforces{:.0f}'.format(i+1)+'\n{\n')

				for k, v in self.forces.items():
					f.write('\t{0} {1};\n'.format(k, v))

				# Patches
				f.write('\tpatches (' + self.patchList[i] + ');\n')

				# Density
				f.write('\t' + 'rhoInf' + ' {};\n'.format(self.rho))

				# Origo
				f.write('\tCofR ({0} {1} {2});\n'.format(self.CofR[0], self.CofR[1], self.CofR[2]))

				f.write('}')

				# ----- force coefficients -------------------------
				f.write('\n\nforceCoeffs{:.0f}'.format(i+1)+'\n{\n')

				for k, v in self.forceCoeffs.items():
					f.write('\t{0} {1};\n'.format(k, v))

				# Patches
				f.write('\tpatches (' + self.patchList[i] + ');\n')

				# Density
				f.write('\t' + 'rhoInf' + ' {};\n'.format(self.rho))

				# Origo
				f.write('\tCofR ({0} {1} {2});\n'.format(self.CofR[0], self.CofR[1], self.CofR[2]))

				# Directions
				f.write('\tliftDir ({0} {1} {2});\n'.format(self.liftDir[0], self.liftDir[1], self.liftDir[2]))
				f.write('\tdragDir ({0} {1} {2});\n'.format(self.dragDir[0], self.dragDir[1], self.dragDir[2]))
				f.write('\tpitchAxis ({0} {1} {2});\n'.format(self.pitchAxis[0], self.pitchAxis[1], self.pitchAxis[2]))

				# Ref values
				f.write('\tmagUInf {};\n'.format(self.U))
				f.write('\tlRef {};\n'.format(self.L))
				f.write('\tAref {};\n'.format(self.A))

				f.write('}')
		else:
			f.write('\n\nforces\n{\n')

			for k, v in self.forces.items():
				f.write('\t{0} {1};\n'.format(k, v))

			# Patches
			f.write('\tpatches (')
			for patch in self.patchList:
				f.write(' '+patch)
			f.write(' );\n')

			# Density
			f.write('\t' + 'rhoInf' + ' {};\n'.format(self.rho))

			# Origo
			f.write('\tCofR ({0} {1} {2});\n'.format(self.CofR[0], self.CofR[1], self.CofR[2]))

			f.write('}')

			# ----- force coefficients -------------------------
			f.write('\n\nforceCoeffs\n{\n')

			for k, v in self.forceCoeffs.items():
				f.write('\t{0} {1};\n'.format(k, v))

			# Patches
			f.write('\tpatches (')
			for patch in self.patchList:
				f.write(' '+patch)
			f.write(' );\n')

			# Density
			f.write('\t' + 'rhoInf' + ' {};\n'.format(self.rho))

			# Origo
			f.write('\tCofR ({0} {1} {2});\n'.format(self.CofR[0], self.CofR[1], self.CofR[2]))

			# Directions
			f.write('\tliftDir ({0} {1} {2});\n'.format(self.liftDir[0], self.liftDir[1], self.liftDir[2]))
			f.write('\tdragDir ({0} {1} {2});\n'.format(self.dragDir[0], self.dragDir[1], self.dragDir[2]))
			f.write('\tpitchAxis ({0} {1} {2});\n'.format(self.pitchAxis[0], self.pitchAxis[1], self.pitchAxis[2]))

			# Ref values
			f.write('\tmagUInf {};\n'.format(self.U))
			f.write('\tlRef {};\n'.format(self.L))
			f.write('\tAref {};\n'.format(self.A))

			f.write('}')

		f.close()