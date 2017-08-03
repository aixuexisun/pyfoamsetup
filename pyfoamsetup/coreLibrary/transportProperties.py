import numpy as np
from collections import OrderedDict
import os
import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Dict(FoamFile.Dict):
	def __init__(self, nu, rho=1, freeSurface=False):
		super().__init__()

		self.FoamFile = OrderedDict()
		self.FoamFile['version']  = 2.0
		self.FoamFile['format']   ='ascii'
		self.FoamFile['class']    ='dictionary'
		self.FoamFile['object']   ='transportProperties;'

		self.freeSurface = freeSurface
		self.nu          = nu
		self.rho         = rho

		if self.freeSurface:
			self.water = OrderedDict()
			self.water['transportModel'] = 'Newtonian'
			self.water['nu']             = 'nu [ 0 2 -1 0 0 0 0 ]  {:.3e}'.format(self.nu[0])
			self.water['rho']            = 'rho [ 1 -3 0 0 0 0 0 ] {:.3f}'.format(self.rho[0])

			self.air = OrderedDict()
			self.air['transportModel'] = 'Newtonian'
			self.air['nu']             = 'nu [ 0 2 -1 0 0 0 0 ]  {:.3e}'.format(self.nu[1])
			self.air['rho']            = 'rho [ 1 -3 0 0 0 0 0 ] {:.3f}'.format(self.rho[1])

		else:
			self.fluid = OrderedDict()
			self.fluid['transportModel'] = 'Newtonian'
			self.fluid['nu']             = 'nu [ 0 2 -1 0 0 0 0 ]  {:.3e}'.format(self.nu)
		
	def write(self, folder, ending=''):
		f = open(folder + 'transportProperties' + ending, 'w')

		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		if self.freeSurface:
			f.write('phases (water air);\n\n')
			self.writeDict(f, self.water, 'water')
			self.writeDict(f, self.air, 'air')

			f.write('sigma           sigma [ 1 0 -2 0 0 0 0 ] 0;\n')
		else:
			self.writeDict(f, self.fluid, 'fluid', wrapper=False)


		f.close()