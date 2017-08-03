from collections import OrderedDict

import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Block():
	def __init__(self, numberOfCells, topology=[0, 1, 2, 3, 4, 5, 6, 7], grading=[1, 1, 1]):
		self.topology      = topology
		self.gradingType   = 'simpleGrading'
		self.grading       = grading
		self.numberOfCells = numberOfCells
		self.type          = 'hex'

class Boundary():
	def __init__(self, name, boundaryType, faces, extraArgument):
		self.name          = name
		self.type          = boundaryType
		self.faces         = faces
		self.extraArgument = extraArgument

class Dict(FoamFile.Dict):
	def __init__(self):
		super().__init__()

		self.FoamFile = OrderedDict()
		self.FoamFile['version'] = 2.0
		self.FoamFile['format']  ='ascii'
		self.FoamFile['class']   ='dictionary'
		self.FoamFile['object']  ='blockMeshDict'

		self.convertToMeters = 1

		self.vertices   = []
		self.blocks     = []
		self.boundaries = []

	def addVertex(self, vertex):
		self.vertices.append(vertex)

	def addBlock(self, numberOfCells, topology=[0, 1, 2, 3, 4, 5, 6, 7], grading=[1, 1, 1]):
		self.blocks.append(Block(numberOfCells, topology, grading))

	def addBoundary(self, name, boundaryType, faces, extraArgument=''):
		self.boundaries.append(Boundary(name, boundaryType, faces, extraArgument))

	def write(self, folder):
		filePath = folder + 'blockMeshDict'

		f = open(filePath, 'w')

		# Write FoamFile header
		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')

		# Convert to meters line
		f.write('convertToMeters {};\n\n'.format(self.convertToMeters))

		# Vertices
		f.write('vertices\n(\n')
		for i in range(len(self.vertices)):
			x = self.vertices[i][0]
			y = self.vertices[i][1]
			z = self.vertices[i][2]

			f.write('\t({0} {1} {2})\n'.format(x, y, z))
		f.write(');\n\n')

		# Blocks
		f.write('blocks\n(\n')
		for i in range(len(self.blocks)):
			block = self.blocks[i]
			f.write('\t')
			f.write(block.type)
			f.write(' (')
			for j in range(len(block.topology)):
				index = block.topology[j]
				if j == (len(block.topology) - 1):
					f.write('{:.0f}) '.format(index))
				else:
					f.write('{:.0f} '.format(index))
			f.write('({:.0f} {:.0f} {:.0f}) '.format(block.numberOfCells[0], block.numberOfCells[1], block.numberOfCells[2]))
			f.write(block.gradingType)
			f.write(' ({0} {1} {2})\n'.format(block.grading[0], block.grading[1], block.grading[2]))
		f.write(')\n\n')

		# Edges
		f.write('edges\n(\n);\n\n')

		# Boundaries
		f.write('boundary\n(\n')
		for i in range(len(self.boundaries)):
			boundary = self.boundaries[i]

			f.write('\t{}\n'.format(boundary.name))
			f.write('\t{\n')
			f.write('\t\ttype {};\n'.format(boundary.type))
			if boundary.type == 'cyclic' or boundary.type == 'cyclicAMI':
				f.write('\t\tneighbourPatch {};\n'.format(boundary.extraArgument))
				f.write('\t\tmatchTolerance 0.01;\n')
			f.write('\t\tfaces\n\t\t(\n')
			for j in range(len(boundary.faces)):
				face = boundary.faces[j]
				f.write('\t\t\t({0} {1} {2} {3})\n'.format(face[0], face[1], face[2], face[3]))
			f.write('\t\t);\n')
			f.write('\t}\n')
		f.write(');\n')

		f.close()