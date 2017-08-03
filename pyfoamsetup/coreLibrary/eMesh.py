import numpy as np

class eMesh():
	def __init__(self, filePath):
		self.nrPoints = 0
		self.nrEdges  = 0

		self.name  = ''

		f = open(filePath, 'r')
		lines = f.readlines()
		f.close()

		nrLines = len(lines)

		foundNrPoints = False
		pointLine = 0
		edgeLine  = 0

		# Find name of object and number of points and edges
		for i in range(nrLines):
			line = lines[i].strip().split()

			if len(line) > 0:
				if line[0] == 'object':
					self.name = line[1].strip(';')
				elif line[0] == '(':
					if foundNrPoints:
						self.nrEdges = int(lines[i-1].strip())
						edgeLine = i+1
					else:
						self.nrPoints = int(lines[i-1].strip())
						foundNrPoints = True
						pointLine = i+1

		self.points = np.zeros((self.nrPoints, 3))
		self.edges  = np.zeros((self.nrEdges,  2), dtype=np.int)

		for i in range(self.nrPoints):
			line = lines[pointLine+i].strip().split()
			self.points[i, 0] = float(line[0][1:])
			self.points[i, 1] = float(line[1])
			self.points[i, 2] = float(line[2].strip(')'))

		for i in range(self.nrEdges):
			line = lines[edgeLine+i].strip().split()
			self.edges[i, 0] = int(line[0][1:])
			self.edges[i, 1] = int(line[1].strip(')'))

	def removeEdgesAbove(self, zMax=0):
		deleteIndices = np.zeros(self.nrEdges, dtype=np.int)
		nrDeleteEdges = 0

		for i in range(self.nrEdges):
			i1 = self.edges[i, 0]
			i2 = self.edges[i, 1]

			z1 = self.points[i1, 2]
			z2 = self.points[i2, 2]

			if z1 > zMax and z2 > zMax:
				deleteIndices[nrDeleteEdges] = i
				nrDeleteEdges += 1

		deleteIndices = deleteIndices[0:nrDeleteEdges]

		self.edges = np.delete(self.edges, deleteIndices, axis=0)
		self.nrEdges -= nrDeleteEdges



	def write(self, filePath):
		f = open(filePath, 'w')

		# Write header
		f.write('FoamFile\n{\n')
		f.write('\tversion    2.0;\n')
		f.write('\tformat     ascii;\n')
		f.write('\tclass      featureEdgeMesh;\n')
		f.write('\tlocation   constant/triSurface;\n')
		f.write('\tobject     {0};\n'.format(self.name))
		f.write('}\n\n\n')

		# Points
		f.write('{:.0f}\n(\n'.format(self.nrPoints))

		for i in range(self.nrPoints):
			f.write('({:.6f} {:.6f} {:.6f})\n'.format(self.points[i, 0], self.points[i, 1], self.points[i, 2]))

		f.write(')\n\n\n')

		# Edges
		f.write('{:.0f}\n(\n'.format(self.nrEdges))

		for i  in range(self.nrEdges):
			f.write('({:.0f} {:.0f})\n'.format(self.edges[i, 0], self.edges[i, 1]))

		f.write(')\n\n\n')

		f.close()
