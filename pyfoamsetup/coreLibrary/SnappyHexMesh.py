import numpy as np
from collections import OrderedDict
import os
import pyfoamsetup.coreLibrary.FoamFile as FoamFile

class Geometry():
	def __init__(self, name, geometryType, argumentDict):
		self.name         = name
		self.geometryType = geometryType
		self.argumentDict = argumentDict

class RefinementSurface():
	def __init__(self, name, lmin, lmax, nrLayers, extraArgumentDict=None):
		self.name     = name
		self.lmin     = lmin
		self.lmax     = lmax
		self.nrLayers = nrLayers
		self.extraArgumentDict = extraArgumentDict

class Feature():
	def __init__(self, fileName, level):
		self.fileName = fileName
		self.level    = level

class RefinementRegion():
	def __init__(self, name, mode, levels):
		self.name   = name
		self.mode   = mode
		self.levels = levels

class Dict(FoamFile.Dict):
	def __init__(self):
		super().__init__()
		self.FoamFile = OrderedDict()
		self.FoamFile['version'] = 2.0
		self.FoamFile['format']  ='ascii'
		self.FoamFile['class']   ='dictionary'
		self.FoamFile['object']  ='snappyHexMeshDict'

		self.FoamFileFeature = OrderedDict()
		self.FoamFileFeature['version'] = 2.0
		self.FoamFileFeature['format']  ='ascii'
		self.FoamFileFeature['class']   ='dictionary'
		self.FoamFileFeature['object']  ='surfaceFeatureExtractDict'

		self.actionControl = OrderedDict()
		self.actionControl['castellatedMesh'] = 'true'
		self.actionControl['snap']            = 'true'
		self.actionControl['addLayers']       = 'true'

		self.castellatedMeshControls = OrderedDict()
		self.castellatedMeshControls['maxLocalCells']              = 1000000000
		self.castellatedMeshControls['maxGlobalCells']             = 2000000000
		self.castellatedMeshControls['minRefinementCells']         = 10
		self.castellatedMeshControls['maxLoadUnbalance:']          = 0.10
		self.castellatedMeshControls['nCellsBetweenLevels']        = 10
		self.castellatedMeshControls['resolveFeatureAngle']        = 60
		self.castellatedMeshControls['locationInMesh']             = '(-2.02 3.001 -2.03)'
		self.castellatedMeshControls['allowFreeStandingZoneFaces'] = 'true'

		self.snapControls = OrderedDict()
		self.snapControls['nSmoothPatch']           = 3
		self.snapControls['tolerance']              = 4.0
		self.snapControls['nSolveIter']             = 50
		self.snapControls['nRelaxIter']             = 5
		self.snapControls['nSmoothPatch']           = 3
		self.snapControls['nSmoothInternal']        = 3
		self.snapControls['detectNearSurfaceSnap']  = 'true'
		self.snapControls['nFeatureSnapIter']       = 10
		self.snapControls['implicitFeatureSnap']    = 'true'
		self.snapControls['explicitFeatureSnap']    = 'true'
		self.snapControls['multiRegionFeatureSnap'] = 'true'

		self.addLayersControls = OrderedDict()
		self.addLayersControls['relativeSizes']             = 'true'
		self.addLayersControls['finalLayerThickness']       = 0.5
		self.addLayersControls['expansionRatio']            = 1.2
		self.addLayersControls['minThickness']              = 0.01
		self.addLayersControls['nGrow']                     = 0
		self.addLayersControls['featureAngle']              = 60
		self.addLayersControls['nRelaxIter']                = 5
		self.addLayersControls['nSmoothSurfaceNormals']     = 1
		self.addLayersControls['nSmoothNormals']            = 3
		self.addLayersControls['nSmoothThickness']          = 10
		self.addLayersControls['maxFaceThicknessRatio']     = 0.5
		self.addLayersControls['maxThicknessToMedialRatio'] = 0.3
		self.addLayersControls['minMedianAxisAngle']        = 90
		self.addLayersControls['nBufferCellsNoExtrude']     = 0
		self.addLayersControls['nLayerIter']                = 50
		self.addLayersControls['nMedialAxisIter']           = 10

		self.meshQualityControls = OrderedDict()
		self.meshQualityControls['maxNonOrtho']         =  70 # Default 65
		self.meshQualityControls['maxBoundarySkewness'] =  20
		self.meshQualityControls['maxInternalSkewness'] =  4
		self.meshQualityControls['maxConcave']          =  80
		self.meshQualityControls['minVol']              =  1e-13
		self.meshQualityControls['minTetQuality']       = -1e30
		self.meshQualityControls['minArea']             = -1
		self.meshQualityControls['minTwist']            =  0.01 
		self.meshQualityControls['minDeterminant']      =  0.001
		self.meshQualityControls['minFaceWeight']       =  0.0
		self.meshQualityControls['minVolRatio']         =  0.01
		self.meshQualityControls['minTriangleTwist']    = -1
		self.meshQualityControls['nSmoothScale']        =  4
		self.meshQualityControls['errorReduction']      =  0.75

		self.meshQualityControls['relaxed'] = OrderedDict()
		self.meshQualityControls['relaxed']['maxNonOrtho'] = 80

		self.finalOptions = OrderedDict()
		self.finalOptions['mergeTolerance'] = 1e-6

		self.geometryList          = []
		self.refinementSurfaceList = []
		self.featureList           = []
		self.refinementRegionList  = []

	def addGeometry(self, name, geometryType, argumentDict):
		self.geometryList.append(Geometry(name, geometryType, argumentDict))

	def addRefinementSurface(self, name, lmin, lmax, nrLayers, extraArgumentDict=None):
		self.refinementSurfaceList.append(RefinementSurface(name, lmin, lmax, nrLayers, extraArgumentDict=extraArgumentDict))

	def addFeature(self, name, level):
		self.featureList.append(Feature(name, level))

	def addRefinementRegion(self, name, mode, levels):
		self.refinementRegionList.append(RefinementRegion(name, mode, levels))

	def write(self, folder, ending=''):
		filePath = folder + 'snappyHexMeshDict' + ending

		f = open(filePath, 'w')

		# Write FoamFile header
		f.write(self.header)
		self.writeDict(f, self.FoamFile, 'FoamFile')
		f.write('\n')
		self.writeDict(f, self.actionControl, 'actionControl', wrapper=False)
		f.write('\n')

		# Write geometry
		f.write('geometry\n{\n')

		for i in range(len(self.geometryList)):
			g = self.geometryList[i]
			f.write('\t{0}\n'.format(g.name))
			f.write('\t{\n')
			f.write('\t\ttype {0};\n'.format(g.geometryType))
			for k, v in g.argumentDict.items():
				f.write('\t\t{0} {1};\n'.format(k, v))

			f.write('\t}\n\n')

		f.write('}\n\n')

		# Write castellatedMeshControls
		f.write('castellatedMeshControls\n{\n')

		for k, v in self.castellatedMeshControls.items():
			f.write('\t{0} {1};\n'.format(k, v))
		f.write('\n')

		# Write features
		f.write('\tfeatures\n\t(\n')
		for i in range(len(self.featureList)):
			f.write('\t\t{\n')
			r = self.featureList[i]
			f.write('\t\t\tfile \"{0}\";\n'.format(r.fileName))
			f.write('\t\t\tlevel {:.0f};\n'.format(r.level))
			f.write('\t\t}')
		f.write('\n\t);\n\n')

		# Write refinementSurfaces
		f.write('\trefinementSurfaces\n\t{\n')
		for i in range(len(self.refinementSurfaceList)):
			r = self.refinementSurfaceList[i]
			f.write('\t\t{0}\n'.format(r.name))
			f.write('\t\t{\n')
			f.write('\t\t\tlevel ({:.0f} {:.0f});\n'.format(r.lmin, r.lmax))
			if r.extraArgumentDict != None:
				for k, v in r.extraArgumentDict.items():
					f.write('\t\t\t{0} {1};\n'.format(k, v))
			f.write('\t\t}\n\n')

		f.write('\t}\n\n')

		# Write refinementRegions
		f.write('\trefinementRegions\n\t{\n')
		for i in range(len(self.refinementRegionList)):
			r = self.refinementRegionList[i]
			f.write('\t\t{0}\n'.format(r.name))
			f.write('\t\t{\n')
			f.write('\t\t\tmode {0};\n'.format(r.mode))
			f.write('\t\t\tlevels (')
			if len(r.levels.shape) > 1:
				for j in range(len(r.levels)):
					f.write('({:.6f} {:.0f})'.format(r.levels[j, 0], r.levels[j, 1]))
					if j != len(r.levels) - 1:
						f.write(' ')
			else:
				f.write('({:.6f} {:.0f})'.format(r.levels[0], r.levels[1]))

			f.write(');\n')
			f.write('\t\t}\n\n')

		f.write('\t}\n')

		f.write('}\n\n')

		# Write snap controls
		self.writeDict(f, self.snapControls, 'snapControls')

		# Write layer controls
		f.write('addLayersControls\n{\n')

		# Mesh shrink solver settings
		f.write('\tmeshShrinker displacementMotionSolver;\n')
		f.write('\tsolver       displacementLaplacian;\n')
		f.write('\tdisplacementLaplacianCoeffs\n\t{\n')
		f.write('\t\tdiffusivity    quadratic inverseDistance 1(wall);\n\t}\n\n')

		f.write('\tlayers\n\t{\n')
		for i in range(len(self.refinementSurfaceList)):
			r = self.refinementSurfaceList[i]
			f.write('\t\t{0}\n'.format(r.name))
			f.write('\t\t{\n')
			f.write('\t\t\tnSurfaceLayers {:.0f};\n'.format(r.nrLayers))
			f.write('\t\t}\n')

		f.write('\t}\n\n')

		for k, v in self.addLayersControls.items():
			f.write('\t{0} {1};\n'.format(k, v))

		f.write('}\n\n')

		# Write meshQualityControls
		self.writeDict(f, self.meshQualityControls, 'meshQualityControls')

		# Write final options
		for k, v in self.finalOptions.items():
			f.write('{0} {1};\n'.format(k, v))

		f.close()

	def writeSurfaceFeatureExtractDict(self, folder, geometryFile):
		# --------- Write surfaceFeatureExtract ------------------------------------------
		currentFolder = os.path.abspath(os.path.dirname(__file__))

		f = open(currentFolder+'/surfaceFeatureExtractSettings', 'r')
		settgingsLines = f.readlines()
		f.close()

		filePath = folder + 'surfaceFeatureExtractDict'

		f = open(filePath, 'w')

		# Write FoamFile header
		f.write(self.header)
		self.writeDict(f, self.FoamFileFeature, 'FoamFile')

		if isinstance(geometryFile, list):
			for i in range(len(geometryFile)):
				f.write('{0}\n'.format(geometryFile[i]))
				f.write('{\n')
				for line in settgingsLines:
					f.write(line)
				f.write('}\n\n')
		else:
			f.write('{0}\n'.format(geometryFile))
			f.write('{\n')
			for line in settgingsLines:
				f.write(line)
			f.write('}\n\n')


		f.close()

