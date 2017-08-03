from setuptools import setup

import os

# --------- Additional package files -------------------------------
setupFilePath = os.path.dirname(os.path.realpath(__file__))

def baseFolderPackageFiles(applicationName):
	includePaths = []

	# 0 folder
	searchFolder = setupFilePath + '/pyfoamsetup/' + applicationName + '/BaseFolder/0'

	filesToInclude = os.listdir(searchFolder)

	for file in filesToInclude:
		includePaths.append("BaseFolder/0/"+file)

	# constant
	searchFolder = setupFilePath + '/pyfoamsetup/' + applicationName + '/BaseFolder/constant'

	filesToInclude = os.listdir(searchFolder)

	for file in filesToInclude:
		includePaths.append("BaseFolder/constant/"+file)

	# triSurface
	searchFolder = setupFilePath + '/pyfoamsetup/' + applicationName + '/BaseFolder/constant/triSurface'

	filesToInclude = os.listdir(searchFolder)

	for file in filesToInclude:
		includePaths.append("BaseFolder/constant/triSurface"+file)

	# system
	searchFolder = setupFilePath + '/pyfoamsetup/' + applicationName + '/BaseFolder/system'

	filesToInclude = os.listdir(searchFolder)

	for file in filesToInclude:
		includePaths.append("BaseFolder/system/"+file)

	return includePaths

towingTankPaths          = baseFolderPackageFiles('TowingTank')
foilSimulationPaths      = baseFolderPackageFiles('FoilSimulation')
wingSimulationPaths      = baseFolderPackageFiles('WingSimulation')
propellerSimulationPaths = baseFolderPackageFiles('PropellerSimulation')

# --------- Additional package files for the core library ----------------------------
coreLibraryPaths = ["fileHeader.txt", "LES_Settings", "surfaceFeatureExtractSettings"]

setup(name="pyfoamsetup",
	  version="0.1",
	  description="A library for setting up OpenFoam simulations",
	  author="Jarle A. Kramer",
	  author_email="jarlekramer@gmail.com",
	  license="MIT",
	  packages=["pyfoamsetup.coreLibrary", "pyfoamsetup.TowingTank", "pyfoamsetup.FoilSimulation", "pyfoamsetup.WingSimulation", "pyfoamsetup.PropellerSimulation"],
	  package_data={"pyfoamsetup.TowingTank":towingTankPaths,
	  				"pyfoamsetup.FoilSimulation":foilSimulationPaths,
	  				"pyfoamsetup.WingSimulation":wingSimulationPaths,
	  				"pyfoamsetup.PropellerSimulation":propellerSimulationPaths, 
	  				"pyfoamsetup.coreLibrary":coreLibraryPaths},
	  install_requires=["numpy",],
	  include_package_data=True,
)