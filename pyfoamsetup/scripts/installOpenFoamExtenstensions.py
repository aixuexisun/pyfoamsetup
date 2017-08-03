import os
import sys
import shutil
import subprocess as sp

''' Script that installes the extensions to the OpenFOAM library that is necessary to use some of the functions in the PyFoamSetup library 
Current extensions are as folows:
- twoDimMotion                 - A boundary that allows for arbitrary motion in the xy-plane, using the built in dynamicMeshLibrary of OpenFOAM
- GoldsteinActuationDiskSource - A actuation disk using a volume source to model a Goldstein optimum propeller'''

foamFolder = os.environ['FOAM_RUN']
foamFolder = foamFolder[0:len(foamFolder)-4]
scriptPath    = os.path.dirname(os.path.realpath(__file__))

foamExtensionSourceFolder = scriptPath + '/Extensions'
foamExtensionDestFolder   = foamFolder + '/PyFoamSetupExtensions'

if not os.path.exists(foamExtensionDestFolder):
	os.makedirs(foamExtensionDestFolder)

# Move the source files to a separate folder in the OpenFOAM run folder
destFolder   = foamExtensionDestFolder   +'/twoDimMotion'
sourceFolder = foamExtensionSourceFolder + '/twoDimMotion'
if os.path.exists(destFolder):
	shutil.rmtree(destFolder)

shutil.copytree(sourceFolder, destFolder)
sp.call('cd {}\nwmake'.format(destFolder), shell=True)


destFolder   = foamExtensionDestFolder   +'/goldsteinActuationDiskSource'
sourceFolder = foamExtensionSourceFolder + '/goldsteinActuationDiskSource'
if os.path.exists(destFolder):
	shutil.rmtree(destFolder)

shutil.copytree(sourceFolder, destFolder)
sp.call('cd {}\nwmake'.format(destFolder), shell=True)

