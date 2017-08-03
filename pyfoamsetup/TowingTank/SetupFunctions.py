import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import os
import sys
import copy
import subprocess as sp

import pyfoamsetup.TowingTank.TowingTank as TowingTank
import pyfoamsetup.coreLibrary.eMesh     as eMesh

import polymesh.mesh             as Mesh
import polymesh.hydrostatic      as Hydrostatic
import ShipDatabase.ShipDatabase as ShipDatabase

def WaveResistance(shipName, L, U, rudder=False, freeFloating=False, cellLengthFactor=None, timeStepFactor=None, folderName=''):
	''' Function that sets up simulations for steady state conditions, with different speeds, drift and heel angles '''
	ship = ShipDatabase.Ship(shipName)
	
	nu = 1.14e-6
	g  = 9.81

	Ls = ship.Lpp
	Ts = ship.T
	Bs = ship.B
	As = ship.S

	Fr = U/np.sqrt(L*g)

	scale  = L/Ls
	T      = Ts*scale
	B      = Bs*scale
	A      = As*scale**2
	
	if freeFloating:
		VolumeS = ship.Volume
		Volume  = VolumeS*scale**3
		x_bS    = ship.centroid[0]
		x_b     = x_bS*scale

	shipMesh = ship.geometry()

	# Scale ship mesh to right model scale
	shipMesh.scale(scale, scale, scale)
	shipMesh.translate(0, 0, -T)
	shipMesh.calculateFaceData()

	if rudder:
		rudderMesh = ship.rudderGeometry()

		rudderMesh.scale(scale, scale, scale)
		rudderMesh.translate(0, 0, -T)
		rudderMesh.calculateFaceData()

	# Split ship mesh into hull and deck
	cutHeight = np.array([0, 0, L*0.015])
	hull = Hydrostatic.extractWetSurface(shipMesh, p0=cutHeight)
	deck = Hydrostatic.extractWetSurface(shipMesh, p0=cutHeight, up=np.array([0, 0, -1.0]))

	caseNameList = []

	for i in range(len(U)):
		name = folderName + shipName + '/Fr{:.3f}_U{:.3f}'.format(Fr[i], U[i])
		caseNameList.append(name)
	
		tank = TowingTank.TowingTank(name, L, T, B, A, U[i], freeSurface=True, symmetry=True, rudder=rudder, defaultSettingsType='waveResistance')
		if cellLengthFactor != None:
			tank.changeCellLengths(cellLengthFactor)
		if timeStepFactor != None:
			tank.changeTimeStep(timeStepFactor)

		tank.setTwoPartShip()
		tank.setKelvinWake()
		tank.setViscousWake(L, 0)

		if freeFloating:
			tank.setFreeHeaveAndPitch(Volume, x_b)

		tank.writeCaseFiles()
		tank.writeScripts()

		shipMesh.exportObj(tank.geometryFolder+'ship.obj')
		
		if rudder:
			rudderMesh.exportObj(tank.geometryFolder+'rudder.obj')

		
		hull.exportObj(tank.geometryFolder+'hull.obj')
		deck.exportObj(tank.geometryFolder+'deck.obj')

		currentPath = os.getcwd()
		os.chdir(tank.runFolder)
		sp.call('surfaceFeatureExtract', shell=True)
		os.chdir(currentPath)

		edgeMesh = eMesh.eMesh(tank.geometryFolder+'ship.eMesh')
		edgeMesh.removeEdgesAbove(zMax = 2*tank.freeSurfaceWidth)
		edgeMesh.write(tank.geometryFolder+'ship.eMesh')

	# Write run script
	writeRunScripts(caseNameList, tank)

def DriftSimulation(shipName, L, U, alphaDriftDeg, alphaHeelDeg, freeSurface=True, cellLengthFactor=None, timeStepFactor=None, steadyState=False, folderName=''):
	''' Function that sets up simulations for steady state conditions, with different speeds, drift and heel angles '''
	ship = ShipDatabase.Ship(shipName)
	
	nu = 1.14e-6
	g  = 9.81

	Ls      = ship.Lpp
	Ts      = ship.T
	Bs      = ship.B
	
	As = Ts*Ls

	scale  = L/Ls
	T      = Ts*scale
	B      = Bs*scale
	A      = As*scale**2

	alphaDrift = alphaDriftDeg*np.pi/180
	alphaHeel  = alphaHeelDeg*np.pi/180

	shipMesh = ship.geometry()

	# Scale ship mesh to right model scale
	shipMesh.scale(scale, scale, scale)
	shipMesh.translate(0, 0, -T)
	shipMesh.calculateFaceData()

	if rudder:
		rudderMesh = ship.rudderGeometry()

		rudderMesh.scale(scale, scale, scale)
		rudderMesh.translate(0, 0, -T)
		rudderMesh.calculateFaceData()

	if freeSurface:
		# Split ship mesh into hull and deck
		cutHeight = np.array([0, 0, L*0.015])
		hull = Hydrostatic.extractWetSurface(shipMesh, p0=cutHeight)
		deck = Hydrostatic.extractWetSurface(shipMesh, p0=cutHeight, up=np.array([0, 0, -1.0]))

	caseNameList = []

	for k in range(len(U)):
		for i in range(len(alphaHeel)):
			for j in range(len(alphaDrift)):
				name = folderName + shipName + '/U{:.3f}_drift{:.1f}_heel{:.1f}'.format(U[k], alphaDriftDeg[j], alphaHeelDeg[i])
				caseNameList.append(name)
			
				if alphaDrift[j] == 0 and alphaHeel[i] == 0:
					symmetry=True
				else:
					symmetry = False

				tank = TowingTank.TowingTank(name, L, T, B, A, U[k], freeSurface=freeSurface, symmetry=symmetry, rudder=rudder, defaultSettingsType='drift')

				if cellLengthFactor != None:
					tank.changeCellLengths(cellLengthFactor)
				if timeStepFactor != None:
					tank.changeTimeStep(timeStepFactor)

				if tank.freeSurface:
					tank.setTwoPartShip()
					tank.setKelvinWake()
				else:
					if steadyState:
						tank.setSolver('simpleFoam')
						tank.endTime = 4000

					y0 = -tank.L*np.sin(alphaDrift[j]) 
					x0 =  tank.L*np.cos(alphaDrift[j])

					tank.setViscousWake(x0, y0)

				tank.writeCaseFiles()
				tank.writeScripts()

				# Rotate geoemtry to right drift angle and export to case folder
				shipMesh_c = Mesh.copyMesh(shipMesh)
				if rudder:
					rudderMesh_c = Mesh.copyMesh(rudderMesh)

				if alphaHeel[i] != 0:
					shipMesh_c.rotate(-alphaHeel[i], 0, 0)

					if rudder:
						rudderMesh_c.rotate(-alphaHeel[i], 0, 0)

				if alphaDrift[j] != 0:
					shipMesh_c.rotate(0, 0, -alphaDrift[j])

					if rudder:
						rudderMesh_c.rotate(0, 0, -alphaDrift[j])

				shipMesh_c.exportObj(tank.geometryFolder+'ship.obj')
				
				if rudder:
					rudderMesh_c.exportObj(tank.geometryFolder+'rudder.obj')

				if freeSurface:
					hull_c = Mesh.copyMesh(hull)
					deck_c = Mesh.copyMesh(deck)

					if alphaHeel[i] != 0:
						hull_c.rotate(-alphaHeel[i], 0, 0)
						deck_c.rotate(-alphaHeel[i], 0, 0)
					if alphaDrift[j] != 0:
						hull_c.rotate(0, 0, -alphaDrift[j])
						deck_c.rotate(0, 0, -alphaDrift[j])

					hull_c.exportObj(tank.geometryFolder+'hull.obj')
					deck_c.exportObj(tank.geometryFolder+'deck.obj')

					currentPath = os.getcwd()
					os.chdir(tank.runFolder)
					sp.call('surfaceFeatureExtract', shell=True)
					os.chdir(currentPath)

					edgeMesh = eMesh.eMesh(tank.geometryFolder+'ship.eMesh')
					edgeMesh.removeEdgesAbove(zMax = 2*tank.freeSurfaceWidth)
					edgeMesh.write(tank.geometryFolder+'ship.eMesh')

	# Write run script
	writeRunScripts(caseNameList, tank)

def SpeedBoatSimulation(shipName, L, Us, weight=None, freeFloating=True, cellLengthFactor=None, timeStepFactor=None, folderName='', centreOfRotation=None, ):
	''' Function that sets up simulations for steady state conditions, for a speed boat simulation'''
	ship     = ShipDatabase.Ship(shipName)
	shipMesh = ship.geometry()
	
	# Environmental data
	nu    = 1.14e-6
	rho_s = 1025
	rho   = 999.1
	g     = 9.81

	Ls = ship.Lpp
	Bs = ship.B

	Fr = Us/np.sqrt(Ls*g)
	U  = Fr*np.sqrt(L*g)
	Re = U*L/nu

	z_bS    = ship.centerOfGravity[2]
	
	if weight == None:
		As      = ship.S
		Ts      = ship.T
		VolumeS = ship.Volume
		x_bS    = ship.centerOfGravity[0]
		
		shipMesh.translate(0, 0, -Ts)
	else:
		def solveFunc(T):
			targetVolume = weight/rho_s
    
			volume = shipVolume(T, shipName)
    
			return (volume - targetVolume)/targetVolume

		Ts = optimize.newton(solveFunc, ship.T, tol=1.0e-06, maxiter=100)
		print(Ts)

		shipMesh.translate(0, 0, -Ts)
		wetMesh   = Hydrostatic.extractWetSurface(shipMesh)
		VolumeS   = Hydrostatic.calculateVolume(wetMesh)
		As        = Hydrostatic.calculateSurface(wetMesh)
		centroidS = Hydrostatic.calculateVolumeCentroid(wetMesh)
		x_bS      = centroidS[0]

	weight_s = VolumeS*rho_s
	
	scale  = L/Ls
	T      = Ts*scale
	B      = Bs*scale
	A      = As*scale**2
	Volume = VolumeS*scale**3
	x_b    = x_bS*scale
	z_b    = z_bS*scale
	
	# Scale ship mesh to right model scale
	shipMesh.scale(scale, scale, scale)
	shipMesh.calculateFaceData()


	# Split ship mesh into hull and deck
	cutHeight = np.array([0, 0, L*0.015])
	hull = Hydrostatic.extractWetSurface(shipMesh, p0=cutHeight)
	deck = Hydrostatic.extractWetSurface(shipMesh, p0=cutHeight, up=np.array([0, 0, -1.0]))

	caseNameList = []

	for i in range(len(U)):
		name = folderName + shipName + '/Us{:.3f}_Weight{:.1f}_Scale{:.3f}'.format(Us[i]/0.5144444444444, weight_s/1e3, scale)

		caseNameList.append(name)

		tank = TowingTank.TowingTank(name, L, T, B, A, U[i], freeSurface=True, symmetry=True, defaultSettingsType='waveResistance')

		if cellLengthFactor != None:
			print('cellLengthFactor', cellLengthFactor)
			tank.changeCellLengths(cellLengthFactor)
		if timeStepFactor != None:
			tank.changeTimeStep(timeStepFactor)

		tank.setTwoPartShip()
		tank.setKelvinWake()
		tank.endTime /= 2

		if freeFloating:
			tank.setFreeHeaveAndPitch(Volume, x_b, z_b=z_b)
			if type(centreOfRotation) == np.ndarray:
				tank.centreOfRotation = centreOfRotation
			else:
				tank.centreOfRotation = np.array([L, 0, -T])

		tank.writeCaseFiles()
		tank.writeScripts()

		shipMesh.exportObj(tank.geometryFolder+'ship.obj')

		hull.exportObj(tank.geometryFolder+'hull.obj')
		deck.exportObj(tank.geometryFolder+'deck.obj')

		currentPath = os.getcwd()
		os.chdir(tank.runFolder)
		sp.call('surfaceFeatureExtract', shell=True)
		os.chdir(currentPath)

		edgeMesh = eMesh.eMesh(tank.geometryFolder+'ship.eMesh')
		edgeMesh.removeEdgesAbove(zMax = 2*tank.freeSurfaceWidth)
		edgeMesh.write(tank.geometryFolder+'ship.eMesh')

	# Write run script
	writeRunScripts(caseNameList, tank)

def writeRunScripts(caseNameList, tank, folderName=''):
	# Write run script
	if tank.superComputer:
		foamPath = os.environ['FOAM_RUN']

		filePath = foamPath + '/TowingTank/' + folderName + 'viljeRun.sh'
		foalderNameLength = len(folderName)
		f = open(filePath, 'w')

		f.write('#!/bin/bash\n')
		f.write('#PBS -A nn4040k\n')
		f.write('#PBS -N driftAngles\n')
		f.write('#PBS -l walltime={:.0f}:00:00\n'.format(24*len(caseNameList)))
		f.write('#PBS -l select={:.0f}:ncpus=32:mpiprocs=16\n\n'.format(tank.nCPUs/16))
		f.write('module load gcc/4.9.1\n')
		f.write('module load mpt/2.13\n')
		f.write('module load openfoam/3.0+\n\n')
		f.write('cd $PBS_O_WORKDIR\n\n')

		for i in range(len(caseNameList)):
			f.write('cd {0}\n'.format(caseNameList[i][foalderNameLength:]))
			f.write('bash runSim.sh\n')
			f.write('cd $PBS_O_WORKDIR\n\n')
	else:
		if sys.platform == 'darwin':
			f = open(tank.foamPath + '/run.sh', 'w')
		else:
			f = open('run.sh', 'w')

		f.write('#!/bin/bash\n\n')
		if sys.platform == 'darwin':
			f.write('cd $HOME/workingDir/OpenFOAM/run/TowingTank\n\n')
		else:
			f.write('cd $FOAM_RUN/TowingTank\n\n')

		for i in range(len(caseNameList)):

			f.write('cd {0}\n'.format(caseNameList[i]))
			f.write('bash mesh.sh\n')
			f.write('bash runSim.sh\n')
			if sys.platform == 'darwin':
				f.write('cd $HOME/workingDir/OpenFOAM/run/TowingTank\n\n')
			else:
				f.write('cd $FOAM_RUN/TowingTank\n\n')

	f.close()

	if tank.superComputer:
		f = open('localMeshing.sh', 'w')
		f.write('#!/bin/bash\n\n')
		f.write('cd $FOAM_RUN/TowingTank\n\n')

		for i in range(len(caseNameList)):
			f.write('cd {0}\n'.format(caseNameList[i]))
			f.write('bash mesh.sh\n')
			f.write('cd $FOAM_RUN/TowingTank\n\n')

		f.close()

def PostProcessCaseFolder(caseFolder):
	# Find solver 
	simInfo         = FileHandling.readSimulationInfo(caseFolder)
	t, F, F_p, F_v  = FileHandling.readForces(caseFolder)
	t2, M, M_p, M_v = FileHandling.readMoments(caseFolder)
	U, A, L, rho    = FileHandling.findDimensions(caseFolder)

	CT = F[:, 0]/(0.5*rho*A*U**2)
	Cp = F_p[:, 0]/(0.5*rho*A*U**2)
	CF = F_v[:, 0]/(0.5*rho*A*U**2)

	CL = F[:, 1]/(0.5*rho*A*U**2)
	CM = M[:, 2]/(0.5*rho*A*L*U**2)

	if simInfo['solver'] == 'simpleFoam':
		return CT[-1], Cp[-1], CF[-1], CL[-1], CM[-1]

	elif simInfo['solver'] == 'interFoam':
		def forceModelFunction(t, A0, A1, omega, epsilon, k):
			return A0 + (np.exp(-k*t))*A1*np.sin(omega*t + epsilon)

		def findPeriod(t, data):
			meanData = np.trapz(data, x=t)/(t[-1] - t[0])

			data = data - meanData

			n = 10000

			t_r = np.linspace(t[0], t[-1], n)
			dt = t_r[1] - t_r[0]

			spl = interpolate.splrep(t, data)
			data_r = interpolate.splev(t_r, spl)

			ps = np.abs(np.fft.rfft(data_r))**2

			freqs = np.fft.rfftfreq(data_r.size, dt)

			keepIndices = np.where(freqs<5.0)

			freqs = freqs[keepIndices]
			ps    = ps[keepIndices]

			return 1/freqs[np.argmax(ps)]

		def findMeanValue(t, Force):
			meanIndices = np.where(t>t[-1]*(1 - 1/2))

			T_est = findPeriod(t[meanIndices], Force[meanIndices])
			omega_est = 2*np.pi/T_est
			A1_est = 0.5*(np.max(np.max(Force[meanIndices]) - np.min(Force[meanIndices])))
			A0_est = np.mean(Force[meanIndices])
			popt, pcov = optimize.curve_fit(forceModelFunction, t[meanIndices], Force[meanIndices], p0=[A0_est, A1_est, omega_est, 0, 0])

			Force_mean = np.copy(popt[0])

			return Force_mean

		CT_mean = findMeanValue(t, CT)
		Cp_mean = findMeanValue(t, Cp)
		CF_mean = findMeanValue(t, CF)
		CL_mean = findMeanValue(t, CL)
		CM_mean = findMeanValue(t, CM)

		return CT_mean, Cp_mean, CF_mean, CL_mean, CM_mean

	else:
		meanIndices = np.where(t>t[-1]*(1 - 1/2))

		CT_mean = np.trapz(CT[meanIndices], x=t[meanIndices])/(t[meanIndices][-1] - t[meanIndices][0])
		Cp_mean = np.trapz(Cp[meanIndices], x=t[meanIndices])/(t[meanIndices][-1] - t[meanIndices][0])
		CF_mean = np.trapz(CF[meanIndices], x=t[meanIndices])/(t[meanIndices][-1] - t[meanIndices][0])
		CL_mean = np.trapz(CL[meanIndices], x=t[meanIndices])/(t[meanIndices][-1] - t[meanIndices][0])
		CM_mean = np.trapz(CM[meanIndices], x=t[meanIndices])/(t[meanIndices][-1] - t[meanIndices][0])

		return CT_mean, Cp_mean, CF_mean, CL_mean, CM_mean
	
def shipVolume(T, shipName):
	'''Function that calculates the volume of a ship as a function of the submergence/depth, T'''
	ship     = ShipDatabase.Ship(shipName)
	shipMesh = ship.geometry()

	shipMesh.translate(0, 0, -T)
    
	wetMesh = Hydrostatic.extractWetSurface(shipMesh)
    
	Volume  = Hydrostatic.calculateVolume(wetMesh)
    
	return Volume