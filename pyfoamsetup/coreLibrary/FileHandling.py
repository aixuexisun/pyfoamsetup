import numpy as np
import os
import shutil
import copy
from collections import OrderedDict

def readForceCoeffs(caseFolder, forceCoeffsName='forceCoeffs'):
	folder = caseFolder + '/postProcessing/' + forceCoeffsName + '/'

	timeFolderList = os.listdir(folder)
	nrFolders = len(timeFolderList)
	timeFolderListNr = np.zeros(nrFolders)
		
	# Sort time folder list
	for i in range(nrFolders):
		try:
			timeFolderListNr[i] = float(timeFolderList[i])
		except:
			timeFolderListNr[i] = 999

	index = np.argsort(timeFolderListNr)

	t  = np.array([])
	Cm = np.array([])
	Cd = np.array([])
	Cl = np.array([])
		
	for i in range(nrFolders):
		timeFolder = timeFolderList[index[i]]
		try:
			f = open(folder+timeFolder+'/'+'forceCoeffs.dat', 'r')
		except FileNotFoundError:
			f = open(folder+timeFolder+'/'+'coefficient.dat', 'r')
		lines = f.readlines()
		f.close()

		nrLines = len(lines)

		# Search for start line
		startLineFound = False
		j = 0
		while not(startLineFound) and j < nrLines:
			if lines[j][0] != '#':
				startLineFound = True
				dataStartLine = copy.deepcopy(j)

			j += 1
	
		nrData = nrLines - dataStartLine

		t  = np.zeros(nrData)
		Cm = np.zeros(nrData)
		Cd = np.zeros(nrData)
		Cl = np.zeros(nrData)

		for j in range(dataStartLine, nrLines):
			data = lines[j].strip().split()

			t[j-dataStartLine]  = float(data[0])
			Cm[j-dataStartLine] = float(data[1])
			Cd[j-dataStartLine] = float(data[2])
			Cl[j-dataStartLine] = float(data[3])

	return (t, Cm, Cd, Cl)

def readForces(caseFolder, forceName='forces'):
	folder = caseFolder + '/postProcessing/'+ forceName + '/'

	timeFolderList = os.listdir(folder)
	nrFolders = len(timeFolderList)
	timeFolderListNr = np.zeros(nrFolders)
		
	# Sort time folder list
	for i in range(nrFolders):
		try:
			timeFolderListNr[i] = float(timeFolderList[i])
		except:
			timeFolderListNr[i] = 999

	index = np.argsort(timeFolderListNr)

	# Use first folder
	timeFolder = timeFolderList[index[0]]

	f = open(folder+timeFolder+'/'+'force.dat', 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)
	startLine = 4
	nrDataPoints = nrLines - startLine

	t   = np.zeros(nrDataPoints)
	F   = np.zeros((nrDataPoints, 3))
	F_p = np.zeros((nrDataPoints, 3))
	F_v = np.zeros((nrDataPoints, 3))

	for i in range(nrDataPoints):
		t[i] = float(lines[i+startLine].strip().split()[0])
		
		line = lines[i+startLine].strip().split('(')

		F_string   = line[1][0:-2]
		F_p_string = line[2][0:-2]
		F_v_string = line[3][0:-2]

		for j in range(3):
			F[i, j] = float(F_string.split()[j])
			F_p[i, j] = float(F_p_string.split()[j])
			F_v[i, j] = float(F_v_string.split()[j])

	return t, F, F_p, F_v

def readMoments(caseFolder, forceName='forces'):
	folder = caseFolder + '/postProcessing/' + forceName + '/'

	timeFolderList = os.listdir(folder)
	nrFolders = len(timeFolderList)
	timeFolderListNr = np.zeros(nrFolders)
		
	# Sort time folder list
	for i in range(nrFolders):
		try:
			timeFolderListNr[i] = float(timeFolderList[i])
		except:
			timeFolderListNr[i] = 999

	index = np.argsort(timeFolderListNr)

	# Use first folder
	timeFolder = timeFolderList[index[0]]

	f = open(folder+timeFolder+'/'+'moment.dat', 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)
	startLine = 4
	nrDataPoints = nrLines - startLine

	t   = np.zeros(nrDataPoints)
	M   = np.zeros((nrDataPoints, 3))
	M_p = np.zeros((nrDataPoints, 3))
	M_v = np.zeros((nrDataPoints, 3))

	for i in range(nrDataPoints):
		t[i] = float(lines[i+startLine].strip().split()[0])
		
		line = lines[i+startLine].strip().split('(')

		M_string   = line[1][0:-2]
		M_p_string = line[2][0:-2]
		M_v_string = line[3][0:-2]

		for j in range(3):
			M[i, j] = float(M_string.split()[j])
			M_p[i, j] = float(M_p_string.split()[j])
			M_v[i, j] = float(M_v_string.split()[j])

	return t, M, M_p, M_v

def changeLine(dictFilePath, lineKeyWord, newLine):
	f = open(dictFilePath, 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)

	f = open(dictFilePath, 'w')

	for i in range(nrLines):
		line = lines[i].strip().split()
		writeNewLine = False

		if len(line) > 0: 
			if line[0] == lineKeyWord:
				f.write(newLine)
				f.write('\n')
				writeNewLine = True
		
		if not(writeNewLine):
			f.write(lines[i])

	f.close()

def deleteLine(dictFilePath, lineKeyWord):
	f = open(dictFilePath, 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)

	f = open(dictFilePath, 'w')

	for i in range(nrLines):
		line = lines[i].strip().split()
		skipLine = False

		if len(line) > 0: 
			if line[0] == lineKeyWord:
				skipLine = True
		
		if not(skipLine):
			f.write(lines[i])

	f.close()

def changeWord(dictFilePath, oldWord, newWord):
	f = open(dictFilePath, 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)
	oldWord_len = len(oldWord)

	f = open(dictFilePath, 'w')

	for i in range(nrLines):
		line = lines[i]

		writeNewLine = False

		if len(line) > 0:
			for j in range(len(line)):
				if line[j:j+oldWord_len] == oldWord:
					writeNewLine = True
					newWordIndex = j
					break

		if writeNewLine:
			j = 0
			while j < len(line):
				if j == newWordIndex:
					f.write(newWord)
					j += oldWord_len
				else:
					f.write(line[j])
					j += 1
		else:
			f.write(line)

	f.close()

def adjustSystemFolder(caseFolder, solver):
	systemFolder = caseFolder + '/system/'
	
	folderList = os.listdir(caseFolder)

	for i in range(len(folderList)):
		folderName = folderList[i]
		folderNameList = folderName.split('.')

		if len(folderNameList) > 1:
			folderPath = caseFolder + '/' + folderName + '/'
			# If right system folder
			if folderNameList[0] == 'system' and folderNameList[1] == solver:
				fileList = os.listdir(folderPath)

				for j in range(len(fileList)):
					shutil.move(folderPath + fileList[j], systemFolder + fileList[j])

				shutil.rmtree(folderPath)
			elif folderNameList[0] == 'system' and not(folderNameList[1] == solver):
				shutil.rmtree(folderPath)

def copyBoundary(boundaryFilePath, boundaryName, newBoundaryName):
	f = open(boundaryFilePath, 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)

	# Search for boundary
	i = 0
	boundaryFound       = False
	boundaryEndingFound = False
	boundaryStartIndex  = 0
	boundaryStopIndex   = 0
	boundaryNrLines     = 0

	while i < nrLines and not(boundaryFound):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == boundaryName:
				boundaryStartIndex = i
				boundaryFound = True

		i += 1

	# Find out how many lines the boundary consists of
	i = boundaryStartIndex
	while i < nrLines and not(boundaryEndingFound):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == '}':
				boundaryStopIndex = i  
				boundaryEndingFound = True

		i += 1

	boundaryNrLines = boundaryStopIndex - boundaryStartIndex

	f = open(boundaryFilePath, 'w')

	for i in range(boundaryStopIndex+1):
		f.write(lines[i])

	f.write('\t'+newBoundaryName+'\n')
	for i in range(boundaryNrLines-1):
		f.write(lines[boundaryStartIndex + 1 + i])

	for i in range(boundaryStopIndex, nrLines):
		f.write(lines[i])

	f.close()

def changeBoundary(boundaryFilePath, boundaryName, newString):
	f = open(boundaryFilePath, 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)

	# Search for boundary
	i = 0
	boundaryFound       = False
	boundaryEndingFound = False
	boundaryStartIndex  = 0
	boundaryStopIndex   = 0
	boundaryNrLines     = 0

	while i < nrLines and not(boundaryFound):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == boundaryName:
				boundaryStartIndex = i
				boundaryFound = True

		i += 1

	# Find out how many lines the boundary consists of
	i = boundaryStartIndex
	while i < nrLines and not(boundaryEndingFound):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == '}':
				boundaryStopIndex = i  
				boundaryEndingFound = True

		i += 1

	boundaryNrLines = boundaryStopIndex - boundaryStartIndex

	f = open(boundaryFilePath, 'w')

	for i in range(boundaryStartIndex+2):
		f.write(lines[i])

	f.write(newString)

	for i in range(boundaryStopIndex, nrLines):
		f.write(lines[i])

	f.close()

def deleteBoundary(boundaryFilePath, boundaryName):
	f = open(boundaryFilePath, 'r')
	lines = f.readlines()
	f.close()

	nrLines = len(lines)

	# Search for boundary
	i = 0
	boundaryFound       = False
	boundaryEndingFound = False
	boundaryStartIndex  = 0
	boundaryStopIndex   = 0
	boundaryNrLines     = 0

	while i < nrLines and not(boundaryFound):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == boundaryName:
				boundaryStartIndex = i
				boundaryFound = True

		i += 1

	# Find out how many lines the boundary consists of
	i = boundaryStartIndex
	while i < nrLines and not(boundaryEndingFound):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == '}':
				boundaryStopIndex = i  
				boundaryEndingFound = True

		i += 1

	boundaryNrLines = boundaryStopIndex - boundaryStartIndex

	f = open(boundaryFilePath, 'w')

	i = 0
	while i < nrLines:
		line = lines[i].strip().split()
		skipLines = False

		if len(line) > 0:
			if line[0] == boundaryName:
				skipLines = True

		if skipLines:
			i += boundaryNrLines+1
		else:
			f.write(lines[i])
			i += 1

def readPatchPressure(caseFolder, sampleName):
	folder = caseFolder + '/postProcessing/' + sampleName + '/'

	timeFolderList = os.listdir(folder)

	nTime = len(timeFolderList)

	f = open(folder + timeFolderList[0] + '/p_walls.raw', 'r')
	lines = f.readlines()
	f.close()

	nCells = len(lines)-2

	x = np.zeros((nTime, nCells))
	y = np.zeros((nTime, nCells))
	z = np.zeros((nTime, nCells))
	p = np.zeros((nTime, nCells))

	for i in range(nTime):
		f = open(folder + timeFolderList[i] + '/p_walls.raw', 'r')
		lines = f.readlines()
		f.close()

		for j in range(2, nCells):
			line = lines[j].strip().split()

			x[i, j-2] = float(line[0])
			y[i, j-2] = float(line[1])
			z[i, j-2] = float(line[2])
			p[i, j-2] = float(line[3])

	return x, y, z, p

def readInletVelocity(caseFolder):
	U = 0

	folder = caseFolder + '/0/include/'

	f = open(folder+'initialConditions', 'r')
	lines = f.readlines()
	f.close()

	n = len(lines)

	for i in range(n):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == 'Umean':
				U = float(line[1].strip(';'))

	return U

def readInitialConditions(caseFolder):
	folder = caseFolder + '/0/include/'
	f = open(folder+'/initialConditions', 'r')
	lines = f.readlines()
	f.close()

	line1 = lines[0].split()
	temp = line1[1]
	xVelocity = float(temp[1:])
	yVelocity = float(line1[2])
	temp = line1[3]
	zVelocity = float(temp[0:len(temp)-2])

	temp = (lines[1].split())[1]
	Umean = float(temp[0:len(temp)-1])
	temp = (lines[2].split())[1]
	turbulentKE = float(temp[0:len(temp)-1])
	temp = (lines[3].split())[1]
	turbulentOmega = float(temp[0:len(temp)-1])
	temp = (lines[4].split())[1]
	turbulentEpsilon = float(temp[0:len(temp)-1])
	temp = (lines[5].split())[1]
	turbulentNut = float(temp[0:len(temp)-1])
	temp = (lines[6].split())[1]
	inputMode = temp[0:len(temp)-1]

	return xVelocity, yVelocity, zVelocity

def readFoilMotion(caseFolderPath, fileName='motion_foil'):
	f = open(caseFolderPath + '/' + fileName, 'r')
	lines=f.readlines()
	f.close()

	nLines = len(lines)
	n = nLines - 2

	t = np.zeros(n)
	heave = np.zeros(n)
	pitch = np.zeros(n)

	for i in range(1, nLines-1):
		line = lines[i].strip().strip(')')[1:].split('(')

		t[i-1] = float(line[0])

		lineData = line[1].split()

		heave[i-1] = float(lineData[1])
		pitch[i-1] = float(lineData[2])

	return t, heave, pitch

def findDimensions(caseFolder):
	filePath = caseFolder + '/system/forces'

	f = open(filePath, 'r')
	lines = f.readlines()
	f.close()

	nLines = len(lines)

	for i in range(nLines):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == 'magUInf':
				U = float(line[1].strip(';'))
			elif line[0] == 'Aref':
  				A = float(line[1].strip(';'))
			elif line[0] == 'lRef':
				L = float(line[1].strip(';'))
			elif line[0] =='rhoInf':
				rho = float(line[1].strip(';'))

	return U, A, L, rho

def readSimulationInfo(caseFolder):
	filePath = caseFolder + '/simulationInfo.txt'

	f = open(filePath, 'r')
	lines = f.readlines()
	f.close()

	simInfo = OrderedDict()

	for i in range(len(lines)):
		line = lines[i].strip().split()

		if len(line) > 0:
			if line[0] == 'smallestSize/L':
				simInfo['smallestSize/L'] = float(line[1])
			elif line[0] == 'solver':
				simInfo['solver'] = line[1]

	return simInfo