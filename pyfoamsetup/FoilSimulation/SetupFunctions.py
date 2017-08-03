import numpy as np

import Geometry.Mesh as Mesh
import FoilDatabase.FoilDatabase as FoilDatabase

import PyFoamSetup.FoilSimulation.FoilSimulation as FoilSimulation

def FoilMotionSimulation(x, y, foilMotion, x_cp=0.4, measurePressure=False, storeTimeSteps=False, fluid='water', meshSetting='medium', timeStepSetting='medium'):
	minIterations = 2000

	mesh = FoilDatabase.create3DMesh(x, y, 5)
	mesh.rotate(0, 0, foilMotion.initialRotation())
	mesh.translate(foilMotion.initialSurge(), foilMotion.initialHeave(), 0)

	x0 = []
	y0 = []

	x0.append(np.max(mesh.verts[:, 0]))
	y0.append(np.min(mesh.verts[:, 1]))

	t, surge, heave, pitch = foilMotion.exportCFDMotion()

	origin = np.array([x_cp, 0, 0])

	name = 'motion/k{:.6f}_Ah{:.6f}_Ap{:.6f}_alpha{:.6f}_epsilonPitch{:.6f}_customHeave{:.6f}__customPitch{:.6f}'.format(foilMotion.k, foilMotion.A_heave, foilMotion.A_pitch*180/np.pi, foilMotion.alpha0*180/np.pi, foilMotion.epsilon_pitch*180/np.pi, np.sum(foilMotion.heave_offset), np.sum(foilMotion.pitch_offset))

	foilSim = FoilSimulation.FoilSimulation(name, foilMotion.U, meshSetting=meshSetting, fluid=fluid, timeStepSetting=timeStepSetting)
	foilSim.addMotion(t, surge, heave, pitch, origin)
	if measurePressure:
		foilSim.samples.samplePressureOnPatch('foilPressure', 'foil', writeInterval = 2)

	foilSim.addViscousWake(x0, y0)

	foilSim.endTime       = foilMotion.endTime
	foilSim.writeInterval = foilMotion.T/72
	if storeTimeSteps:
		foilSim.purgeWrite = 0
	else:
		foilSim.purgeWrite = 2
		
	dt_min                = foilMotion.T/minIterations
	foilSim.deltaT        = min(foilSim.deltaT, dt_min)
	foilSim.maxDeltaT     = min(foilSim.maxDeltaT, dt_min)

	foilSim.writeCaseFiles()
	foilSim.writeScripts()

	mesh.exportObj(foilSim.geometryFolder + 'foil.obj')

	foilMotion.writeMotionData(foilSim.runFolder, 'foilMotionData.txt')
	foilMotion.saveToFile(foilSim.runFolder+'/harmonicMotionParameters.txt')

	return name, foilSim
