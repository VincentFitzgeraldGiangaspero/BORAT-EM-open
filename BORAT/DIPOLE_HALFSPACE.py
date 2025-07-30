"""
26/03/23 @Vincent

NEW MODULAR SBR RAYTRACING

"""

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os

import SourceFunctions
import DomainFunctions
import EikonalSolver
import SolutionPLK
import Functions
import POIntegral

import POPlots
import CreateVTK

from functools import partial
from p_tqdm import p_map, p_umap, p_imap, p_uimap
from tqdm.auto import tqdm

from joblib import Parallel, delayed
###################################################################

#constants
epsilon0 = 8.85e-12  # permittivity of free space
rad = np.pi/180   # radiants conversion
c = 299792458.0  # Speed of light in m/s
_cutoff = 0.001

###################################################################
"""
Parallel Function
"""
#region


def EikonalParallel(nRay):
    y = EikonalSolver.Solve(
        Reflections, Cutoff, nSteps, InitialConditions[nRay], BoundarySurfaces, BoundaryPEC)
    print("Integration finished for Ray n  = ", nRay)
    return y


#endregion
###################################################################
"""
New Vecchi and Kim ray tube approach

"""

CaseName = 'Exomars_antenna_4thesis'
parallel = True
nProcess = 4
#*************


rho = 1e-2
a = 0.3


freq = 400e6



# [1.253 , 0 , 1.067]
# a_location = [1.253 , 0 , 1.067]
a_location = [0, 0, 0]

rotationAngle = 45
Aperture = 70
rotationAxis = [1, 0, 0]
R_exit = 16   # *_lambda

totalRays = 5000
nDivisionSource = 5
nDivisionES = 40

# observation points
tstart_obs = 0
tstop_obs = 360
delt = 4

pstart_obs = 0
pstop_obs = 0
delp = 1

iPol = 1


Reflections = False
Cutoff = False
nSteps = 10

###################################################################


_lambda = c / freq  # Wavelength in m

k = 2 * np.pi / _lambda  # Wave number
ka = k*a
# R_exit=R_exit*_lambda
R_exit = R_exit

p = 4*np.pi*epsilon0/(k**2)

"""
NEW EO AMPLITUDE
"""

# E0_amplitude = p*(k**2)/(epsilon0*4*np.pi)*1/rho
E0_amplitude = p*(k**2)/(epsilon0*4*np.pi)

outputFolder = CaseName


if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)


print('outputFolder : ', outputFolder)
print('Frequency : ', freq/1e9, ' GHz \n')
print('E0 : ', E0_amplitude, '\n')
print('𝝺 : ', _lambda, '\n')
print('Wave Number : ', k, '\n')
print('Ka HF : ', ka, '\n')
print('R exit : ', R_exit, '\n')
print('Source Location : ', a_location, '\n')

###################################################################
"""
Exit Aperture
"""
#region

# BoundarySurfaces = DomainFunctions.ESsurfaceExomars(12,a_location, nDivisionES)
# _boundarySurface.compute_normals
# BoundarySurfaces = pv.PolyData('3DMeshes/Sphere_capsule_Exomars.vtk')
BoundarySurfaces = pv.PolyData(
        '3DMeshes/Sphere_Exomars.vtk')
BoundarySurfaces.save(outputFolder + '/BoundarySurface.vtk')

#endregion
###################################################################
"""
Definition Sources
"""
#region 

InitialConditions, Ray0, RayTube0 = SourceFunctions.Antenna_rotated(
    rho, a_location, E0_amplitude, totalRays,rotationAngle, rotationAxis,Aperture, iPol)

Ray0.save(outputFolder + '/' +'Ray0_' + CaseName + '.vtk')
RayTube0.save(outputFolder + '/' + 'RayTube0_' + CaseName + '.vtk')


totalRays = Ray0.n_points
# always true, by definition of the tubes
totalTubes = RayTube0.n_faces


#endregion
###################################################################
"""
Definition PEC
"""
#region 

if Reflections:

    BoundaryPEC = DomainFunctions.ESsurfaceSphere(R_exit, nDivisionES)

    BoundaryPEC.save(outputFolder + '/BoundaryPEC.vtk')
else:
    #dummy filler
    BoundaryPEC = BoundarySurfaces
    
#endregion
###################################################################
"""
Eikonal Integration
"""
#region

print('Antenna modeling: ', Ray0.n_points, ' total Rays - ', RayTube0.n_faces , ' total Tubes \n')

if parallel:
    
    EikonalSolutions = p_map(partial(EikonalSolver.Solve, Reflections, Cutoff, nSteps ,BoundarySurfaces, BoundaryPEC),
                             InitialConditions)
    
else:
    
    EikonalSolutions = []
    counter = 0

    for y0 in InitialConditions:
        print('* Ray n: ', counter+1)
        counter = counter+1

        EikonalSolutions.append(EikonalSolver.Solve(Reflections, Cutoff, 
                                                    nSteps,BoundarySurfaces, BoundaryPEC,y0))
        
    print('\n')


RaySolutions = np.zeros(shape=(12, nSteps, totalRays))

for nRay in range(totalRays):
    RaySolutions[:, :,nRay] = EikonalSolutions[nRay]
    

print('Saving  Ray Tracing solution in PLK file: ',
      outputFolder + '/RayTracing.plk', ' \n')

SolutionPLK.exportSolutionList(
    outputFolder+'/RayTracing.plk', RaySolutions)

#endregion
###################################################################
"""
Saving Eikonal Solutions
"""
#region

# Save Just final points
RayFinalStep=CreateVTK.FinalRay(RaySolutions,Ray0,rho)
RayFinalStep.save(outputFolder + '/' + 'RaySolutionsFinal_' + CaseName + '.vtk')

# Save all points
# RayCompleteSteps=CreateVTK.CompleteRay(RaySolutions,Ray0,rho)
# RayCompleteSteps.save(outputFolder + '/' + 'RaySolutionsAll_' + CaseName + '.vtm')


#endregion
###################################################################
"""
Tubes Calculation
"""
#region 

FinalRayTube=CreateVTK.FinalTubes(RaySolutions,RayTube0,rho)
FinalRayTube.save(outputFolder + '/' + 'RayTubeSolution_' + CaseName + '.vtk')

#endregion
###################################################################
""" 
#Aperture Integration
"""
#region


Ethscat, thetadeg, phideg = POIntegral.ApertureIntegrationKim(
    FinalRayTube, tstart_obs, tstop_obs, delt, pstart_obs, pstop_obs, delp, k,rho)


# step1 = timeit.default_timer()

# Ethscat, thetadeg, phideg = POIntegral_Parallel.ApertureParallel1(
#     FinalRayTube, tstart_obs, tstop_obs, delt, pstart_obs, pstop_obs, delp, k,rho)
# step2 = timeit.default_timer()
# print('Parallel: ', step2 - step1, 'seconds')


Solutions2Save = [thetadeg, phideg, Ethscat, _lambda, k, freq]


print('Saving  SBR solution in PLK file: ',
      outputFolder + '/SBR.plk', ' \n')

SolutionPLK.exportSolutionList(outputFolder+'/SBR.plk', Solutions2Save)

#endregion
###################################################################
""" 
#PLots
"""
phi4plot=0

POPlots.SimplePlot(
    thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder)

POPlots.PlotApertureExoMars(
    thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder)
    
###################################################################
""" 
#Saving solution
"""
#region



#endregion
###################################################################
exit()
