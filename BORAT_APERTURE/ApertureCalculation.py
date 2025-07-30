import POIntegral
import SolutionPLK
import pyvista as pv
import POPlots
import numpy as np
import os 

#########################################################################
# region

c = 299792458.0  # Speed of light in m/s
Frequency = 400e6
Rho = 1e-1

outputFolder = 'THESIS_NEW_85'

if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
    
#endregion
#########################################################################
CaseName = 'NEW_85'

FileTubes = 'ExoMars_LARGEsphere_3D_6M_85s_DEBUG_CLOUDS_8k_ION/RayTubeSolution_Exomars.vtk'

tstart_obs=0
tstop_obs= 180
delt=1
pstart_obs=-180
pstop_obs=180
delp=1


FinalRayTube= pv.PolyData(FileTubes)


FinalRayTube.rotate_vector([0, 1, 0], 45,
                           inplace=True, transform_all_input_vectors=True)

Lambda = c/Frequency
K=2*np.pi/Lambda
#########################################################################
""" 
#Aperture Integration
"""
# region

Ethscat, thetadeg, phideg = POIntegral.ApertureIntegrationKim(
    FinalRayTube, 
    tstart_obs, 
    tstop_obs, 
    delt, 
    pstart_obs, 
    pstop_obs, 
    delp, 
    K, 
    Rho)

EthTOT = np.sqrt(abs(Ethscat[0, :, :])**2 + abs(Ethscat[1, :, :])**2)

# endregion
#########################################################################
""" 
#Saving solution
"""
# region


Solutions2Save = [thetadeg, phideg, Ethscat,Frequency]

SBRFile = outputFolder + '/' + CaseName + 'SBR.plk'

print('Saving  SBR solution in PLK file: ',
      SBRFile, ' \n')

SolutionPLK.exportSolutionList(SBRFile, Solutions2Save)

# endregion
#########################################################################
""" 
#PLots
"""
# region

phi4plot = 0

POPlots.SimplePlot(
    thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder)

POPlots.PlotApertureExoMars(
    thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder)


# endregion
#########################################################################
