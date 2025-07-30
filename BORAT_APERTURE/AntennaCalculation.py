import POIntegral4Blackout
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


# outputFolder= './50s'
# outputFolder= './DEBUG'
outputFolder = './73s30aperture'


if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
    
#endregion
#########################################################################
CaseName = '73s_30aperture'


FileTubes = '73s30aperture/RayTubeSolution_Exomars.vtk'


tstart_obs = 0
tstop_obs = 180
delt = 4
pstart_obs = -180
pstop_obs = 180
delp = 4


FinalRayTube= pv.PolyData(FileTubes)


FinalRayTube.rotate_vector([0,1,0], 45,
                           inplace=True, transform_all_input_vectors=True)

# FinalRayTube.save('RotatedTubes.vtk')
# FinalRayTube=FinalRayTube.compute_normals()
# FinalRayTube.cell_data['Normal'] = FinalRayTube.face_normals
    

Lambda = c/Frequency
K=2*np.pi/Lambda
#########################################################################
""" 
#Aperture Integration
"""
# region

Ethscat, thetadeg, phideg = POIntegral4Blackout.ApertureIntegrationKim(
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

SBRFile = outputFolder + '/' + CaseName + 'SBR_phi.plk'

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
    thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder,CaseName)

POPlots.PlotApertureExoMars(
    thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder,CaseName)


# endregion
#########################################################################

