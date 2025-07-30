"""
-- Version 1.F March 2023
"""
#########################################################################
"""
Import
"""
# region


#########################################################################


#########################################################################


# endregion
#########################################################################
# region
"""
Global Variables
"""


# endregion
#########################################################################
"""
Parallel Function
"""
# region




import matplotlib.pyplot as plt
import CreateVTK
import numpy as np
import timeit
import tecplot as tp
import pyvista as pv
import ConfigBorat
import Domain
import SourceFunctions
import SolutionPLK
import EikonalSolver
import Functions
import POIntegral
import POPlots
import DomainFunctions
import ExportTecplot
import GlobalData
from functools import partial
from p_tqdm import p_map, p_umap, p_imap, p_uimap
from tqdm.auto import tqdm
from scipy.spatial import KDTree
from multiprocessing import Pool
from itertools import repeat
from joblib import Parallel, delayed
def EikonalParallel(nRay):
    y = EikonalSolver.Solve(
        Solver, BoundarySurfaces, BoundaryPEC, InitialConditions[nRay])
    print("Integration finished for Ray n  = ", nRay)
    return y


# endregion
#########################################################################


if __name__ == '__main__':

    ConfigFile = 'ExoMarsDebug.cfg'

    #########################################################################
    # Start

    Solver = ConfigBorat.Read(ConfigFile)

    outputFolder = Solver.OutputFolder
    CaseName = Solver.Case

    print("\n")
    print("********************************************\n")
    print("Start BORAT3DP0 MAR 2023: v1.F")
    print("\n")
    print('Input FIle Name : ', ConfigFile)
    print("\n")
    print("********************************************\n")

    start = timeit.default_timer()

    #########################################################################
    """
    Create Domain
    """
    # region

    pvMesh = Domain.create(Solver)

    # BoundarySurfaces = pvMesh.extract_surface()
    MeshPoints = np.array(pvMesh.points)

    RefractiveIndexArray = np.zeros(shape=(4, pvMesh.n_points))
    RefractiveIndexArray[0, :] = np.array(pvMesh['RefractiveIndex'])
    RefractiveIndexArray[1, :] = np.array(pvMesh['gradientX'])
    RefractiveIndexArray[2, :] = np.array(pvMesh['gradientY'])
    RefractiveIndexArray[3, :] = np.array(pvMesh['gradientZ'])

    kdtree = KDTree(
        np.c_[MeshPoints[:, 0], MeshPoints[:, 1], MeshPoints[:, 2]])

    BoundaryPEC = pv.PolyData('3DMeshes/ExoMars_Capsule_polydata.vtk')
    BoundaryPEC = BoundaryPEC.compute_normals(flip_normals=True)

    # BoundarySurfaces = DomainFunctions.ESsurfaceExomars(5, BoundaryPEC.center, 50)
    # BoundarySurfaces = pv.PolyData('3DMeshes/CylinderES.vtk')
    BoundarySurfaces = pv.PolyData(
        '3DMeshes/Sphere_capsule_Exomars.vtk')

    BoundarySurfacesCFD = pvMesh.extract_surface()

    MeshPointsBoundaryCFD = np.array(BoundarySurfacesCFD.points)

    kdtreeBoundaryCFD = KDTree(
        np.c_[MeshPointsBoundaryCFD[:, 0], MeshPointsBoundaryCFD[:, 1], MeshPointsBoundaryCFD[:, 2]])

    BoundarySurfacesCFD.save(outputFolder + '/BoundarySurfacesCFD.vtk')

    step_domain = timeit.default_timer()
    print('Creation Data Run Time: ', step_domain - start, 'seconds')
    print("********************************************\n")

    # endregion
    #########################################################################
    """
    Export Domain
    """
    # region

    # pvMesh.save(outputFolder + '/CFD.vtk')
    BoundarySurfaces.save(outputFolder + '/BoundarySurface.vtk')
    BoundaryPEC.save(outputFolder + '/BoundaryPEC.vtk')

    # endregion
    #########################################################################
    """
    Global Variables
    """
    # region

    GlobalData.pvMesh = pvMesh
    GlobalData.MeshPoints = MeshPoints
    GlobalData.RefractiveIndexArray = RefractiveIndexArray
    GlobalData.kdtree = kdtree
    GlobalData.BoundaryPEC = BoundaryPEC
    GlobalData.BoundarySurfaces = BoundarySurfaces
    GlobalData.BoundarySurfacesCFD = BoundarySurfacesCFD
    GlobalData.kdtreeBoundaryCFD = kdtreeBoundaryCFD
    GlobalData.MeshPointsBoundaryCFD = MeshPointsBoundaryCFD
    GlobalData.Solver = Solver
    # endregion
    #########################################################################
    """
    Create Antenna
    """
    # region

    InitialConditions, Ray0, RayTube0 = SourceFunctions.Antenna(Solver)

    Ray0.save(outputFolder + '/' + 'Ray0.vtk')
    RayTube0.save(outputFolder + '/' + 'RayTube0.vtk')

    # endregion
    #########################################################################
    """
    Eikonal Integration
    """
    # region

    totalRay = Ray0.n_points

    print('Starting Ray Tracing simulation  \n')
    print('Antenna modeling: ', Ray0.n_points,
          ' total Rays - ', RayTube0.n_faces, ' total Tubes \n')

    if Solver.Parallel:

        print('Parallel SImulation')

        if Solver.Parallel_Method == "joblib":

            joblib_start = timeit.default_timer()

            EikonalSolutions = Parallel(n_jobs=Solver.nProcess)(delayed(EikonalParallel)(k)
                                                                for k in range(0, totalRay))
            joblib_stop = timeit.default_timer()

            print('Joblib Run Time: ', joblib_stop - joblib_start, 'seconds')

        elif Solver.Parallel_Method == "pool":

            pool_start = timeit.default_timer()

            with Pool(processes=Solver.nProcess) as pool:

             EikonalSolutions = pool.starmap(
                 EikonalSolver.Solve, zip(repeat(Solver), repeat(BoundarySurfaces), repeat(BoundaryPEC), InitialConditions))

            pool_stop = timeit.default_timer()

            print('Pool Run Time: ', pool_stop - pool_start, 'seconds')

        elif Solver.Parallel_Method == "p_map":

            EikonalSolutions = p_map(partial(EikonalSolver.Solve, Solver),
                                     InitialConditions)
            # EikonalSolutions=p_map(EikonalSolver.Solve,repeat(Solver),
            #                        InitialConditions, repeat(BoundarySurfaces), repeat(BoundaryPEC))

    else:

        print('Serial SImulation')

        EikonalSolutions = []
        counter = 0

        for y0 in tqdm(InitialConditions):

            # for y0 in InitialConditions:

            resultIntegration = EikonalSolver.Solve(
                Solver, y0)
            EikonalSolutions.append(resultIntegration)

        print('\n')

    RaySolutions = np.zeros(shape=(12, Solver.raysSampling, totalRay))

    for nRay in range(totalRay):
        RaySolutions[:, :, nRay] = EikonalSolutions[nRay]

    # endregion
    #########################################################################
    """
    Saving Eikonal Solutions
    """
    # region

    # Save Just final points
    RayFinalStep = CreateVTK.FinalRay(RaySolutions, Ray0, Solver.Rho)
    RayFinalStep.save(outputFolder + '/' +
                      'RaySolutionsFinal_' + CaseName + '.vtk')

    print('Saving  Ray Tracing solution in PLK file: ',
          outputFolder + '/RayTracing.plk', ' \n')

    SolutionPLK.exportSolution(
        outputFolder+'/RayTracing.plk', RaySolutions)

    ExportTecplot.export_ray(Solver.raySolutionFileTecplot, RaySolutions)

    # ExportTecplot.export_ray(Solver.raySolutionFileTecplot, solution_rays)

    # Save all points
    # RayCompleteSteps=CreateVTK.CompleteRay(RaySolutions,Ray0,rho)
    # RayCompleteSteps.save(outputFolder + '/' + 'RaySolutionsAll_' + CaseName + '.vtm')

    # endregion
    #########################################################################
    """
    Tubes Calculation
    """
    # region

    FinalRayTube = CreateVTK.FinalTubes(RaySolutions, RayTube0, Solver.Rho)
    FinalRayTube.save(outputFolder + '/' +
                      'RayTubeSolution_' + CaseName + '.vtk')

    # endregion
    #########################################################################
    """ 
    #Aperture Integration
    """
    # region

    Ethscat, thetadeg, phideg = POIntegral.ApertureIntegrationKim(
        FinalRayTube, Solver.tstart_obs, Solver.tstop_obs, Solver.delt, Solver.pstart_obs, Solver.pstop_obs, Solver.delp, Solver.K, Solver.Rho)

    # endregion
    #########################################################################
    """ 
    #Saving solution
    """
    # region

    Solutions2Save = [thetadeg, phideg, Ethscat,
                      Solver.Lambda, Solver.K, Solver.Frequency]

    print('Saving  SBR solution in PLK file: ',
          outputFolder + '/SBR.plk', ' \n')

    SolutionPLK.exportSolutionList(outputFolder+'/SBR.plk', Solutions2Save)

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

    stop = timeit.default_timer()

    print('\n')
    print('Total Run Time: ', stop - start, ' seconds')
    print('End of simulation ')

    #########################################################################

    exit()
