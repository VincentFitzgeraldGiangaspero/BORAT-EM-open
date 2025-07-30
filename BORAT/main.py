"""
BORAT-EM: Main execution script for ray-tracing based EM wave propagation in plasma.
"""

import numpy as np
import timeit
import pyvista as pv
import matplotlib.pyplot as plt

# BORAT internal modules
import CreateVTK
import ConfigBorat
import SourceFunctions
import SolutionPLK
import EikonalSolver
import POPlots
import ExportTecplot
import GlobalData
import POIntegral4Blackout

# Parallelization utilities
from functools import partial
from p_tqdm import p_map
from tqdm.auto import tqdm
from scipy.spatial import KDTree
from multiprocessing import Pool
from itertools import repeat
from joblib import Parallel, delayed


def EikonalParallel(nRay):
    y = EikonalSolver.Solve(Solver, BoundarySurfaces, BoundaryPEC, InitialConditions[nRay])
    print("Integration finished for Ray n =", nRay)
    return y


if __name__ == "__main__":

    ConfigFile = "config/ExoMars.cfg"

    Solver = ConfigBorat.Read(ConfigFile)
    outputFolder = Solver.OutputFolder
    CaseName = Solver.Case

    SolutionPLK.exportSolution(outputFolder + "/SolverConfig.plk", Solver)

    print("\n********************************************")
    print("Start BORAT-EM: v1")
    print("Input File Name:", ConfigFile)
    print("********************************************\n")

    start = timeit.default_timer()

    ##############################
    # Domain setup
    ##############################

    # for tecplot
    # pvMesh = Domain.create(Solver)
    # pvMesh.save(outputFolder + 'CFD.vtk')
    # BoundarySurfaces = pvMesh.extract_surface()

    # Load CFD mesh and field data
    pvMesh = pv.read(Solver.CFDsolution)
    MeshPoints = np.array(pvMesh.points)

    RefractiveIndexArray = np.vstack(
        [
            np.array(pvMesh["RefractiveIndex"]),
            np.array(pvMesh["Gradient X Import"]),
            np.array(pvMesh["Gradient Y Import"]),
            np.array(pvMesh["Gradient Z Import"]),
        ]
    )

    kdtree = KDTree(MeshPoints)

    # Load boundary geometry
    BoundaryPEC = pv.PolyData("3DMeshes/ExoMars_Capsule_polydata.vtk").compute_normals(flip_normals=True)
    BoundaryES = pv.PolyData("3DMeshes/Sphere20_capsule_coarse.vtk")
    BoundarySurfacesCFD = pv.PolyData("3DMeshes/BoundaryCFD.vtk")

    MeshPointsBoundaryCFD = np.array(BoundarySurfacesCFD.points)
    kdtreeBoundaryCFD = KDTree(MeshPointsBoundaryCFD)

    BoundarySurfacesCFD.save(outputFolder + "/BoundarySurfacesCFD.vtk")
    BoundarySurfaces = BoundaryES

    step_domain = timeit.default_timer()
    print("Creation Data Run Time:", step_domain - start, "seconds")

    ##############################
    # Export Domain
    ##############################

    pvMesh.save(outputFolder + "/CFD.vtk")
    BoundaryES.save(outputFolder + "/BoundarySurfaceES.vtk")
    BoundaryPEC.save(outputFolder + "/BoundaryPEC.vtk")
    BoundarySurfaces.save(outputFolder + "/BoundarySurfaceSingle.vtk")

    ##############################
    # Global Variables
    ##############################

    GlobalData.pvMesh = pvMesh
    GlobalData.MeshPoints = MeshPoints
    GlobalData.RefractiveIndexArray = RefractiveIndexArray
    GlobalData.kdtree = kdtree
    GlobalData.BoundaryPEC = BoundaryPEC
    GlobalData.BoundaryES = BoundaryES
    GlobalData.BoundarySurfacesCFD = BoundarySurfacesCFD
    GlobalData.BoundarySurfaces = BoundarySurfaces
    GlobalData.kdtreeBoundaryCFD = kdtreeBoundaryCFD
    GlobalData.MeshPointsBoundaryCFD = MeshPointsBoundaryCFD
    GlobalData.Solver = Solver

    ##############################
    # Create Antenna Rays
    ##############################

    InitialConditions, Ray0, RayTube0 = SourceFunctions.Antenna(Solver)
    Ray0.save(outputFolder + "/Ray0.vtk")
    RayTube0.save(outputFolder + "/RayTube0.vtk")

    ##############################
    # Ray Tracing (Eikonal Solver)
    ##############################

    totalRay = Ray0.n_points
    print("Antenna modeling:", totalRay, "Rays,", RayTube0.n_faces, "Tubes")

    if Solver.Parallel:
        print("Parallel Simulation")
        if Solver.Parallel_Method == "joblib":
            EikonalSolutions = Parallel(n_jobs=Solver.nProcess)(delayed(EikonalParallel)(k) for k in range(totalRay))
        elif Solver.Parallel_Method == "pool":
            with Pool(processes=Solver.nProcess) as pool:
                EikonalSolutions = pool.starmap(
                    EikonalSolver.Solve, zip(repeat(Solver), repeat(BoundaryES), repeat(BoundaryPEC), InitialConditions)
                )
        elif Solver.Parallel_Method == "p_map":
            EikonalSolutions = p_map(partial(EikonalSolver.Solve, Solver), InitialConditions)
    else:
        print("Serial Simulation")
        EikonalSolutions = [EikonalSolver.Solve(Solver, y0) for y0 in tqdm(InitialConditions)]

    ##############################
    # Store Eikonal Solutions
    ##############################

    RaySolutions = np.zeros((12, Solver.raysSampling, totalRay))
    for nRay in range(totalRay):
        RaySolutions[:, :, nRay] = EikonalSolutions[nRay]

    RayFinalStep = CreateVTK.FinalRay(RaySolutions, Ray0, Solver.Rho)
    RayFinalStep.save(outputFolder + "/RaySolutionsFinal_" + CaseName + ".vtk")
    SolutionPLK.exportSolution(outputFolder + "/RayTracing.plk", RaySolutions)
    ExportTecplot.export_ray(Solver.raySolutionFileTecplot, RaySolutions)

    ##############################
    # Create Ray Tubes
    ##############################

    FinalRayTube = CreateVTK.FinalTubes(RaySolutions, RayTube0, Solver.Rho)
    FinalRayTube.save(outputFolder + "/RayTubeSolution_" + CaseName + ".vtk")

    ##############################
    # Aperture Integration
    ##############################

    FinalRayTube.rotate_vector([0, 1, 0], 45, inplace=True, transform_all_input_vectors=True)

    Ethscat, thetadeg, phideg = POIntegral4Blackout.ApertureIntegrationKim(
        FinalRayTube,
        Solver.tstart_obs,
        Solver.tstop_obs,
        Solver.delt,
        Solver.pstart_obs,
        Solver.pstop_obs,
        Solver.delp,
        Solver.K,
        Solver.Rho,
    )

    ##############################
    # Save Final Output
    ##############################

    Solutions2Save = [thetadeg, phideg, Ethscat, Solver.Lambda, Solver.K, Solver.Frequency]
    SolutionPLK.exportSolutionList(outputFolder + "/SBR.plk", Solutions2Save)

    ##############################
    # Plots
    ##############################

    phi4plot = 0
    POPlots.SimplePlot(thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder)
    POPlots.PlotApertureExoMars(thetadeg[phi4plot], Ethscat[:, phi4plot, :], outputFolder)

    stop = timeit.default_timer()
    print("\nTotal Run Time:", stop - start, "seconds")
    print("End of simulation")
