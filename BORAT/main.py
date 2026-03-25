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
import pytecio
import RefractiveIndex

# Parallelization utilities
from functools import partial
from p_tqdm import p_map
from tqdm.auto import tqdm
from scipy.spatial import KDTree


if __name__ == "__main__":

    ConfigFile = "../config/ExoMars.cfg"

    Solver = ConfigBorat.Read(ConfigFile)

    # Initialize TecIO library with path from config
    if Solver.tecioPath:
        pytecio.initialize_tecio(Solver.tecioPath)

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

    # Load CFD mesh and field data
    if Solver.CFDsolution.lower().endswith(".vtk"):
        pvMesh = pv.read(Solver.CFDsolution)
    else:
        # Assume it's a Tecplot file
        pvMesh = Domain.create(Solver)
        BoundarySurfaces = pvMesh.extract_surface()

    MeshPoints = np.array(pvMesh.points)

    # Build RefractiveIndexArray based on configuration
    if Solver.precomputedRefractive:
        # Load precomputed refractive index and gradients from mesh
        RefractiveIndexArray = np.array([np.array(pvMesh[Solver.RefractiveIndex_VarName])])

        if Solver.precomputedGrad:
            RefractiveIndexArray = np.vstack(
                [
                    RefractiveIndexArray,
                    np.array(pvMesh[Solver.Gradient_X_VarName]),
                    np.array(pvMesh[Solver.Gradient_Y_VarName]),
                    np.array(pvMesh[Solver.Gradient_Z_VarName]),
                ]
            )
    else:
        # Compute refractive index and gradients from CFD data
        RefractiveIndexArray = RefractiveIndex.compute_refractive_index(pvMesh, Solver)

    kdtree = KDTree(MeshPoints)

    # Load boundary geometry
    BoundaryPEC = pv.PolyData(Solver.BoundaryPEC_File).compute_normals(flip_normals=True)
    BoundaryES = pv.PolyData(Solver.BoundaryES_File)
    BoundarySurfacesCFD = pv.PolyData(Solver.BoundarySurfacesCFD_File)

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
    print("Antenna modeling:", totalRay, "Rays,", RayTube0.n_cells, "Tubes")

    if Solver.Parallel:
        print("Parallel Simulation")
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
