
import CreateVTK
import SolutionPLK

Rho = 1e-2

FileTubes = 'ExoMars_LARGEsphere_140k_NO_ION_REF/RayTubeSolution_Exomars.vtk'


    FinalRayTube = CreateVTK.FinalTubes(RaySolutions, RayTube0, Solver.Rho)
    FinalRayTube.save(outputFolder + '/' +
                      'RayTubeSolution_' + CaseName + '.vtk')