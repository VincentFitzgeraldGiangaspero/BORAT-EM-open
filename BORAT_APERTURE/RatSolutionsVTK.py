import CreateVTK
import SolutionPLK
import pyvista as pv
import ExportTecplot

RaySolutionsFile = 'ExoMarsDebug_SPHERE_ionized_80_coarse20k/RayTracing.plk'
outputFolder = 'ExoMarsDebug_SPHERE_ionized_80_coarse20k'
Ray0File = 'ExoMarsDebug_SPHERE_ionized_80_coarse20k/Ray0.vtk'
rho= 0.1

RaySolutions = SolutionPLK.importSolutionList(RaySolutionsFile)[0]
Ray0 = pv.PolyData('ExoMarsDebug_SPHERE_ionized_80_coarse20k/Ray0.vtk')

ExportTecplot.export_ray(
    outputFolder+'/' + 'ExoMarsDebug_SPHERE_NONionized_80_coarse20k.szplt', RaySolutions)


RayCompleteSteps=CreateVTK.CompleteRay(RaySolutions,Ray0,rho)
RayCompleteSteps.save(outputFolder + '/' + 'RaySolutionsAll.vtm')
