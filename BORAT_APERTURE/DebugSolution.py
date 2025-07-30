import ExportTecplot
import SolutionPLK


importFile = 'ExoMars3D_50s_rotated_NEW_N0/RayTracing.plk'
solution = SolutionPLK.importSolutionList(importFile)

RaySolutions = solution[0]

raySolutionFileTecplot = 'ExoMars3D_50s_rotated_NEW_N0/NewRaySolutionTecplot.szplt'
ExportTecplot.export_ray(raySolutionFileTecplot, RaySolutions)
