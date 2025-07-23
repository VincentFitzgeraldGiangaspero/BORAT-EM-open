# BORAT-EM
Python solver for EM propagation in plasma based on geometrical optics ray tracing

## How to run:

0) Import BORAT-EM virtual environment in anaconda for modules dependencies (located in folder env/BORAT-EM_Linux.yaml). In case this doesn't work, create your own virtual env and install the essential modules:

pip install pyvista scipy p-tqdm joblib

1) change Tecio Path in pytecio.py (unfortunately hard coded, to be fixed). If you have troubles with this, deactivate the module and do I/O only trough vtk

## Cases to run: 

1) ExoMars2D_NoIon.cfg -> no ionisation, fast, rays straight lines
2) ExoMars2D.cfg -> ionisation, slower, rays bending
3) ExoMars3D_NoIon.cfg -> no ionisation, fast, rays straight lines
4) ExoMars3D.cfg -> ionisation, slower, rays bending

## Visualization:

1) Tecplot: original CFD file in szplt (CFD folder) + ray solutions in szplt (Output/Case subfolder). Visualize rays through mesh activation
2) Paraview: processed CFD solution generated during run (Output/Case subfolder) + ray solutionns in vtm format.

## Input Files:

1) config folder
2) CFD folder for cfd solutions in vtk and szplt format

## Output Files:

1) Output folder with subfolder for the specific testcase
2) rays in vtm format and szplt

## Manual files conversion: 

1) PARAVIEW -> TECPLOT format VTU
2) TECPLOT -> PARAVIEW format plt 2009



