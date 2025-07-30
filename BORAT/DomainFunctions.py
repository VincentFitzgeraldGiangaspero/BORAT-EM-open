import pyvista as pv
import numpy as np
import meshzoo
import math


###################################################################

def PECsphere(Rin, nDivisionDomainPec):

    points, cells = meshzoo.icosa_sphere(nDivisionDomainPec)

    points = points*Rin
    arrayNodeMap = np.zeros(shape=(cells.shape[0], 4), dtype='<i4')
    arrayNodeMap[:, 0] = 3
    arrayNodeMap[:, 1:4] = cells

    boundaryPEC = pv.PolyData(points, arrayNodeMap.ravel())

    return boundaryPEC

###################################################################


def PECplate(size, N, nDivisionDomainPec):

    L = size

    boundaryPEC = pv.Plane(i_size=L, j_size=L, direction=N,
                           i_resolution=nDivisionDomainPec, j_resolution=nDivisionDomainPec)

    boundaryPEC = boundaryPEC.triangulate()

    return boundaryPEC


###################################################################

def ESsurfaceBOX(L, nDivisionDomainPec):

    boundaryPEC = pv.Box(bounds=(-L, L, -L, L, -L, L))

    boundaryPEC = boundaryPEC.triangulate()

    return boundaryPEC


###################################################################


def ESsurfaceSphere(Rout, nDivisionDomain):

    points, cells = meshzoo.icosa_sphere(nDivisionDomain)

    points = points*Rout
    arrayNodeMap = np.zeros(shape=(cells.shape[0], 4), dtype='<i4')
    arrayNodeMap[:, 0] = 3
    arrayNodeMap[:, 1:4] = cells

    boundarySurface = pv.PolyData(points, arrayNodeMap.ravel())

    return boundarySurface

###################################################################


def ESsurfaceExomars(Rout, center, nDivisionDomain):

    
    
    points, cells = meshzoo.icosa_sphere(nDivisionDomain)

    points = points*Rout
    arrayNodeMap = np.zeros(shape=(cells.shape[0], 4), dtype='<i4')
    arrayNodeMap[:, 0] = 3
    arrayNodeMap[:, 1:4] = cells

    boundarySurface = pv.PolyData(points, arrayNodeMap.ravel())
    boundarySurface.translate(
        (center[0], center[1], center[2]), inplace=True,transform_all_input_vectors=True)
    
    return boundarySurface

###################################################################
