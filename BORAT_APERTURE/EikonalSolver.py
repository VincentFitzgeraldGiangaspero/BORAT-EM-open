"""
27/03/23 @Vincent

Eikonal Integrator

"""

import numpy as np
import pyvista as pv
from scipy import integrate
import math

from scipy.spatial.transform import Rotation as Rot
from scipy.interpolate import griddata

from scipy.interpolate import LinearNDInterpolator
# from scipy.spatial import GlobalData.kdtree


# from BORAT3D_global import pvMesh,GlobalData.MeshPoints,GlobalData.RefractiveIndexArray,GlobalData.kdtree,GlobalData.BoundaryPEC,GlobalData.BoundarySurfaces,GlobalData.BoundarySurfacesCFD
import GlobalData


def ComputeDistance(xyz):

    # get indices of closest three points to use as vetice to calculate distance to
    d, idx_of_point_in_mesh = GlobalData.kdtreeBoundaryCFD.query(xyz, 3)

    # anchor point to span a plane from
    anchor_points = GlobalData.MeshPointsBoundaryCFD[idx_of_point_in_mesh[0], :]

    # use next two nearest points to span a plane with two vectors
    # from anchor point
    plane_points = GlobalData.MeshPointsBoundaryCFD[idx_of_point_in_mesh[1:], :]
    plane_vecs = np.array(plane_points)
    plane_vecs[0, :] -= anchor_points
    plane_vecs[1, :] -= anchor_points

    # calculate normal vectors of the planes
    normals = np.cross(plane_vecs[0, :], plane_vecs[1, :])

    # distance from each point to its anchor point for spanning a plane
    PQ = anchor_points - xyz
    # distance is dot product between normal and PQ vector
    # since normals and PQ are arrays of vectors
    # use einsum to calc dot product along first dimension
    # dists = np.einsum('ij,ij->i', PQ, normals)
    dists = np.dot(PQ, normals)

    return dists


###################################################################
"""
Event Functions
"""
# region
###################################################################


def BoundaryEvent(t, y):

    point = pv.PointSet(np.array(y[0:3], dtype=float))
    dist = point.compute_implicit_distance(GlobalData.BoundarySurfaces)
    flag = dist['implicit_distance'][0]

    return flag


BoundaryEvent.direction = 1
BoundaryEvent.terminal = True

###################################################################


def ReflectionEvent(t, y):

    # method 1 for 3D, faster
    point = pv.PointSet(np.array(y[0:3], dtype=float))
    dist = point.compute_implicit_distance(GlobalData.BoundaryPEC)
    flag = dist['implicit_distance'][0]

    return flag


ReflectionEvent.direction = 0
ReflectionEvent.terminal = True

###################################################################


def CutoffEvent(t, y):

    flag = 1

    n, _, _, _ = RefractiveIndexInterpolation([y[0], y[1], y[2]])

    if n < _CutoffValue:
        flag = 0

    return flag


CutoffEvent.terminal = True
CutoffEvent.direction = 0

# endregion
###################################################################
"""
Reflection Functions
"""
# region
###################################################################


def ReflectionFresnel(xyz):

    # here you flip???? CHECK IT
    # GlobalData.BoundaryPEC.compute_normals(flip_normals=True)
    idCell = GlobalData.BoundaryPEC.find_closest_cell(xyz[0:3])
    idPoint = GlobalData.BoundaryPEC.find_closest_point(xyz[0:3])

    n = GlobalData.BoundaryPEC.point_normals[idPoint]  # this looks better

    # wave direction
    ki = [xyz[3], xyz[4], xyz[5]]

    ki = ki/np.linalg.norm(ki)

    """
    CHECK THIS EXPONENTIAL WHERE DOES IT GO
    """
    kr = ki-2 * (np.dot(ki, n) / np.linalg.norm(n)**2) * n

    kr = kr/np.linalg.norm(kr)  # to normalize

    Ei = np.array([xyz[6], xyz[7], xyz[8]])

    """
    ##### TO BE FIXED
    """
    # PEC
    Rp = 1
    Rn = -1

    Rp = 1
    Rt = -1

    # 2 #diadic components, works the same way

    ei_t = np.cross(n, ki)/np.linalg.norm(np.cross(n, ki))
    er_t = ei_t

    ei_p = np.cross(ki, ei_t)
    er_p = np.cross(kr, er_t)

    Rdiadic = Rp*np.outer(er_p, ei_p) + Rt*np.outer(er_t, ei_t)

    Er2 = np.dot(Rdiadic, Ei)

    # 3) just a formula

    Er3 = -Ei+2*np.cross(np.cross(n, Ei), n)

    # 1) normal and parallel components of polarization
    ei_n = np.cross(ki, n)
    ei_n = ei_n/np.linalg.norm(ei_n)

    ei_p = np.cross(ei_n, ki)

    er_p = -ei_p+2 * (np.dot(n, ei_p)) * n
    er_p = er_p/np.linalg.norm(er_p)

    Erp_x = ei_p[0]*er_p[0] + ei_p[1]*er_p[0] + ei_p[2]*er_p[0]
    Erp_y = ei_p[0]*er_p[1] + ei_p[1]*er_p[1] + ei_p[2]*er_p[1]
    Erp_z = ei_p[0]*er_p[2] + ei_p[1]*er_p[2] + ei_p[2]*er_p[2]

    Ern_x = ei_n[0]*ei_n[0] + ei_n[0]*ei_n[1] + ei_n[0]*ei_n[2]
    Ern_y = ei_n[1]*ei_n[0] + ei_n[1]*ei_n[1] + ei_n[1]*ei_n[2]
    Ern_z = ei_n[2]*ei_n[0] + ei_n[2]*ei_n[1] + ei_n[2]*ei_n[2]

    Ern = np.array([Ern_x, Ern_y, Ern_z])
    Erp = np.array([Erp_x, Erp_y, Erp_z])

    Er1 = (Ern*Rn + Erp*Rp) * Ei

    # Er2 more precise, Er3 doesn't change the sign properly\
    Er = Er2

    if np.linalg.norm(np.cross(n, ki)) < 1e-3:
        # if np.dot(ki, n)>0:
        #     n=-n
        Er = -Ei+2*np.cross(np.cross(n, Ei), n)
        # check this second equation

        ERr_method2 = -Ei+np.dot(2*n, Ei)

    return kr[0], kr[1], kr[2], Er[0], Er[1], Er[2]


# endregion
###################################################################
"""
Eikonal Equation 4 polarization
"""
# region
###################################################################
"""
change to interpolation
"""


def EikonalEquation(t, x_y_z):
    x, y, z, u, v, w, ex, ey, ez, _, _ = x_y_z

    if _Solver.Ionized_solution:

        n, nx, ny, nz = RefractiveIndexInterpolation([x, y, z])

    else:

        n, nx, ny, nz = RefractiveIndexH([x, y, z])

    yp = np.zeros(shape=(11,))

    if n==0:
        n
    
    yp[0] = u/(n**2)
    yp[1] = v/(n**2)
    yp[2] = w/(n**2)

    yp[3] = nx/n
    yp[4] = ny/n
    yp[5] = nz/n

    yp[6] = - (ex*nx/n + ey*ny/n + ez*nz/n) * u / ((u**2 + v**2 + w**2))
    yp[7] = - (ex*nx/n + ey*ny/n + ez*nz/n) * v / ((u**2 + v**2 + w**2))
    yp[8] = - (ex*nx/n + ey*ny/n + ez*nz/n) * w / ((u**2 + v**2 + w**2))

    yp[9] = n
    yp[10] = np.sqrt(yp[0]**2+yp[1]**2+yp[2]**2)

    return yp


# endregion
###################################################################
"""
Eikonal Integration
"""
# region
###################################################################


def IntegrateEikonal(y0):

    solution = []

    t0 = 0
    tf = 10000
    maxStep = _Solver.maxStep
    methodOde = _Solver.method
    Rtol = _Solver.Tolerance

    isInside = True

    if _Reflections and _CutoffFlag:
        eventsIntegration = (BoundaryEvent, ReflectionEvent, CutoffEvent)

    elif not _Reflections and _CutoffFlag:

        eventsIntegration = (BoundaryEvent, CutoffEvent)

    elif _Reflections and not _CutoffFlag:

        eventsIntegration = (BoundaryEvent, ReflectionEvent)

    else:
        eventsIntegration = (BoundaryEvent)

    while isInside:

        y_ivp = integrate.solve_ivp(EikonalEquation, [
            t0, tf], y0, max_step=maxStep, dense_output=True, events=eventsIntegration, method=methodOde, rtol=Rtol)

        solution.append(y_ivp)

        # print(y_ivp.nfev, ' steps required.')

        if y_ivp.status == 1:
            if y_ivp.t_events[0].size != 0:

                # print("Out of Domain")
                isInside = False

            if _CutoffFlag:

                # print("Cutoff")
                isInside = False

            if _Reflections:

                if y_ivp.t_events[1].size != 0:

                    # print("Reflection")

                    t0 = y_ivp.t_events[1][0]

                    yR = y_ivp.y_events[1][0].tolist()

                    v1, v2, v3, er1, er2, er3 = ReflectionFresnel(yR[0:9])

                    # move away from reflection point
                    yR[0] = yR[0] + 0.01*v1
                    yR[1] = yR[1] + 0.01*v2
                    yR[2] = yR[2] + 0.01*v3

                    yR[3] = v1
                    yR[4] = v2
                    yR[5] = v3

                    yR[6] = er1
                    yR[7] = er2
                    yR[8] = er3

                    y0 = yR

        elif y_ivp.status == 0:
            # print("End of steps \n")
            isInside = False

        else:
            # print("Integration Failed \n")
            isInside = False

    return solution

###################################################################


def Solve(Solver, InitialConditions):

    global _nSteps
    global _Reflections
    global _CutoffValue
    global _CutoffFlag
    global _Solver

    _nSteps = Solver.raysSampling
    _Reflections = Solver.Reflections
    _CutoffValue = Solver.CutoffValue
    _CutoffFlag = Solver.CutoffFlag

    _Solver = Solver

    # storage of Eikonal Solution
    RaySolutions = np.zeros(shape=(12, _nSteps))

    solution = IntegrateEikonal(InitialConditions)

    Steps = _nSteps//len(solution)

    for nSolution in range(len(solution)):

        tspan = np.linspace(
            solution[nSolution].t[0], solution[nSolution].t[-1], Steps)

        RaySolutions[:11, nSolution*Steps:nSolution *
                     Steps+Steps] = solution[nSolution].sol(tspan)

    RaySolutions[11, :] = np.sqrt(RaySolutions[3, :]**2 +
                                  RaySolutions[4, :]**2 +
                                  RaySolutions[5, :]**2)

    return RaySolutions


# endregion
###################################################################
"""
Refractive Index Functions
"""
# region
###################################################################


def RefractiveIndexH(xyz):

    n = 1
    nx = 0
    ny = 0
    nz = 0

    return [n, nx, ny, nz]

###################################################################


def RefractiveIndexLuneburg(xyz):
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    R_xyz = np.sqrt((x**2 + y**2 + z**2))
    R = 2
    if R_xyz > R:
        n = 1
        nx = 0
        ny = 0
        nz = 0
    else:
        n = np.sqrt(2-(R_xyz/R)**2)
        nx = (-2*x/(R**2))/(2*np.sqrt(2-(R_xyz/R)**2))
        ny = (-2*y/(R**2))/(2*np.sqrt(2-(R_xyz/R)**2))
        nz = (-2*z/(R**2))/(2*np.sqrt(2-(R_xyz/R)**2))

    return n, nx, ny, nz


###################################################################

def RefractiveIndexHelical(xyz):

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    R_xy = np.sqrt((x**2 + y**2))
    R = 2
    n0 = 2

    n = n0/(1+(R_xy/R)**2)
    nx = (-2*n0*x/(R**2))/((1+(R_xy/R)**2)**2)
    ny = (-2*n0*y/(R**2))/((1+(R_xy/R)**2)**2)
    nz = 0

    return n, nx, ny, nz


###################################################################
"""
def RefractiveIndexInterpolation

"""


def RefractiveIndexInterpolationFULL(xyz):

    if GlobalData.Solver.Ionized_solution:

        point = pv.PointSet(np.array(xyz, dtype=float))
        dist = point.compute_implicit_distance(GlobalData.BoundarySurfacesCFD)
        flag = dist['implicit_distance'][0]

        if flag < 0:

            _, index = GlobalData.kdtree.query(xyz, k=8)

            n = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                               GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[0, index], xyz, method='nearest')[0]
            nx = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                                GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[1, index], xyz, method='nearest')[0]
            ny = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                                GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[2, index], xyz, method='nearest')[0]
            nz = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                                GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[3, index], xyz, method='nearest')[0]

            if n < GlobalData.Solver.CutoffValue:
                n = GlobalData.Solver.CutoffValue-0.0000001

        else:

            n = 1
            nx = 0
            ny = 0
            nz = 0
    else:

        n = 1
        nx = 0
        ny = 0
        nz = 0

    return n, nx, ny, nz


def RefractiveIndexInterpolation(xyz):

    if GlobalData.Solver.Ionized_solution:

        _, index = GlobalData.kdtree.query(xyz, k=8)

        n = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                            GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[0, index], xyz, method='nearest')[0]
        nx = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                            GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[1, index], xyz, method='nearest')[0]
        ny = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                            GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[2, index], xyz, method='nearest')[0]
        nz = griddata(np.c_[GlobalData.MeshPoints[index, 0], GlobalData.MeshPoints[index, 1],
                            GlobalData.MeshPoints[index, 2]], GlobalData.RefractiveIndexArray[3, index], xyz, method='nearest')[0]

        if n < GlobalData.Solver.CutoffValue:
            n = GlobalData.Solver.CutoffValue-0.0000001

    
    else:

        n = 1
        nx = 0
        ny = 0
        nz = 0

    return n, nx, ny, nz

# endregion
###################################################################
