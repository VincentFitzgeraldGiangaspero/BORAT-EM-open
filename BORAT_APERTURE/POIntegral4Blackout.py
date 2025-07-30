"""
26/03/23 @Vincent

Integration of Physical Optics on the Equivalence Surface Boundary at the end of the ray tracing integration

"""

import numpy as np
import math
import itertools

import ShapeFactor

from p_tqdm import p_map, p_umap, p_imap, p_uimap
from tqdm.auto import tqdm

###################################################################

# constants
epsilon0 = 8.85e-12  # permittivity of free space
rad = np.pi/180   # radiants conversion
c = 299792458.0  # Speed of light in m/s

# region
###################################################################


def read_model_coordinates(points):

    # read PEC nodes
    xpts = points[:, 0]
    ypts = points[:, 1]
    zpts = points[:, 2]
    nverts = len(xpts)
    return xpts, ypts, zpts, nverts

###################################################################


def generate_transpose_matrix(facets):
    node1 = facets[:, 0]
    node2 = facets[:, 1]
    node3 = facets[:, 2]
    return node1, node2, node3

###################################################################
# endregion


""" 
Aperture Integration with Kim
"""
# region


def ApertureIntegrationKim(RayTubes, tstart_obs, tstop_obs, delt, pstart_obs, pstop_obs, delp, k, rho):

    totalTubes = RayTubes.n_faces

    vind = np.zeros(shape=(totalTubes, 3), dtype='int')
    beta = np.zeros(totalTubes)
    alpha = np.zeros(totalTubes)

    Efield = np.zeros(shape=(3, totalTubes), dtype=complex)
    Hfield = np.zeros(shape=(3, totalTubes), dtype=complex)
    ki = np.zeros(shape=(3, totalTubes))
    modE = np.zeros(totalTubes)
    Phase = np.zeros(totalTubes)

    AthetaDebug = np.zeros(totalTubes,dtype=complex)
    AphiDebug = np.zeros(totalTubes, dtype=complex)
    ShapeDebug = np.zeros(totalTubes, dtype=complex)

    ############################
    # Set up observation angles

    if tstart_obs == tstop_obs:
        thr0 = tstart_obs*rad

    if pstart_obs == pstop_obs:
        phr0 = pstart_obs*rad

    # steps in theta and in phi for the observation direction

    it = math.floor((tstop_obs-tstart_obs)/delt) + 1
    ip = math.floor((pstop_obs-pstart_obs)/delp) + 1

    phi = np.zeros(shape=(ip, it))
    theta = np.zeros(shape=(ip, it))

    ############################
    # Tubes quantities

    xiyizi = RayTubes.cell_centers()
    Normal = RayTubes['Normal']
    Area = RayTubes['Area']

    xpts, ypts, zpts, nverts = read_model_coordinates(RayTubes.points)

    node1, node2, node3 = generate_transpose_matrix(
        RayTubes.faces.reshape(len(RayTubes.faces)//4, 4)[:, 1:4])

    for i in range(totalTubes):
        pts = np.array([node1[i], node2[i], node3[i]])
        vind[i, :] = pts

    x = xpts
    y = ypts
    z = zpts

    for i in range(totalTubes):

        """ 
        turn polydata normals if K*N is negative (inward-outward)
        """

        dot = np.dot(Normal[i], RayTubes['K'][i])
        if dot < 0:
            Normal[i] = -Normal[i]

        beta[i] = np.arccos(Normal[i, 2])
        alpha[i] = np.arctan2(Normal[i, 1], Normal[i, 0])

    ############################
    # EM field quantities

    Z0 = 1/(epsilon0*c)

    for m in range(totalTubes):

        ki[:, m] = RayTubes['K'][m]
        # modE[m] = RayTubes['|E|'][m]
        Phase[m] = RayTubes['phase'][m]

        # we are transforming in complex vector
        Efield[:, m] = RayTubes['Efield'][m]

        Efield[np.isnan(Efield[:, m]),m]=0
        Efield[np.isinf(Efield[:, m]),m]=0
        
        Hfield[:, m] = 1/Z0 * np.cross(ki[:, m], Efield[:, m])

    ############################
    # Start of angles Loop

    totSteps = it*ip
    counterSteps = 0

    EScat = np.zeros(shape=(2, ip, it), dtype=complex)

    # for i1 in tqdm(range(ip)):

    #     for i2 in tqdm(range(it)):
    
    for i1, i2 in tqdm(itertools.product(range(ip), range(it)), total=ip * it, desc="Processing", leave=True):

            counterSteps = counterSteps + 1

            phi[i1, i2] = pstart_obs + i1 * delp
            theta[i1, i2] = tstart_obs + i2 * delt

            phr = phi[i1, i2] * rad
            thr = theta[i1, i2] * rad

            # print('*** step = ', counterSteps, ' / ', totSteps,
            #       ' -- (Theta : ', theta[i1, i2], ' Phy : ', phi[i1, i2], ') \n')

            st = math.sin(thr)
            ct = math.cos(thr)
            cp = math.cos(phr)
            sp = math.sin(phr)

            # observation direction
            r_versor = np.array([st*cp, st*sp, ct])
            theta_versor = np.array([ct*cp, ct*sp, -st])
            phi_versor = np.array([-sp, cp, 0])

            theta_versor = -theta_versor
            phi_versor = -phi_versor

            # r versor components
            u = st*cp
            v = st*sp
            w = ct

            # theta versor components
            uu = ct*cp
            vv = ct*sp
            ww = -st

            # Transformation matrix
            Tspherical = [[st*cp, st*sp, ct],
                          [ct*cp, ct*sp, -st],
                          [-sp, cp, 0]]

            # observation direction
            R_obs = np.array([u, v, w])

            sumAtheta = 0
            sumAphi = 0

            ############################
            """
            #Tubes Aperture
            """

            for iTube in range(totalTubes):

                _xiyizi = xiyizi.points[iTube]
                _ki = ki[:, iTube]
                _Efield = Efield[:, iTube]
                _Phase = Phase[iTube]
                _modE = modE[iTube]
                _Normal = Normal[iTube, :]
                _vind = vind[iTube, :]
                _Area = Area[iTube]

                # Tube direction
                Ri = np.array([_ki[0], _ki[1], _ki[2]])

                Ic = ShapeFactor.GisbonShapeFactorModified(
                    k, x, y, z, _vind, -Ri, R_obs, _Area, _xiyizi)
                
                # Ic = ShapeFactor.GisbonShapeFactor(
                #     k, x, y, z, _vind, -Ri, R_obs, _Area)
                
                _Shape = complex(np.real(Ic), -np.imag(Ic))
                if np.real(_Shape) == 0:
                    _Shape
                    
                # _Shape = Ic

                # _Shape = Ic
                # if _Area == 0:
                #     _Shape =0
                # else:
                #     _Shape=_Shape/_Area
                radiusTube = np.sqrt(_xiyizi[0]**2+_xiyizi[1]**2+_xiyizi[2]**2)
                
                if radiusTube < 14:
                    _Shape = 0
                    
                    
                _jK4pi = 1j*k/(4*np.pi)

                _PhaseTerm = np.exp(1j*k*(np.dot(R_obs, _xiyizi)))

                _Ei = _Efield*np.exp(-1j*k*_Phase)*np.exp(-1j*k*rho)

                _Hi = np.cross(_ki, _Ei)

                _aTheta = np.dot(np.cross(-phi_versor, _Ei) +
                                 np.cross(theta_versor, _Hi), _Normal)

                _aPhi = np.dot(np.cross(theta_versor, _Ei) +
                               np.cross(phi_versor, _Hi), _Normal)

                _Atheta = _jK4pi*_PhaseTerm*_aTheta*_Shape
                _Aphi = _jK4pi*_PhaseTerm*_aPhi*_Shape

                sumAtheta = sumAtheta+_Atheta
                sumAphi = sumAphi+_Aphi
                
                ShapeDebug[iTube] = _Shape
                AthetaDebug[iTube] = _Atheta
                AphiDebug[iTube] = _Aphi

            EScat[0, i1, i2] = sumAtheta
            EScat[1, i1, i2] = sumAphi

    return EScat, theta, phi

# endregion


###################################################################
""" 
# Radar Cross Section
"""


def RCS(RayTubes, tstart_obs, tstop_obs, delt, pstart_obs, pstop_obs, delp, k, rho):

    totalTubes = RayTubes.n_faces

    _lambda = 2*np.pi/k

    vind = np.zeros(shape=(totalTubes, 3), dtype='int')
    beta = np.zeros(totalTubes)
    alpha = np.zeros(totalTubes)

    Efield = np.zeros(shape=(3, totalTubes), dtype=complex)
    Hfield = np.zeros(shape=(3, totalTubes), dtype=complex)
    ki = np.zeros(shape=(3, totalTubes))
    modE = np.zeros(totalTubes)
    Phase = np.zeros(totalTubes)

    ############################
    # Set up observation angles

    if tstart_obs == tstop_obs:
        thr0 = tstart_obs*rad

    if pstart_obs == pstop_obs:
        phr0 = pstart_obs*rad

    # steps in theta and in phi for the observation direction

    it = math.floor((tstop_obs-tstart_obs)/delt) + 1
    ip = math.floor((pstop_obs-pstart_obs)/delp) + 1

    phi = np.zeros(shape=(ip, it))
    theta = np.zeros(shape=(ip, it))

    ############################
    # Tubes quantities

    xiyizi = RayTubes.cell_centers()
    Normal = RayTubes['Normal']
    Area = RayTubes['Area']

    xpts, ypts, zpts, nverts = read_model_coordinates(RayTubes.points)

    node1, node2, node3 = generate_transpose_matrix(
        RayTubes.faces.reshape(len(RayTubes.faces)//4, 4)[:, 1:4])

    for i in range(totalTubes):
        pts = np.array([node1[i], node2[i], node3[i]])
        vind[i, :] = pts

    x = xpts
    y = ypts
    z = zpts

    for i in range(totalTubes):

        """ 
        turn polydata normals if K*N is negative (inward-outward)
        """

        dot = np.dot(Normal[i], RayTubes['K'][i])
        if dot < 0:
            Normal[i] = -Normal[i]

        beta[i] = np.arccos(Normal[i, 2])
        alpha[i] = np.arctan2(Normal[i, 1], Normal[i, 0])

    ############################
    # EM field quantities

    Z0 = 1/(epsilon0*c)

    for m in range(totalTubes):

        ki[:, m] = RayTubes['K'][m]
        # modE[m] = RayTubes['|E|'][m]
        Phase[m] = RayTubes['phase'][m]

        # we are transforming in complex vector

        Efield[:, m] = RayTubes['Efield'][m]
        Hfield[:, m] = 1/Z0 * np.cross(ki[:, m], Efield[:, m])

    ############################
    # Start of angles Loop

    totSteps = it*ip
    counterSteps = 0

    EScat = np.zeros(shape=(2, ip, it), dtype=complex)
    RCSth = np.zeros(shape=(2, ip, it))

    for i1 in range(ip):

        for i2 in range(it):

            counterSteps = counterSteps + 1

            phi[i1, i2] = pstart_obs + i1 * delp
            theta[i1, i2] = tstart_obs + i2 * delt

            phr = phi[i1, i2] * rad
            thr = theta[i1, i2] * rad

            print('*** step = ', counterSteps, ' / ', totSteps,
                  ' -- (Theta : ', theta[i1, i2], ' Phy : ', phi[i1, i2], ') \n')

            st = math.sin(thr)
            ct = math.cos(thr)
            cp = math.cos(phr)
            sp = math.sin(phr)

            # observation direction
            r_versor = np.array([st*cp, st*sp, ct])
            theta_versor = np.array([ct*cp, ct*sp, -st])
            phi_versor = np.array([-sp, cp, 0])

            theta_versor = -theta_versor
            phi_versor = -phi_versor

            # r versor components
            u = st*cp
            v = st*sp
            w = ct

            # theta versor components
            uu = ct*cp
            vv = ct*sp
            ww = -st

            # Transformation matrix
            Tspherical = [[st*cp, st*sp, ct],
                          [ct*cp, ct*sp, -st],
                          [-sp, cp, 0]]

            # observation direction
            R_obs = np.array([u, v, w])

            sumAtheta = 0
            sumAphi = 0

            ############################
            """
            #Tubes Aperture
            """

            for iTube in range(totalTubes):

                _xiyizi = xiyizi.points[iTube]
                _ki = ki[:, iTube]
                _Efield = Efield[:, iTube]
                _Phase = Phase[iTube]
                _modE = modE[iTube]
                _Normal = Normal[iTube, :]
                _vind = vind[iTube, :]
                _Area = Area[iTube]

                # Tube direction
                Ri = np.array([_ki[0], _ki[1], _ki[2]])

                Ic = ShapeFactor.GisbonShapeFactorModified(
                    k, x, y, z, _vind, -Ri, R_obs, _Area, _xiyizi)

                _Shape = complex(np.real(Ic), -np.imag(Ic))

                # _Shape=_Shape/_Area

                _jK4pi = 1j*k/(4*np.pi)

                _PhaseTerm = np.exp(1j*k*(np.dot(R_obs, _xiyizi)))

                _Ei = _Efield*np.exp(-1j*k*_Phase)*np.exp(-1j*k*rho)

                _Hi = np.cross(_ki, _Ei)

                _aTheta = np.dot(np.cross(-phi_versor, _Ei) +
                                 np.cross(theta_versor, _Hi), _Normal)

                _aPhi = np.dot(np.cross(theta_versor, _Ei) +
                               np.cross(phi_versor, _Hi), _Normal)

                _Atheta = _jK4pi*_PhaseTerm*_aTheta*_Shape
                _Aphi = _jK4pi*_PhaseTerm*_aPhi*_Shape

                sumAtheta = sumAtheta+_Atheta
                sumAphi = sumAphi+_Aphi

            EScat[0, i1, i2] = sumAtheta
            EScat[1, i1, i2] = sumAphi
            RCSth[0, i1, i2] = 10*np.log10(4*np.pi*(abs(sumAtheta)**2)+1e-10)
            RCSth[1, i1, i2] = 10*np.log10(4*np.pi*(abs(sumAphi)**2)+1e-10)
            # RCSth[0, i1, i2] = 10*np.log10(abs(sumAtheta)+1e-10)
            # RCSth[1, i1, i2] = 10*np.log10(abs(sumAphi)+1e-10)

    return EScat, RCSth, theta, phi

###################################################################


def RCSPO(RayTubes, tstart_obs, tstop_obs, delt, pstart_obs, pstop_obs, delp, k, rho):

    totalTubes = RayTubes.n_faces

    _lambda = 2*np.pi/k

    vind = np.zeros(shape=(totalTubes, 3), dtype='int')
    beta = np.zeros(totalTubes)
    alpha = np.zeros(totalTubes)

    Efield = np.zeros(shape=(3, totalTubes), dtype=complex)
    Hfield = np.zeros(shape=(3, totalTubes), dtype=complex)
    ki = np.zeros(shape=(3, totalTubes))
    modE = np.zeros(totalTubes)
    Phase = np.zeros(totalTubes)

    ############################
    # Set up observation angles

    if tstart_obs == tstop_obs:
        thr0 = tstart_obs*rad

    if pstart_obs == pstop_obs:
        phr0 = pstart_obs*rad

    # steps in theta and in phi for the observation direction

    it = math.floor((tstop_obs-tstart_obs)/delt) + 1
    ip = math.floor((pstop_obs-pstart_obs)/delp) + 1

    phi = np.zeros(shape=(ip, it))
    theta = np.zeros(shape=(ip, it))

    ############################
    # Tubes quantities

    xiyizi = RayTubes.cell_centers()
    Normal = RayTubes['Normal']
    Area = RayTubes['Area']

    xpts, ypts, zpts, nverts = read_model_coordinates(RayTubes.points)

    node1, node2, node3 = generate_transpose_matrix(
        RayTubes.faces.reshape(len(RayTubes.faces)//4, 4)[:, 1:4])

    for i in range(totalTubes):
        pts = np.array([node1[i], node2[i], node3[i]])
        vind[i, :] = pts

    x = xpts
    y = ypts
    z = zpts

    for i in range(totalTubes):

        """ 
        turn polydata normals if K*N is negative (inward-outward)
        """

        dot = np.dot(Normal[i], RayTubes['K'][i])
        if dot < 0:
            Normal[i] = -Normal[i]

        beta[i] = np.arccos(Normal[i, 2])
        alpha[i] = np.arctan2(Normal[i, 1], Normal[i, 0])

    ############################
    # EM field quantities

    Z0 = 1/(epsilon0*c)

    for m in range(totalTubes):

        ki[:, m] = RayTubes['K'][m]
        # modE[m] = RayTubes['|E|'][m]
        Phase[m] = RayTubes['phase'][m]

        # we are transforming in complex vector

        Efield[:, m] = RayTubes['Efield'][m]
        Hfield[:, m] = 1/Z0 * np.cross(ki[:, m], Efield[:, m])

    ############################
    # Start of angles Loop

    totSteps = it*ip
    counterSteps = 0

    EScat = np.zeros(shape=(2, ip, it), dtype=complex)
    RCSth = np.zeros(shape=(2, ip, it))

    for i1 in range(ip):

        for i2 in range(it):

            counterSteps = counterSteps + 1

            phi[i1, i2] = pstart_obs + i1 * delp
            theta[i1, i2] = tstart_obs + i2 * delt

            phr = phi[i1, i2] * rad
            thr = theta[i1, i2] * rad

            print('*** step = ', counterSteps, ' / ', totSteps,
                  ' -- (Theta : ', theta[i1, i2], ' Phy : ', phi[i1, i2], ') \n')

            st = math.sin(thr)
            ct = math.cos(thr)
            cp = math.cos(phr)
            sp = math.sin(phr)

            # observation direction
            r_versor = np.array([st*cp, st*sp, ct])
            theta_versor = np.array([ct*cp, ct*sp, -st])
            phi_versor = np.array([-sp, cp, 0])

            # theta_versor = -theta_versor
            # phi_versor = -phi_versor

            # r versor components
            u = st*cp
            v = st*sp
            w = ct

            # theta versor components
            uu = ct*cp
            vv = ct*sp
            ww = -st

            # Transformation matrix
            Tspherical = [[st*cp, st*sp, ct],
                          [ct*cp, ct*sp, -st],
                          [-sp, cp, 0]]

            # observation direction
            R_obs = np.array([u, v, w])

            sumAtheta = 0
            sumAphi = 0

            ############################
            """
            #Tubes Aperture
            """

            for iTube in range(totalTubes):

                _xiyizi = xiyizi.points[iTube]
                _ki = ki[:, iTube]
                _Efield = Efield[:, iTube]
                _Phase = Phase[iTube]
                _modE = modE[iTube]
                _Normal = Normal[iTube, :]
                _vind = vind[iTube, :]
                _Area = Area[iTube]

                # Tube direction
                Ri = np.array([_ki[0], _ki[1], _ki[2]])

                Ic = ShapeFactor.GisbonShapeFactorModified(
                    k, x, y, z, _vind, -Ri, R_obs, _Area, _xiyizi)

                # Ic= ShapeFactor.GisbonShapeFactor(k, x, y, z, _vind, -Ri, R_obs, _Area)

                radiusTube = np.sqrt(_xiyizi[0]**2+_xiyizi[1]**2+_xiyizi[2]**2)
                
                if radiusTube < 12:
                    Ic=0
                    
                
                _Shape = complex(np.real(Ic), -np.imag(Ic))


                _jK4pi = 1j*k/(4*np.pi)

                _PhaseTerm = np.exp(1j*k*(np.dot(R_obs, _xiyizi)))

                _Ei = _Efield*np.exp(-1j*k*_Phase)*np.exp(-1j*k*rho)

                _Hi = np.cross(_ki, _Ei)

                _aTheta = np.dot(np.cross(-phi_versor, _Ei) +
                                 np.cross(theta_versor, _Hi), _Normal)

                _aPhi = np.dot(np.cross(theta_versor, _Ei) +
                               np.cross(phi_versor, _Hi), _Normal)

                _Atheta = _jK4pi*_PhaseTerm*_aTheta*_Shape
                _Aphi = _jK4pi*_PhaseTerm*_aPhi*_Shape

                sumAtheta = sumAtheta+_Atheta
                sumAphi = sumAphi+_Aphi

            EScat[0, i1, i2] = sumAtheta
            EScat[1, i1, i2] = sumAphi
            RCSth[0, i1, i2] = 10*np.log10(4*np.pi*(abs(sumAtheta)**2)+1e-10)
            RCSth[1, i1, i2] = 10*np.log10(4*np.pi*(abs(sumAphi)**2)+1e-10)
            # RCSth[0, i1, i2] = 10*np.log10(4*np.pi*(abs(sumAtheta)**2)/(_lambda**2)+1e-10)
            # RCSth[1, i1, i2] = 10*np.log10(4*np.pi*(abs(sumAphi)**2)/(_lambda**2)+1e-10)
            # RCSth[0, i1, i2] = 10*np.log10(abs(sumAtheta)**2)
            # RCSth[1, i1, i2] = 10*np.log10(abs(sumAphi)**2)

    return EScat, RCSth, theta, phi
