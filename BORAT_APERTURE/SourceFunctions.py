import pyvista as pv
import numpy as np
import meshzoo
import math

import Functions
import EikonalSolver
import pyacvd

import GlobalData


###################################################################

def Antenna(Solver):

    a_location = Solver.Antenna_Location
    rho = Solver.Rho
    iPol = Solver.iPol
    E0_amplitude = Solver.E0_Amplitude
    rotationAngleDeg = Solver.Antenna_Direction
    Aperture = Solver.Antenna_Aperture
    TotalRays = Solver.TotalRays
    rotationAxis = Solver.AntennaAxis

    div = TotalRays//20
    InitialConditions = []
    e0 = []
    e0_spherical = []
    E0_sperical = []
    E0_rho = []

    n0, _, _, _ = EikonalSolver.RefractiveIndexInterpolation(
        [a_location[0], a_location[1], a_location[2]])

    points, cells = meshzoo.icosa_sphere(TotalRays)

    arrayNodeMap = np.zeros(shape=(cells.shape[0], 4), dtype='<i4')
    arrayNodeMap[:, 0] = 3
    arrayNodeMap[:, 1:4] = cells

    icoSphere0 = pv.PolyData(points, arrayNodeMap.ravel())

    originCut = [0, 0, np.sin((90-Aperture)*np.pi/180)]

    icoSphere0 = icoSphere0.clip(normal='z', origin=originCut, invert=False)

    icoSphere0 = icoSphere0.rotate_vector((0, 1, 0), 90, inplace=True)

    icoSphere0 = icoSphere0.scale([rho, rho, rho], inplace=True)

    # clus = pyacvd.Clustering(icoSphere0)
    # # mesh is not dense enough for uniform remeshing
    # clus.subdivide(4)
    # clus.cluster(TotalRays)

    # # plot clustered cow mesh
    # # clus.plot()

    # icoSphere = clus.create_mesh()
    
    icoSphere = icoSphere0
    
    points = icoSphere.points

    cells = icoSphere.faces

    Normal = icoSphere.points/rho

    totRay = icoSphere.n_points
    RayDir = np.zeros(shape=(icoSphere.n_points, 3))

    for nRay in range(totRay):

        _, theta, phi = Functions.Cart2sphere(
            Normal[nRay, 0], Normal[nRay, 1], Normal[nRay, 2])

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        RayDir[nRay, 0] = r_versor[0]*n0
        RayDir[nRay, 1] = r_versor[1]*n0
        RayDir[nRay, 2] = r_versor[2]*n0

        if iPol == 1:
            e = np.array([np.cos(theta)*np.cos(phi),
                          np.cos(theta) * np.sin(phi),
                          -np.sin(theta)])
        else:

            e = np.array([-np.sin(phi),
                          np.cos(phi),
                          0])

        e_spherical = np.array([np.dot(e, r_versor),
                                np.dot(e, theta_versor),
                                np.dot(e, phi_versor)])

        e0.append(e)
        e0_spherical.append(e_spherical)

        Edir = -e*np.sin(theta)

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        # E0_amplitude_rays=E0_amplitude*1/rho
        E0_amplitude_rays = E0_amplitude*1/rho

        E0_rho.append(E0_amplitude_rays*Edir)

        E0_sperical.append(E0_amplitude_rays*Edir_spherical)

    E0_rho = np.array(E0_rho)
    e0 = np.array(e0)
    e0_spherical = np.array(e0_spherical)
    E0_sperical = np.array(E0_sperical)

    Ray0 = pv.PolyData(points)

    Ray0.point_data['kx'] = RayDir[:, 0]
    Ray0.point_data['ky'] = RayDir[:, 1]
    Ray0.point_data['kz'] = RayDir[:, 2]

    Ray0.point_data['ex'] = e0[:, 0]
    Ray0.point_data['ey'] = e0[:, 1]
    Ray0.point_data['ez'] = e0[:, 2]

    Ray0['er'] = e0_spherical[:, 0]
    Ray0['et'] = e0_spherical[:, 1]
    Ray0['ep'] = e0_spherical[:, 2]

    Ray0['K'] = np.column_stack([Ray0['kx'], Ray0['ky'], Ray0['kz']])

    Ray0['e0'] = np.column_stack([Ray0['ex'], Ray0['ey'], Ray0['ez']])
    Ray0['e0_spherical'] = np.column_stack(
        [Ray0['er'], Ray0['et'], Ray0['ep']])

    Ray0['E0'] = np.column_stack(
        [(E0_rho[:, 0]), (E0_rho[:, 1]), (E0_rho[:, 2])])

    Ray0['E0x'] = E0_rho[:, 0]
    Ray0['E0y'] = E0_rho[:, 1]
    Ray0['E0z'] = E0_rho[:, 2]

    Ray0['E0r'] = E0_sperical[:, 0]
    Ray0['E0t'] = E0_sperical[:, 1]
    Ray0['E0p'] = E0_sperical[:, 2]

    Ray0['|E0|'] = np.sqrt(Ray0['E0x']**2+Ray0['E0y']**2+Ray0['E0z']**2)

    Ray0['E0_spherical'] = np.column_stack(
        [Ray0['E0r'], Ray0['E0t'], Ray0['E0p']])

    RayTube0 = pv.PolyData(points, cells)

    totTubes = cells.shape[0]//4

    e0_tube = []
    e0_tube_spherical = []
    E0_tube = []
    E0_tube_spherical = []
    ATube0 = []
    thetaTube = []
    phiTube = []

    rVersor_cell = RayTube0.face_normals
    RayTubes0 = RayTube0.compute_cell_sizes(
        length=False, area=True, volume=False)

    rho_centers = RayTube0.cell_centers()

    for nTubes in range(totTubes):

        ATube0.append(RayTubes0['Area'][nTubes])

        _, theta, phi = Functions.Cart2sphere(
            rVersor_cell[nTubes][0], rVersor_cell[nTubes][1], rVersor_cell[nTubes][2])

        thetaTube.append(theta)
        phiTube.append(phi)

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        e = theta_versor
        e_spherical = np.array([np.dot(e, r_versor), np.dot(
            e, theta_versor), np.dot(e, phi_versor)])

        e0_tube.append(e)
        e0_tube_spherical.append(e_spherical)

        Edir = -e*np.sin(theta)

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        # E0_Tube = E0_amplitude*1/(np.linalg.norm(rho_centers.points[nTubes]))
        E0_Tube = E0_amplitude*1/rho

        E0_tube.append(E0_Tube * Edir)

        E0_tube_spherical.append(E0_Tube*Edir_spherical)

    e0_tube = np.array(e0_tube)
    e0_tube_spherical = np.array(e0_tube_spherical)
    E0_tube = np.array(E0_tube)
    E0_tube_spherical = np.array(E0_tube_spherical)
    ATube0 = np.array(ATube0)
    thetaTube = np.array(thetaTube)
    phiTube = np.array(phiTube)

    RayTube0.cell_data['kx'] = rVersor_cell[:, 0]
    RayTube0.cell_data['ky'] = rVersor_cell[:, 1]
    RayTube0.cell_data['kz'] = rVersor_cell[:, 2]

    RayTube0.cell_data['ex'] = e0_tube[:, 0]
    RayTube0.cell_data['ey'] = e0_tube[:, 1]
    RayTube0.cell_data['ez'] = e0_tube[:, 2]

    RayTube0.cell_data['er'] = e0_tube_spherical[:, 0]
    RayTube0.cell_data['et'] = e0_tube_spherical[:, 1]
    RayTube0.cell_data['ep'] = e0_tube_spherical[:, 2]

    RayTube0.cell_data['E0x'] = E0_tube[:, 0]
    RayTube0.cell_data['E0y'] = E0_tube[:, 1]
    RayTube0.cell_data['E0z'] = E0_tube[:, 2]

    RayTube0.cell_data['E0r'] = E0_tube_spherical[:, 0]
    RayTube0.cell_data['E0t'] = E0_tube_spherical[:, 1]
    RayTube0.cell_data['E0p'] = E0_tube_spherical[:, 2]

    RayTube0.cell_data['K'] = np.column_stack(
        [RayTube0['kx'], RayTube0['ky'], RayTube0['kz']])

    RayTube0.cell_data['e0'] = np.column_stack(
        [RayTube0['ex'], RayTube0['ey'], RayTube0['ez']])
    RayTube0.cell_data['e0_spherical'] = np.column_stack(
        [RayTube0['er'], RayTube0['et'], RayTube0['ep']])

    RayTube0.cell_data['E0'] = np.column_stack(
        [RayTube0.cell_data['E0x'], RayTube0.cell_data['E0y'], RayTube0.cell_data['E0z']])

    RayTube0.cell_data['E0_spherical'] = np.column_stack(
        [RayTube0.cell_data['E0r'], RayTube0.cell_data['E0t'], RayTube0.cell_data['E0p']])

    RayTube0.cell_data['|E0|'] = np.sqrt(
        RayTube0['E0x']**2+RayTube0['E0y']**2+RayTube0['E0z']**2)

    RayTube0.cell_data['Area'] = ATube0

    RayTube0.rotate_vector(rotationAxis, -rotationAngleDeg,
                           inplace=True, transform_all_input_vectors=True)

    Ray0.rotate_vector(rotationAxis, -rotationAngleDeg,
                       inplace=True, transform_all_input_vectors=True)

    RayTube0 = RayTube0.translate(
        (a_location[0], a_location[1], a_location[2]), inplace=True, transform_all_input_vectors=True)

    Ray0 = Ray0.translate(
        (a_location[0], a_location[1], a_location[2]), inplace=True, transform_all_input_vectors=True)

    points = Ray0.points

    for nRay in range(totRay):

        InitialConditions.append([points[nRay][0], points[nRay][1], points[nRay][2],
                                  Ray0['K'][nRay][0], Ray0['K'][nRay][1], Ray0['K'][nRay][2],
                                  Ray0['e0'][nRay][0], Ray0['e0'][nRay][1], Ray0['e0'][nRay][2],
                                  0, 0])     # Change these initial conditions

    return InitialConditions, Ray0, RayTube0

###################################################################

def Antenna_rotated(rho, a_location, E0_amplitude, TotalRays, rotationAngle, rotationAxis, Aperture, iPol):

    # a_location = Solver.Antenna_Location
    # rho = Solver.Rho
    # iPol = Solver.iPol
    # E0_amplitude = Solver.E0_Amplitude
    # NormalAntenna = Solver.Antenna_Direction
    # Aperture = Solver.Antenna_Aperture
    # TotalRays = Solver.TotalRays

    InitialConditions = []
    e0 = []
    e0_spherical = []
    E0_sperical = []
    E0_rho = []

    n0 = 1

    points, cells = meshzoo.icosa_sphere(20)

    arrayNodeMap = np.zeros(shape=(cells.shape[0], 4), dtype='<i4')
    arrayNodeMap[:, 0] = 3
    arrayNodeMap[:, 1:4] = cells

    icoSphere0 = pv.PolyData(points, arrayNodeMap.ravel())

    originCut = [0, 0, np.sin((90-Aperture)*np.pi/180)]

    icoSphere0 = icoSphere0.clip(normal='z', origin=originCut, invert=False)

    icoSphere0 = icoSphere0.rotate_vector((0, 1, 0), 90, inplace=True)

    icoSphere0 = icoSphere0.scale([rho, rho, rho], inplace=True)

    clus = pyacvd.Clustering(icoSphere0)
    # mesh is not dense enough for uniform remeshing
    clus.subdivide(4)
    clus.cluster(TotalRays)

    # plot clustered cow mesh
    # clus.plot()

    icoSphere = clus.create_mesh()

    points = icoSphere.points

    cells = icoSphere.faces

    Normal = icoSphere.points/rho

    totRay = icoSphere.n_points
    RayDir = np.zeros(shape=(icoSphere.n_points, 3))

    for nRay in range(totRay):

        _, theta, phi = Functions.Cart2sphere(
            Normal[nRay, 0], Normal[nRay, 1], Normal[nRay, 2])

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        RayDir[nRay, 0] = r_versor[0]*n0
        RayDir[nRay, 1] = r_versor[1]*n0
        RayDir[nRay, 2] = r_versor[2]*n0

        if iPol == 1:
            e = np.array([np.cos(theta)*np.cos(phi),
                          np.cos(theta) * np.sin(phi),
                          -np.sin(theta)])
        else:

            e = np.array([-np.sin(phi),
                          np.cos(phi),
                          0])

        e_spherical = np.array([np.dot(e, r_versor),
                                np.dot(e, theta_versor),
                                np.dot(e, phi_versor)])

        e0.append(e)
        e0_spherical.append(e_spherical)

        Edir = -e*np.sin(theta)

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        # E0_amplitude_rays=E0_amplitude*1/rho
        E0_amplitude_rays = E0_amplitude*1/rho

        E0_rho.append(E0_amplitude_rays*Edir)

        E0_sperical.append(E0_amplitude_rays*Edir_spherical)

    E0_rho = np.array(E0_rho)
    e0 = np.array(e0)
    e0_spherical = np.array(e0_spherical)
    E0_sperical = np.array(E0_sperical)

    Ray0 = pv.PolyData(points)

    Ray0.point_data['kx'] = RayDir[:, 0]
    Ray0.point_data['ky'] = RayDir[:, 1]
    Ray0.point_data['kz'] = RayDir[:, 2]

    Ray0.point_data['ex'] = e0[:, 0]
    Ray0.point_data['ey'] = e0[:, 1]
    Ray0.point_data['ez'] = e0[:, 2]

    Ray0['er'] = e0_spherical[:, 0]
    Ray0['et'] = e0_spherical[:, 1]
    Ray0['ep'] = e0_spherical[:, 2]

    Ray0['K'] = np.column_stack([Ray0['kx'], Ray0['ky'], Ray0['kz']])

    Ray0['e0'] = np.column_stack([Ray0['ex'], Ray0['ey'], Ray0['ez']])
    Ray0['e0_spherical'] = np.column_stack(
        [Ray0['er'], Ray0['et'], Ray0['ep']])

    Ray0['E0'] = np.column_stack(
        [(E0_rho[:, 0]), (E0_rho[:, 1]), (E0_rho[:, 2])])

    Ray0['E0x'] = E0_rho[:, 0]
    Ray0['E0y'] = E0_rho[:, 1]
    Ray0['E0z'] = E0_rho[:, 2]

    Ray0['E0r'] = E0_sperical[:, 0]
    Ray0['E0t'] = E0_sperical[:, 1]
    Ray0['E0p'] = E0_sperical[:, 2]

    Ray0['|E0|'] = np.sqrt(Ray0['E0x']**2+Ray0['E0y']**2+Ray0['E0z']**2)

    Ray0['E0_spherical'] = np.column_stack(
        [Ray0['E0r'], Ray0['E0t'], Ray0['E0p']])

    RayTube0 = pv.PolyData(points, cells)

    totTubes = cells.shape[0]//4

    e0_tube = []
    e0_tube_spherical = []
    E0_tube = []
    E0_tube_spherical = []
    ATube0 = []
    thetaTube = []
    phiTube = []

    rVersor_cell = RayTube0.face_normals
    RayTubes0 = RayTube0.compute_cell_sizes(
        length=False, area=True, volume=False)

    rho_centers = RayTube0.cell_centers()

    for nTubes in range(totTubes):

        ATube0.append(RayTubes0['Area'][nTubes])

        _, theta, phi = Functions.Cart2sphere(
            rVersor_cell[nTubes][0], rVersor_cell[nTubes][1], rVersor_cell[nTubes][2])

        thetaTube.append(theta)
        phiTube.append(phi)

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        e = theta_versor
        e_spherical = np.array([np.dot(e, r_versor), np.dot(
            e, theta_versor), np.dot(e, phi_versor)])

        e0_tube.append(e)
        e0_tube_spherical.append(e_spherical)

        Edir = -e*np.sin(theta)

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        # E0_Tube = E0_amplitude*1/(np.linalg.norm(rho_centers.points[nTubes]))
        E0_Tube = E0_amplitude*1/rho

        E0_tube.append(E0_Tube * Edir)

        E0_tube_spherical.append(E0_Tube*Edir_spherical)

    e0_tube = np.array(e0_tube)
    e0_tube_spherical = np.array(e0_tube_spherical)
    E0_tube = np.array(E0_tube)
    E0_tube_spherical = np.array(E0_tube_spherical)
    ATube0 = np.array(ATube0)
    thetaTube = np.array(thetaTube)
    phiTube = np.array(phiTube)

    RayTube0.cell_data['kx'] = rVersor_cell[:, 0]
    RayTube0.cell_data['ky'] = rVersor_cell[:, 1]
    RayTube0.cell_data['kz'] = rVersor_cell[:, 2]

    RayTube0.cell_data['ex'] = e0_tube[:, 0]
    RayTube0.cell_data['ey'] = e0_tube[:, 1]
    RayTube0.cell_data['ez'] = e0_tube[:, 2]

    RayTube0.cell_data['er'] = e0_tube_spherical[:, 0]
    RayTube0.cell_data['et'] = e0_tube_spherical[:, 1]
    RayTube0.cell_data['ep'] = e0_tube_spherical[:, 2]

    RayTube0.cell_data['E0x'] = E0_tube[:, 0]
    RayTube0.cell_data['E0y'] = E0_tube[:, 1]
    RayTube0.cell_data['E0z'] = E0_tube[:, 2]

    RayTube0.cell_data['E0r'] = E0_tube_spherical[:, 0]
    RayTube0.cell_data['E0t'] = E0_tube_spherical[:, 1]
    RayTube0.cell_data['E0p'] = E0_tube_spherical[:, 2]

    RayTube0.cell_data['K'] = np.column_stack(
        [RayTube0['kx'], RayTube0['ky'], RayTube0['kz']])

    RayTube0.cell_data['e0'] = np.column_stack(
        [RayTube0['ex'], RayTube0['ey'], RayTube0['ez']])
    RayTube0.cell_data['e0_spherical'] = np.column_stack(
        [RayTube0['er'], RayTube0['et'], RayTube0['ep']])

    RayTube0.cell_data['E0'] = np.column_stack(
        [RayTube0.cell_data['E0x'], RayTube0.cell_data['E0y'], RayTube0.cell_data['E0z']])

    RayTube0.cell_data['E0_spherical'] = np.column_stack(
        [RayTube0.cell_data['E0r'], RayTube0.cell_data['E0t'], RayTube0.cell_data['E0p']])

    RayTube0.cell_data['|E0|'] = np.sqrt(
        RayTube0['E0x']**2+RayTube0['E0y']**2+RayTube0['E0z']**2)

    RayTube0.cell_data['Area'] = ATube0

    # RayTube0.rotate_vector([1,0,0], 90,
    #                        inplace=True, transform_all_input_vectors=True)

    # Ray0.rotate_vector([1,0,0], 90,
    #                    inplace=True, transform_all_input_vectors=True)

    RayTube0.rotate_vector(rotationAxis, -rotationAngle,
                           inplace=True, transform_all_input_vectors=True)

    Ray0.rotate_vector(rotationAxis, -rotationAngle,
                       inplace=True, transform_all_input_vectors=True)

    RayTube0 = RayTube0.translate(
        (a_location[0], a_location[1], a_location[2]), inplace=True, transform_all_input_vectors=True)

    Ray0 = Ray0.translate(
        (a_location[0], a_location[1], a_location[2]), inplace=True, transform_all_input_vectors=True)

    points = Ray0.points

    for nRay in range(totRay):

        InitialConditions.append([points[nRay][0], points[nRay][1], points[nRay][2],
                                  Ray0['K'][nRay][0], Ray0['K'][nRay][1], Ray0['K'][nRay][2],
                                  Ray0['e0'][nRay][0], Ray0['e0'][nRay][1], Ray0['e0'][nRay][2],
                                  0, 0])     # Change these initial conditions

    return InitialConditions, Ray0, RayTube0

###################################################################

def Dipole(rho, a_location, E0_amplitude, nDivisionSource, iPol):

    InitialConditions = []
    e0 = []
    e0_spherical = []
    E0_sperical = []
    E0_rho = []

    points, cells = meshzoo.icosa_sphere(nDivisionSource)

    points = points*rho

    arrayNodeMap = np.zeros(shape=(cells.shape[0], 4), dtype='<i4')
    arrayNodeMap[:, 0] = 3
    arrayNodeMap[:, 1:4] = cells

    icoSphere = pv.PolyData(points, arrayNodeMap.ravel())

    icoSphere = icoSphere.translate(
        (a_location[0], a_location[1], a_location[2]), inplace=True)

    Normal = icoSphere.point_normals

    points = icoSphere.points

    cells = icoSphere.faces

    totRay = points.shape[0]

    for nRay in range(totRay):

        _, theta, phi = Functions.Cart2sphere(
            Normal[nRay, 0], Normal[nRay, 1], Normal[nRay, 2])

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        z_versor = np.array([0, 0, 1])

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        n0 = 1

        RayDir1 = r_versor[0]*n0
        RayDir2 = r_versor[1]*n0
        RayDir3 = r_versor[2]*n0

        if iPol == 1:
            e = np.array([np.cos(theta)*np.cos(phi),
                          np.cos(theta) * np.sin(phi),
                          -np.sin(theta)])
        else:

            e = np.array([-np.sin(phi),
                          np.cos(phi),
                          0])

        e_spherical = np.array([np.dot(e, r_versor),
                                np.dot(e, theta_versor),
                                np.dot(e, phi_versor)])

        e0.append(e)
        e0_spherical.append(e_spherical)

        Edir = e*np.sin(theta)

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        # E0_amplitude_rays=E0_amplitude*1/rho
        E0_amplitude_rays = E0_amplitude

        E0_rho.append(E0_amplitude_rays*Edir)

        E0_sperical.append(E0_amplitude_rays*Edir_spherical)

        InitialConditions.append([points[nRay][0], points[nRay][1], points[nRay][2], RayDir1,
                                  RayDir2, RayDir3, e[0], e[1], e[2], 0, 0])     # Change these initial conditions

    E0_rho = np.array(E0_rho)
    e0 = np.array(e0)
    e0_spherical = np.array(e0_spherical)
    E0_sperical = np.array(E0_sperical)

    Ray0 = pv.PolyData(np.array(InitialConditions)[:, 0:3])

    Ray0['kx'] = np.array(InitialConditions)[:, 3]
    Ray0['ky'] = np.array(InitialConditions)[:, 4]
    Ray0['kz'] = np.array(InitialConditions)[:, 5]

    Ray0['ex'] = np.array(InitialConditions)[:, 6]
    Ray0['ey'] = np.array(InitialConditions)[:, 7]
    Ray0['ez'] = np.array(InitialConditions)[:, 8]

    Ray0['er'] = e0_spherical[:, 0]
    Ray0['et'] = e0_spherical[:, 1]
    Ray0['ep'] = e0_spherical[:, 2]

    Ray0['K'] = np.column_stack([Ray0['kx'], Ray0['ky'], Ray0['kz']])

    Ray0['e0'] = np.column_stack([Ray0['ex'], Ray0['ey'], Ray0['ez']])
    Ray0['e0_spherical'] = np.column_stack(
        [Ray0['er'], Ray0['et'], Ray0['ep']])

    Ray0['E0'] = np.column_stack(
        [(E0_rho[:, 0]), (E0_rho[:, 1]), (E0_rho[:, 2])])

    Ray0['E0x'] = E0_rho[:, 0]
    Ray0['E0y'] = E0_rho[:, 1]
    Ray0['E0z'] = E0_rho[:, 2]

    Ray0['E0r'] = E0_sperical[:, 0]
    Ray0['E0t'] = E0_sperical[:, 1]
    Ray0['E0p'] = E0_sperical[:, 2]

    Ray0['|E0|'] = np.sqrt(Ray0['E0x']**2+Ray0['E0y']**2+Ray0['E0z']**2)

    Ray0['E0_spherical'] = np.column_stack(
        [Ray0['E0r'], Ray0['E0t'], Ray0['E0p']])

    RayTube0 = pv.PolyData(np.array(InitialConditions)[:, 0:3], cells)

    totTubes = cells.shape[0]//4

    e0_tube = []
    e0_tube_spherical = []
    E0_tube = []
    E0_tube_spherical = []
    ATube0 = []
    thetaTube = []
    phiTube = []

    zversor = np.array([0, 0, 1])
    rVersor_cell = RayTube0.face_normals
    RayTubes0 = RayTube0.compute_cell_sizes(
        length=False, area=True, volume=False)

    rho_centers = RayTube0.cell_centers()

    for nTubes in range(totTubes):

        ATube0.append(RayTubes0['Area'][nTubes])

        _, theta, phi = Functions.Cart2sphere(
            rVersor_cell[nTubes][0], rVersor_cell[nTubes][1], rVersor_cell[nTubes][2])

        thetaTube.append(theta)
        phiTube.append(phi)

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        e = -theta_versor
        e_spherical = np.array([np.dot(e, r_versor), np.dot(
            e, theta_versor), np.dot(e, phi_versor)])

        e0_tube.append(e)
        e0_tube_spherical.append(e_spherical)

        Edir = -theta_versor*np.sin(theta)

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        # E0_Tube = E0_amplitude*1/(np.linalg.norm(rho_centers.points[nTubes]))
        E0_Tube = E0_amplitude

        E0_tube.append(E0_Tube * Edir)

        E0_tube_spherical.append(E0_Tube*Edir_spherical)

    e0_tube = np.array(e0_tube)
    e0_tube_spherical = np.array(e0_tube_spherical)
    E0_tube = np.array(E0_tube)
    E0_tube_spherical = np.array(E0_tube_spherical)
    ATube0 = np.array(ATube0)
    thetaTube = np.array(thetaTube)
    phiTube = np.array(phiTube)

    RayTube0.cell_data['kx'] = rVersor_cell[:, 0]
    RayTube0.cell_data['ky'] = rVersor_cell[:, 1]
    RayTube0.cell_data['kz'] = rVersor_cell[:, 2]

    RayTube0.cell_data['ex'] = e0_tube[:, 0]
    RayTube0.cell_data['ey'] = e0_tube[:, 1]
    RayTube0.cell_data['ez'] = e0_tube[:, 2]

    RayTube0.cell_data['er'] = e0_tube_spherical[:, 0]
    RayTube0.cell_data['et'] = e0_tube_spherical[:, 1]
    RayTube0.cell_data['ep'] = e0_tube_spherical[:, 2]

    RayTube0.cell_data['E0x'] = E0_tube[:, 0]
    RayTube0.cell_data['E0y'] = E0_tube[:, 1]
    RayTube0.cell_data['E0z'] = E0_tube[:, 2]

    RayTube0.cell_data['E0r'] = E0_tube_spherical[:, 0]
    RayTube0.cell_data['E0t'] = E0_tube_spherical[:, 1]
    RayTube0.cell_data['E0p'] = E0_tube_spherical[:, 2]

    RayTube0.cell_data['K'] = np.column_stack(
        [RayTube0['kx'], RayTube0['ky'], RayTube0['kz']])

    RayTube0.cell_data['e0'] = np.column_stack(
        [RayTube0['ex'], RayTube0['ey'], RayTube0['ez']])
    RayTube0.cell_data['e0_spherical'] = np.column_stack(
        [RayTube0['er'], RayTube0['et'], RayTube0['ep']])

    RayTube0.cell_data['E0'] = np.column_stack(
        [RayTube0.cell_data['E0x'], RayTube0.cell_data['E0y'], RayTube0.cell_data['E0z']])

    RayTube0.cell_data['E0_spherical'] = np.column_stack(
        [RayTube0.cell_data['E0r'], RayTube0.cell_data['E0t'], RayTube0.cell_data['E0p']])

    RayTube0.cell_data['|E0|'] = np.sqrt(
        RayTube0['E0x']**2+RayTube0['E0y']**2+RayTube0['E0z']**2)

    RayTube0.cell_data['Area'] = ATube0

    return InitialConditions, Ray0, RayTube0

###################################################################


def PlaneWaveSquare(size, K, a_location, E0_amplitude, nDivisionSource, iPol):

    InitialConditions = []
    e0 = []
    e0_spherical = []
    E0_sperical = []
    E0_rho = []

    mesh = pv.Plane(i_size=size, j_size=size, direction=K,
                    i_resolution=nDivisionSource, j_resolution=nDivisionSource)

    mesh = mesh.triangulate()

    mesh = mesh.translate(
        (a_location[0], a_location[1], a_location[2]), inplace=True)

    points = mesh.points

    cells = mesh.faces

    Normal = np.array(K)

    totRay = points.shape[0]

    for nRay in range(totRay):

        _, theta, phi = Functions.Cart2sphere(
            Normal[0], Normal[1], Normal[2])

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        n0 = 1

        RayDir1 = r_versor[0]*n0
        RayDir2 = r_versor[1]*n0
        RayDir3 = r_versor[2]*n0

        if iPol == 1:
            e = np.array([np.cos(theta)*np.cos(phi),
                          np.cos(theta) * np.sin(phi),
                          -np.sin(theta)])
        else:

            e = np.array([-np.sin(phi),
                          np.cos(phi),
                          0])

        e_spherical = np.array([np.dot(e, r_versor),
                                np.dot(e, theta_versor),
                                np.dot(e, phi_versor)])

        e0.append(e)
        e0_spherical.append(e_spherical)

        Edir = e

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        E0_rho.append(E0_amplitude*Edir)
        # E0_rho.append(E0_amplitude*e0*np.exp(1j*k*rho))

        E0_sperical.append(E0_amplitude*Edir_spherical)

        InitialConditions.append([points[nRay][0], points[nRay][1], points[nRay][2], RayDir1,
                                  RayDir2, RayDir3, e[0], e[1], e[2], 0, 0])     # Change these initial conditions

    E0_rho = np.array(E0_rho)
    e0 = np.array(e0)
    e0_spherical = np.array(e0_spherical)
    E0_sperical = np.array(E0_sperical)

    Ray0 = pv.PolyData(np.array(InitialConditions)[:, 0:3])

    Ray0['kx'] = np.array(InitialConditions)[:, 3]
    Ray0['ky'] = np.array(InitialConditions)[:, 4]
    Ray0['kz'] = np.array(InitialConditions)[:, 5]

    Ray0['ex'] = np.array(InitialConditions)[:, 6]
    Ray0['ey'] = np.array(InitialConditions)[:, 7]
    Ray0['ez'] = np.array(InitialConditions)[:, 8]

    Ray0['er'] = e0_spherical[:, 0]
    Ray0['et'] = e0_spherical[:, 1]
    Ray0['ep'] = e0_spherical[:, 2]

    Ray0['K'] = np.column_stack([Ray0['kx'], Ray0['ky'], Ray0['kz']])

    Ray0['e0'] = np.column_stack([Ray0['ex'], Ray0['ey'], Ray0['ez']])
    Ray0['e0_spherical'] = np.column_stack(
        [Ray0['er'], Ray0['et'], Ray0['ep']])

    Ray0['E0'] = np.column_stack(
        [(E0_rho[:, 0]), (E0_rho[:, 1]), (E0_rho[:, 2])])

    Ray0['E0x'] = E0_rho[:, 0]
    Ray0['E0y'] = E0_rho[:, 1]
    Ray0['E0z'] = E0_rho[:, 2]

    Ray0['E0r'] = E0_sperical[:, 0]
    Ray0['E0t'] = E0_sperical[:, 1]
    Ray0['E0p'] = E0_sperical[:, 2]

    Ray0['|E0|'] = np.sqrt(Ray0['E0x']**2+Ray0['E0y']**2+Ray0['E0z']**2)

    Ray0['E0_spherical'] = np.column_stack(
        [Ray0['E0r'], Ray0['E0t'], Ray0['E0p']])

    RayTube0 = pv.PolyData(np.array(InitialConditions)[:, 0:3], cells)

    totTubes = cells.shape[0]//4

    e0_tube = []
    e0_tube_spherical = []
    E0_tube = []
    E0_tube_spherical = []
    ATube0 = []
    thetaTube = []
    phiTube = []

    zversor = np.array([0, 0, 1])
    rVersor_cell = -RayTube0.face_normals
    RayTubes0 = RayTube0.compute_cell_sizes(
        length=False, area=True, volume=False)

    for nTubes in range(totTubes):

        ATube0.append(RayTubes0['Area'][nTubes])

        _, theta, phi = Functions.Cart2sphere(
            Normal[0], Normal[1], Normal[2])

        thetaTube.append(theta)
        phiTube.append(phi)

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        if iPol == 1:
            e = np.array([np.cos(theta)*np.cos(phi),
                          np.cos(theta) * np.sin(phi),
                          -np.sin(theta)])
        else:

            e = np.array([-np.sin(phi),
                          np.cos(phi),
                          0])

        e_spherical = np.array([np.dot(e, r_versor), np.dot(
            e, theta_versor), np.dot(e, phi_versor)])

        e0_tube.append(e)
        e0_tube_spherical.append(e_spherical)

        Edir = e

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        E0_tube.append(E0_amplitude * Edir)

        E0_tube_spherical.append(E0_amplitude*Edir_spherical)

    e0_tube = np.array(e0_tube)
    e0_tube_spherical = np.array(e0_tube_spherical)
    E0_tube = np.array(E0_tube)
    E0_tube_spherical = np.array(E0_tube_spherical)
    ATube0 = np.array(ATube0)
    thetaTube = np.array(thetaTube)
    phiTube = np.array(phiTube)

    RayTube0.cell_data['kx'] = rVersor_cell[:, 0]
    RayTube0.cell_data['ky'] = rVersor_cell[:, 1]
    RayTube0.cell_data['kz'] = rVersor_cell[:, 2]

    RayTube0.cell_data['ex'] = e0_tube[:, 0]
    RayTube0.cell_data['ey'] = e0_tube[:, 1]
    RayTube0.cell_data['ez'] = e0_tube[:, 2]

    RayTube0.cell_data['er'] = e0_tube_spherical[:, 0]
    RayTube0.cell_data['et'] = e0_tube_spherical[:, 1]
    RayTube0.cell_data['ep'] = e0_tube_spherical[:, 2]

    RayTube0.cell_data['E0x'] = E0_tube[:, 0]
    RayTube0.cell_data['E0y'] = E0_tube[:, 1]
    RayTube0.cell_data['E0z'] = E0_tube[:, 2]

    RayTube0.cell_data['E0r'] = E0_tube_spherical[:, 0]
    RayTube0.cell_data['E0t'] = E0_tube_spherical[:, 1]
    RayTube0.cell_data['E0p'] = E0_tube_spherical[:, 2]

    RayTube0.cell_data['K'] = np.column_stack(
        [RayTube0['kx'], RayTube0['ky'], RayTube0['kz']])

    RayTube0.cell_data['e0'] = np.column_stack(
        [RayTube0['ex'], RayTube0['ey'], RayTube0['ez']])
    RayTube0.cell_data['e0_spherical'] = np.column_stack(
        [RayTube0['er'], RayTube0['et'], RayTube0['ep']])

    RayTube0.cell_data['E0'] = np.column_stack(
        [RayTube0.cell_data['E0x'], RayTube0.cell_data['E0y'], RayTube0.cell_data['E0z']])

    RayTube0.cell_data['E0_spherical'] = np.column_stack(
        [RayTube0.cell_data['E0r'], RayTube0.cell_data['E0t'], RayTube0.cell_data['E0p']])

    RayTube0.cell_data['|E0|'] = np.sqrt(
        RayTube0['E0x']**2+RayTube0['E0y']**2+RayTube0['E0z']**2)

    RayTube0.cell_data['Area'] = ATube0

    return InitialConditions, Ray0, RayTube0

###################################################################


def PlaneWaveDisk(Radius, K, a_location, E0_amplitude, nDivisionSource, iPol):

    InitialConditions = []
    e0 = []
    e0_spherical = []
    E0_sperical = []
    E0_rho = []

    points, cells = meshzoo.disk(5, nDivisionSource)

    points = points * Radius

    pointsArray = np.zeros(shape=(points.shape[0], 3))
    pointsArray[:, 0:2] = points
    pointsArray[:, 2] = 0

    arrayNodeMap = np.zeros(shape=(cells.shape[0], 4), dtype='<i4')
    arrayNodeMap[:, 0] = 3
    arrayNodeMap[:, 1:4] = cells

    Disk = pv.PolyData(pointsArray, arrayNodeMap.ravel())

    Disk = Disk.translate(
        (a_location[0], a_location[1], a_location[2]), inplace=True)

    points = Disk.points

    cells = Disk.faces

    Normal = np.array(K)

    totRay = points.shape[0]

    for nRay in range(totRay):

        _, theta, phi = Functions.Cart2sphere(
            Normal[0], Normal[1], Normal[2])

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        n0 = 1

        RayDir1 = r_versor[0]*n0
        RayDir2 = r_versor[1]*n0
        RayDir3 = r_versor[2]*n0

        if iPol == 1:
            e = np.array([np.cos(theta)*np.cos(phi),
                          np.cos(theta) * np.sin(phi),
                          -np.sin(theta)])
        else:

            e = np.array([-np.sin(phi),
                          np.cos(phi),
                          0])

        e_spherical = np.array([np.dot(e, r_versor),
                                np.dot(e, theta_versor),
                                np.dot(e, phi_versor)])

        e0.append(e)
        e0_spherical.append(e_spherical)

        Edir = e

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        E0_rho.append(E0_amplitude*Edir)
        # E0_rho.append(E0_amplitude*e0*np.exp(1j*k*rho))

        E0_sperical.append(E0_amplitude*Edir_spherical)

        InitialConditions.append([points[nRay][0], points[nRay][1], points[nRay][2], RayDir1,
                                  RayDir2, RayDir3, e[0], e[1], e[2], 0, 0])     # Change these initial conditions

    E0_rho = np.array(E0_rho)
    e0 = np.array(e0)
    e0_spherical = np.array(e0_spherical)
    E0_sperical = np.array(E0_sperical)

    Ray0 = pv.PolyData(np.array(InitialConditions)[:, 0:3])

    Ray0['kx'] = np.array(InitialConditions)[:, 3]
    Ray0['ky'] = np.array(InitialConditions)[:, 4]
    Ray0['kz'] = np.array(InitialConditions)[:, 5]

    Ray0['ex'] = np.array(InitialConditions)[:, 6]
    Ray0['ey'] = np.array(InitialConditions)[:, 7]
    Ray0['ez'] = np.array(InitialConditions)[:, 8]

    Ray0['er'] = e0_spherical[:, 0]
    Ray0['et'] = e0_spherical[:, 1]
    Ray0['ep'] = e0_spherical[:, 2]

    Ray0['K'] = np.column_stack([Ray0['kx'], Ray0['ky'], Ray0['kz']])

    Ray0['e0'] = np.column_stack([Ray0['ex'], Ray0['ey'], Ray0['ez']])
    Ray0['e0_spherical'] = np.column_stack(
        [Ray0['er'], Ray0['et'], Ray0['ep']])

    Ray0['E0'] = np.column_stack(
        [(E0_rho[:, 0]), (E0_rho[:, 1]), (E0_rho[:, 2])])

    Ray0['E0x'] = E0_rho[:, 0]
    Ray0['E0y'] = E0_rho[:, 1]
    Ray0['E0z'] = E0_rho[:, 2]

    Ray0['E0r'] = E0_sperical[:, 0]
    Ray0['E0t'] = E0_sperical[:, 1]
    Ray0['E0p'] = E0_sperical[:, 2]

    Ray0['|E0|'] = np.sqrt(Ray0['E0x']**2+Ray0['E0y']**2+Ray0['E0z']**2)

    Ray0['E0_spherical'] = np.column_stack(
        [Ray0['E0r'], Ray0['E0t'], Ray0['E0p']])

    RayTube0 = pv.PolyData(np.array(InitialConditions)[:, 0:3], cells)

    totTubes = cells.shape[0]//4

    e0_tube = []
    e0_tube_spherical = []
    E0_tube = []
    E0_tube_spherical = []
    ATube0 = []
    thetaTube = []
    phiTube = []

    zversor = np.array([0, 0, 1])
    rVersor_cell = -RayTube0.face_normals
    RayTubes0 = RayTube0.compute_cell_sizes(
        length=False, area=True, volume=False)

    for nTubes in range(totTubes):

        ATube0.append(RayTubes0['Area'][nTubes])

        _, theta, phi = Functions.Cart2sphere(
            Normal[0], Normal[1], Normal[2])

        thetaTube.append(theta)
        phiTube.append(phi)

        st = np.sin(theta)
        ct = np.cos(theta)
        cp = np.cos(phi)
        sp = np.sin(phi)

        r_versor = np.array([st*cp, st*sp, ct])
        theta_versor = np.array([ct*cp, ct*sp, -st])
        phi_versor = np.array([-sp, cp, 0])

        if iPol == 1:
            e = np.array([np.cos(theta)*np.cos(phi),
                          np.cos(theta) * np.sin(phi),
                          -np.sin(theta)])
        else:

            e = np.array([-np.sin(phi),
                          np.cos(phi),
                          0])

        e_spherical = np.array([np.dot(e, r_versor), np.dot(
            e, theta_versor), np.dot(e, phi_versor)])

        e0_tube.append(e)
        e0_tube_spherical.append(e_spherical)

        Edir = e

        Edir_spherical = np.array([np.dot(Edir, r_versor),
                                   np.dot(Edir, theta_versor),
                                   np.dot(Edir, phi_versor)])

        E0_tube.append(E0_amplitude * Edir)

        E0_tube_spherical.append(E0_amplitude*Edir_spherical)

    e0_tube = np.array(e0_tube)
    e0_tube_spherical = np.array(e0_tube_spherical)
    E0_tube = np.array(E0_tube)
    E0_tube_spherical = np.array(E0_tube_spherical)
    ATube0 = np.array(ATube0)
    thetaTube = np.array(thetaTube)
    phiTube = np.array(phiTube)

    RayTube0.cell_data['kx'] = rVersor_cell[:, 0]
    RayTube0.cell_data['ky'] = rVersor_cell[:, 1]
    RayTube0.cell_data['kz'] = rVersor_cell[:, 2]

    RayTube0.cell_data['ex'] = e0_tube[:, 0]
    RayTube0.cell_data['ey'] = e0_tube[:, 1]
    RayTube0.cell_data['ez'] = e0_tube[:, 2]

    RayTube0.cell_data['er'] = e0_tube_spherical[:, 0]
    RayTube0.cell_data['et'] = e0_tube_spherical[:, 1]
    RayTube0.cell_data['ep'] = e0_tube_spherical[:, 2]

    RayTube0.cell_data['E0x'] = E0_tube[:, 0]
    RayTube0.cell_data['E0y'] = E0_tube[:, 1]
    RayTube0.cell_data['E0z'] = E0_tube[:, 2]

    RayTube0.cell_data['E0r'] = E0_tube_spherical[:, 0]
    RayTube0.cell_data['E0t'] = E0_tube_spherical[:, 1]
    RayTube0.cell_data['E0p'] = E0_tube_spherical[:, 2]

    RayTube0.cell_data['K'] = np.column_stack(
        [RayTube0['kx'], RayTube0['ky'], RayTube0['kz']])

    RayTube0.cell_data['e0'] = np.column_stack(
        [RayTube0['ex'], RayTube0['ey'], RayTube0['ez']])
    RayTube0.cell_data['e0_spherical'] = np.column_stack(
        [RayTube0['er'], RayTube0['et'], RayTube0['ep']])

    RayTube0.cell_data['E0'] = np.column_stack(
        [RayTube0.cell_data['E0x'], RayTube0.cell_data['E0y'], RayTube0.cell_data['E0z']])

    RayTube0.cell_data['E0_spherical'] = np.column_stack(
        [RayTube0.cell_data['E0r'], RayTube0.cell_data['E0t'], RayTube0.cell_data['E0p']])

    RayTube0.cell_data['|E0|'] = np.sqrt(
        RayTube0['E0x']**2+RayTube0['E0y']**2+RayTube0['E0z']**2)

    RayTube0.cell_data['Area'] = ATube0

    return InitialConditions, Ray0, RayTube0
