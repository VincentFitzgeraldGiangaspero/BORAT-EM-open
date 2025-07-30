import pytecio
import numpy as np


def export_ray(raySolutionFile, raySolution):

    dataset_title = 'raytracing 3D'
    var_names = ['x0', 'x1', 'x2', 'k_x', 'k_y', 'k_z', 'ex', 'ey', 'ez', 'phase', 'path', 'n']

    VALUELOCATION_NODECENTERED = 1
    FD_DOUBLE = 2

    var_data_types = [FD_DOUBLE]*len(var_names)
    value_locations = [VALUELOCATION_NODECENTERED]*len(var_names)

    total_zones = np.shape(raySolution)[2]

    nVars = np.shape(raySolution)[0]

    if len(var_names) != nVars:
        print('Error: more states to save')

    zone_title = 'Ray n '

    file = pytecio.open_file(raySolutionFile, dataset_title, var_names)

    for nzone in range(total_zones):

        zname = zone_title + str(nzone+1)

        total_points = np.shape(raySolution[:,:,nzone])[1]
        # total_points = np.count_nonzero((raySolution[0, :, nzone]))

        I = total_points
        J = 1
        K = 1

        zone = pytecio.create_ordered_zone(
            file, zname, (I, J, K), value_locations=value_locations, var_data_types=var_data_types)

        for var in range(nVars):
            pytecio.zone_write_double_values(
                file, zone, var+1, raySolution[var,:,nzone])  # double vals

    pytecio.close_file(file)

    return


def export_ray_animation(raySolutonFile, rayAnimationFile):

    print('export_ray_animation ----> to be implemented')

    data = pytecio.read(raySolutonFile)

    dataset_title = 'raytracing 3D - animation'
    var_names = ['x0', 'x1', 'x2', 'e_x', 'e_y', 'e_z', 'phase', 'path', 'n']

    file_format = 1   # 0 - .plt;  1 - .szplt
    file_type = 0  # 0 - grid & solution; 1 - grid only; 2 - solution only
    data_type = 2  # 1 - single; 2 - double; ...

    VALUELOCATION_NODECENTERED = 1
    FD_DOUBLE = 2
    # FD_FLOAT = 1
    var_data_types = [FD_DOUBLE]*len(var_names)
    value_locations = [VALUELOCATION_NODECENTERED]*len(var_names)

    total_rays = data.numZones
    nVars = data.numVars
    total_steps = data.zone_info[0]['IJK'][0]

    raySolution = np.zeros(shape=(nVars, total_steps, total_rays))

    for ray in range(data.numZones):
        for var in range(data.numVars):
            for s in range(total_steps):
                raySolution[var, s, ray] = data._read_zone_var(
                    ray+1, var+1)[s][0][0]

    file = pytecio.open_file(rayAnimationFile, dataset_title, var_names)

    nzone = 0
    for nrays in range(total_rays):
        for nsteps in range(total_steps):

            nzone = nzone+1
            zname = 'Ray n ' + str(nrays+1) + '-step ' + str(nsteps)

            I = nsteps+1
            J = 1
            K = 1

            zone = pytecio.create_ordered_zone(
                file, zname, (I, J, K), value_locations=value_locations, var_data_types=var_data_types)
            pytecio.zone_set_solution_time(
                file, zone, strand=1, solution_time=nsteps+1)

            for var in range(nVars):
                pytecio.zone_write_double_values(
                    file, zone, var+1, raySolution[var, 0:nsteps+1, nrays])  # double vals

    pytecio.close_file(file)

    return
