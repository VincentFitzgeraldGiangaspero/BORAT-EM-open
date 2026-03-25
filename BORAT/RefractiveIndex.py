"""
RefractiveIndex.py: Module for computing refractive index and gradients from CFD mesh
"""

import numpy as np


def compute_refractive_index(pvMesh, Solver):
    """
    Compute refractive index and gradients from CFD mesh data.

    This function is called when precomputed refractive index data is not available.
    It routes to the appropriate model based on Solver.Model configuration.

    Args:
        pvMesh: PyVista mesh object containing CFD data
        Solver: ConfigBorat.Read object with solver configuration

    Returns:
        RefractiveIndexArray: numpy array with shape (num_components, num_points)
                             where num_components is 1 if no gradients, 4 if gradients computed
    """

    model = Solver.Model.lower()

    if model == "non_collisional":
        return compute_non_collisional(pvMesh, Solver)
    else:
        raise ValueError(f"Unknown refractive index model: {model}")


def compute_non_collisional(pvMesh, Solver):
    """
    Compute refractive index for non-collisional plasma model.

    This model applies when collisions can be neglected in the plasma.
    The refractive index is computed from electron density and frequency.

    Physics:
        For a non-collisional plasma:
        n = sqrt(1 - X_e)
        where X_e = 80.5 * N_e / f^2

        N_e is the electron number density (m^-3)
        f is the frequency (Hz)

    Args:
        pvMesh: PyVista mesh object containing CFD data
        Solver: ConfigBorat.Read object with solver configuration

    Returns:
        RefractiveIndexArray: numpy array with shape (num_components, num_points)
                             where num_components is 1 if no gradients, 4 if gradients computed
    """

    # Extract electron number density from mesh
    try:
        electron_number_density = np.array(pvMesh[Solver.ElectronNumberDensity_VarName])
    except KeyError:
        raise ValueError(
            f"Electron number density variable '{Solver.ElectronNumberDensity_VarName}' not found in mesh. "
            f"Available variables: {list(pvMesh.array_names)}"
        )

    # Use electron number density directly (assumed to be in m^-3)
    number_density = electron_number_density

    # Get frequency from solver configuration
    frequency = Solver.Frequency

    # Compute X_e parameter
    # X_e = 80.5 * N_e / f^2
    Xe = 80.5 * number_density / (frequency**2)

    # Compute refractive index with clipping
    # ri = sqrt(max(0, 1 - X_e))
    ri_squared = np.maximum(0, 1.0 - Xe)
    refractive_index = np.sqrt(ri_squared)

    # Build output array
    RefractiveIndexArray = np.array([refractive_index])

    if Solver.precomputedGrad:
        # Compute gradients numerically from the mesh
        grad_x, grad_y, grad_z = compute_gradients(pvMesh.points, refractive_index)
        RefractiveIndexArray = np.vstack(
            [
                RefractiveIndexArray,
                grad_x,
                grad_y,
                grad_z,
            ]
        )

    return RefractiveIndexArray


def compute_gradients(mesh_points, scalar_field):
    """
    Compute scalar field gradients using finite differences on unstructured mesh.

    This is a helper function for computing gradients when they're not precomputed.
    Uses least-squares gradient reconstruction for unstructured meshes via PyVista.

    Args:
        mesh_points: numpy array of mesh point coordinates (num_points, 3)
        scalar_field: numpy array of scalar values at mesh points (num_points,)

    Returns:
        grad_x, grad_y, grad_z: gradient components (each shape: num_points,)
    """

    try:
        import pyvista as pv
    except ImportError:
        raise ImportError("PyVista is required for gradient computation. Install it with: pip install pyvista")

    # Create a simple point cloud from mesh points
    cloud = pv.PolyData(mesh_points)
    cloud["scalar_field"] = scalar_field

    # Compute gradients using PyVista's gradient filter
    # This uses the least-squares gradient reconstruction method
    cloud_with_grad = cloud.compute_derivative(scalars="scalar_field")

    # Extract gradient components
    # PyVista returns gradients as a (num_points, 3) array
    gradients = cloud_with_grad["Gradient"]

    grad_x = gradients[:, 0]
    grad_y = gradients[:, 1]
    grad_z = gradients[:, 2]

    return grad_x, grad_y, grad_z
