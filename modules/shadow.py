import pyvista as pv
from shapely.ops import unary_union
from shapely.geometry import Polygon, MultiPolygon
import numpy as np


def unify_shadow_mesh(shadow_mesh: pv.PolyData) -> pv.PolyData:
    """
    Combines all overlapping faces in shadow_mesh into a single polygon.
    Then triangulates and cleans it.

    Parameters
    ----------
    shadow_mesh : pv.PolyData
        Mesh with possibly overlapping faces.

    Returns
    -------
    pv.PolyData
        Clean, unified mesh without overlaps.
    """
    # Extract polygons from the shadow mesh
    faces = shadow_mesh.faces.reshape(-1, 4)[:, 1:]  # take vertex indices
    points = shadow_mesh.points[:, :2]  # only XY coordinates (flat at z=0)

    polygons = []
    for face in faces:
        polygon_coords = points[face]
        poly = Polygon(polygon_coords)
        if poly.is_valid and poly.area > 0:
            polygons.append(poly)

    # Union of all polygons to eliminate overlaps
    unified_poly = unary_union(polygons)

    # Generate a PolyData from the union
    if unified_poly.is_empty:
        # Case where there is no shadow
        return pv.PolyData()

    if unified_poly.geom_type == 'MultiPolygon':
        all_coords = []
        all_faces = []
        face_offset = 0
        for poly in unified_poly.geoms:
            x, y = poly.exterior.coords.xy
            coords_2d = np.column_stack((x, y))
            coords_3d = np.hstack((coords_2d, np.zeros((len(coords_2d), 1))))
            n_pts = len(coords_3d)
            faces_poly = [n_pts] + list(range(face_offset, face_offset + n_pts))
            all_coords.append(coords_3d)
            all_faces.extend(faces_poly)
            face_offset += n_pts
        all_coords = np.vstack(all_coords)
        all_faces = np.array(all_faces, dtype=np.int32)
    else:
        # Simple polygon case
        x, y = unified_poly.exterior.coords.xy
        coords_2d = np.column_stack((x, y))
        coords_3d = np.hstack((coords_2d, np.zeros((len(coords_2d), 1))))
        n_pts = len(coords_3d)
        all_coords = coords_3d
        all_faces = np.array([n_pts] + list(range(n_pts)), dtype=np.int32)

    # Build the resulting mesh
    unified_shadow_mesh = pv.PolyData(all_coords, all_faces)
    unified_shadow_mesh = unified_shadow_mesh.triangulate().clean()
    return unified_shadow_mesh

def polydata_to_shapely(poly: pv.PolyData) -> MultiPolygon:
    """
    Converts a pv.PolyData (2D at z=0) to a Shapely MultiPolygon.
    Takes the faces and points and unifies them.

    Parameters
    ----------
    poly : pv.PolyData
        2D mesh (at z=0) to be converted.

    Returns
    -------
    MultiPolygon
        Resulting polygon or multipolygon.
    """
    if poly.n_cells < 1:
        return MultiPolygon()

    faces = poly.faces.reshape(-1, 4)[:, 1:]
    points = poly.points[:, :2]

    polygons = []
    for face in faces:
        ring_coords = points[face]
        shp_poly = Polygon(ring_coords)
        if shp_poly.is_valid and shp_poly.area > 0:
            polygons.append(shp_poly)

    if not polygons:
        return MultiPolygon()

    return unary_union(polygons)  # This can return Polygon or MultiPolygon

def shapely_to_polydata(shp_geom) -> pv.PolyData:
    """
    Converts a Shapely polygon or multipolygon to a PyVista PolyData.
    Assumes z=0 for all geometry.

    Parameters
    ----------
    shp_geom : shapely.geometry.BaseGeometry
        Shapely polygon or multipolygon.

    Returns
    -------
    pv.PolyData
        Corresponding triangulated and clean mesh.
    """
    if shp_geom.is_empty:
        return pv.PolyData()

    # Ensure it is MultiPolygon
    if shp_geom.geom_type == 'Polygon':
        shp_geom = MultiPolygon([shp_geom])

    all_coords = []
    all_faces = []
    offset = 0

    for poly in shp_geom.geoms:
        exterior_x, exterior_y = poly.exterior.coords.xy
        coords_2d = np.column_stack((exterior_x, exterior_y))
        coords_3d = np.hstack([coords_2d, np.zeros((len(coords_2d), 1))])

        n_pts = len(coords_3d)
        face = [n_pts] + list(range(offset, offset + n_pts))
        offset += n_pts

        all_coords.append(coords_3d)
        all_faces.extend(face)

        # If the polygon has holes, they could be handled separately
        # by adding additional faces, if necessary.

    all_coords = np.vstack(all_coords)
    all_faces = np.array(all_faces, dtype=np.int32)

    pd = pv.PolyData(all_coords, all_faces)
    return pd.triangulate().clean()

def project_mesh_onto_z0(mesh: pv.PolyData, direction: np.ndarray) -> pv.PolyData:
    """
    Projects the 'mesh' onto the z=0 plane using the given 'direction'.

    Parameters
    ----------
    mesh : pv.PolyData
        Mesh to be projected (surface).
    direction : np.ndarray
        Direction vector (dx, dy, dz) for the projection.

    Returns
    -------
    shadow_mesh : pv.PolyData
        New mesh projected onto z=0 with the same connectivity.

    Raises
    ------
    ValueError
        If direction[2] == 0, as there is no intersection with z=0.
    """
    if direction[2] == 0:
        raise ValueError("The z component of the projection vector is zero. "
                         "Cannot project onto z=0.")

    original_points = mesh.points.copy()
    shadow_points = []

    for p in original_points:
        # Parametric equation: p + t * direction
        # We seek z=0 => p.z + t*dz = 0 => t = -p.z/dz
        t = -p[2] / direction[2]
        p_proj = p + t * direction
        shadow_points.append(p_proj)

    shadow_points = np.array(shadow_points)

    # Build the new mesh with the same connectivity
    shadow_mesh = pv.PolyData(shadow_points, mesh.faces)

    # Unify possible overlapping faces
    shadow_mesh = unify_shadow_mesh(shadow_mesh)

    return shadow_mesh

def process_shadows(buildings_combined_mesh: pv.PolyData, sunlight_direction: np.ndarray, building_footprints: MultiPolygon) -> pv.PolyData:
    """
    """

    # Project the entire mesh onto z=0
    shadow_mesh = project_mesh_onto_z0(buildings_combined_mesh, sunlight_direction)
    shadow_mesh = shadow_mesh.triangulate().clean()

    # Convert the shadow to Shapely
    shadow_polygons = polydata_to_shapely(shadow_mesh)
    # Subtract ALL building footprints
    final_shadow = shadow_polygons.difference(building_footprints)

    # Convert the resulting shadow to PolyData
    shadow_mesh_no_bases = shapely_to_polydata(final_shadow)

    return shadow_mesh_no_bases