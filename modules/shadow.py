import pyvista as pv
from shapely.geometry import Polygon, MultiPolygon, mapping
import numpy as np
import json
from shapely.ops import unary_union
from time import time

from .coordinates import convert_multipolygon_coordinates_EPSG_to_4326 
from tqdm import tqdm

def polydata_to_shapely(poly: pv.PolyData) -> MultiPolygon:
    """
    Converts pv.PolyData (2D at z=0) to Shapely MultiPolygon.
    Extracts polygons without merging. Optimized for performance.
    """
    if poly.n_cells == 0:
        return MultiPolygon()

    # Extracts faces and points
    faces = poly.faces.reshape((-1, 4))[:, 1:]
    points = poly.points[:, :2]

    # Creates polygons and filters invalid or zero-area ones
    polygons = [Polygon(points[face]) for face in faces]
    valid_polygons = [p for p in polygons if p.is_valid and p.area > 0]

    return MultiPolygon(valid_polygons) if valid_polygons else MultiPolygon()

def shapely_to_polydata(shp_geom) -> pv.PolyData:
    """
    Converts a Shapely geometry (Polygon or MultiPolygon) to a pv.PolyData.
    Assumes z=0 for the entire geometry.
    """
    if shp_geom.is_empty:
        return pv.PolyData()

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

    all_coords = np.vstack(all_coords)
    all_faces = np.array(all_faces, dtype=np.int32)

    pd = pv.PolyData(all_coords, all_faces)
    return pd.triangulate().clean()

def project_mesh_onto_z0(mesh: pv.PolyData, direction: np.ndarray) -> pv.PolyData:
    """
    Projects 'mesh' onto the z=0 plane using the given direction.
    Does not merge the resulting faces of the projection.

    Parameters
    ----------
    mesh : pv.PolyData
        Mesh to project.
    direction : np.ndarray
        Direction vector (dx, dy, dz) for the projection.

    Returns
    -------
    pv.PolyData
        New mesh projected onto z=0.
    """
    if direction[2] == 0:
        raise ValueError("The z component of the projection vector is zero. Cannot project onto z=0.")

    # Vectorized computation for projection
    # formula: shadow_points = original_points + t * direction
    # where t = -original_points[:, 2] / direction[2]
    # This computes the intersection of the line defined by the direction vector with the z=0 plane.
    # The direction vector is assumed to be normalized.
    original_points = mesh.points
    t = -original_points[:, 2] / direction[2]
    shadow_points = original_points + t[:, np.newaxis] * direction
    

    # Generates the new mesh with the same connectivity without merging the faces.
    shadow_mesh = pv.PolyData(shadow_points, mesh.faces)
    return shadow_mesh

def process_shadows(buildings_combined_mesh: pv.PolyData, sunlight_direction: np.ndarray) -> pv.PolyData:
    """
    Calculates the shadows of 'buildings_combined_mesh' by projecting it onto z=0.
    """
    # Projects the mesh onto z=0
    shadow_mesh = project_mesh_onto_z0(buildings_combined_mesh, sunlight_direction)

    return shadow_mesh

def save_shadows_to_geojson(shadow_mesh, file_path, all_buildings_bases, epsg_source, remove_bases):
    """
    Saves the projected shadows to a GeoJSON file.
    Merges overlapping polygons while preserving inner rings (holes).
    Converts coordinates from EPSG:***** to EPSG:4326.

    Arguments
    ---------
    shadow_mesh : pv.PolyData
        The projected shadows mesh.
    file_path : str
        The path to save the GeoJSON file.
    all_buildings_bases : pv.PolyData
        The base mesh of all buildings.
    epsg_source : str
        The EPSG code of the coordinate system of the input mesh.
    remove_bases : bool
        Whether to remove the base of the buildings from the shadows.
    """
    
    shadow_polygons = polydata_to_shapely(shadow_mesh)

    # Merges overlapping polygons while preserving inner rings (holes)
    merged = unary_union(shadow_polygons)

    if remove_bases and all_buildings_bases:
        merged = merged.difference(all_buildings_bases)

    # Sets the reference system to EPSG:4326 and converts to GeoJSON
    merged_4326 = convert_multipolygon_coordinates_EPSG_to_4326(merged, epsg_source)
    geojson_dict = mapping(merged_4326)

    with open(file_path, 'w') as f:
        json.dump(geojson_dict, f)

    return file_path
