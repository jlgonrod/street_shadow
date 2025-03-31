import pyvista as pv
from shapely.geometry import Polygon, MultiPolygon, mapping
import numpy as np
import json
from shapely.ops import unary_union

from .coordinates import convert_multipolygon_coordinates_25829_to_4326 

def polydata_to_shapely(poly: pv.PolyData) -> MultiPolygon:
    """
    Converts a pv.PolyData (2D at z=0) to a Shapely MultiPolygon.
    Extracts the polygons without merging them.
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

    # Returns the polygons without merging them.
    return MultiPolygon(polygons)

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

    original_points = mesh.points.copy()
    shadow_points = []

    for p in original_points:
        # Parametric equation: p + t * direction, and we solve for t such that p.z + t * dz = 0
        t = -p[2] / direction[2]
        p_proj = p + t * direction
        shadow_points.append(p_proj)

    shadow_points = np.array(shadow_points)
    # Generates the new mesh with the same connectivity without merging the faces.
    shadow_mesh = pv.PolyData(shadow_points, mesh.faces)
    return shadow_mesh

def process_shadows(buildings_combined_mesh: pv.PolyData, sunlight_direction: np.ndarray) -> pv.PolyData:
    """
    Calculates the shadows of 'buildings_combined_mesh' by projecting it onto z=0.
    """
    # Projects the mesh onto z=0
    shadow_mesh = project_mesh_onto_z0(buildings_combined_mesh, sunlight_direction)
    shadow_mesh = shadow_mesh.triangulate().clean()

    return shadow_mesh

def save_shadows_to_geojson(shadow_mesh, file_path, all_buildings_bases, remove_bases):
    """
    Saves the projected shadows to a GeoJSON file.
    Merges overlapping polygons while preserving inner rings (holes).
    Converts coordinates from EPSG:25829 to EPSG:4326.
    """
    shadow_polygons = polydata_to_shapely(shadow_mesh)

    # Merges overlapping polygons while preserving inner rings (holes)
    merged = unary_union(shadow_polygons)

    if remove_bases and all_buildings_bases:
        merged = merged.difference(all_buildings_bases)

    # Sets the reference system to EPSG:4326 and converts to GeoJSON
    merged_4326 = convert_multipolygon_coordinates_25829_to_4326(merged)
    geojson_dict = mapping(merged_4326)

    with open(file_path, 'w') as f:
        json.dump(geojson_dict, f)

    return file_path
