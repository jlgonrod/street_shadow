import pyvista as pv
from shapely.geometry import Polygon, MultiPolygon, mapping
import numpy as np
import json
from shapely.ops import unary_union
from time import time
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from .coordinates import convert_multipolygon_coordinates_EPSG_to_4326
import cupy as cp

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

def parallel_unary_union(multi_polygon, num_processes=None):
    """
    Performs a parallel unary union on a MultiPolygon.

    Parameters
    ----------
    multi_polygon : MultiPolygon
        A Shapely MultiPolygon to merge.
    num_processes : int, optional
        Number of processes to use. Defaults to the number of CPU cores.

    Returns
    -------
    MultiPolygon
        The merged MultiPolygon.
    """
    if multi_polygon.is_empty or len(multi_polygon.geoms) == 0:
        return MultiPolygon()

    if num_processes is None:
        num_processes = cpu_count()

    # Extract individual polygons from the MultiPolygon
    polygons = list(multi_polygon.geoms)

    # Split polygons into chunks for parallel processing
    chunk_size = max(1, len(polygons) // num_processes)
    chunks = [polygons[i:i + chunk_size] for i in range(0, len(polygons), chunk_size)]

    # Perform unary_union on each chunk in parallel
    with Pool(processes=num_processes) as pool:
        partial_results = pool.map(unary_union, chunks)

    # Merge the partial results sequentially
    return unary_union(partial_results)

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
    start_time = time()
    shadow_polygons = polydata_to_shapely(shadow_mesh)
    end_time = time()
    print(f"\tConversion to Shapely took {end_time - start_time:.2f} seconds.")

    # Use parallel unary union
    start_time = time()
    merged = parallel_unary_union(shadow_polygons)
    end_time = time()
    print(f"\tParallel unary union took {end_time - start_time:.2f} seconds.")

    if remove_bases and all_buildings_bases:
        start_time = time()
        merged = merged.difference(all_buildings_bases)
        end_time = time()
        print(f"\tRemoving bases took {end_time - start_time:.2f} seconds.")

    # Sets the reference system to EPSG:4326 and converts to GeoJSON
    start_time = time()
    merged_4326 = convert_multipolygon_coordinates_EPSG_to_4326(merged, epsg_source)
    end_time = time()
    print(f"\tCoordinate conversion took {end_time - start_time:.2f} seconds.")

    start_time = time()
    geojson_dict = mapping(merged_4326)
    end_time = time()
    print(f"\tMapping to GeoJSON took {end_time - start_time:.2f} seconds.")

    start_time = time()
    with open(file_path, 'w') as f:
        json.dump(geojson_dict, f)
    end_time = time()
    print(f"\tSaving GeoJSON took {end_time - start_time:.2f} seconds.")

    return file_path
