import pyvista as pv
from shapely.geometry import Polygon, MultiPolygon, mapping
import numpy as np
import json
from shapely.ops import unary_union
from time import time
import os
from multiprocessing import Pool, cpu_count
from .coordinates import convert_multipolygon_coordinates_EPSG_to_4326
from shapely.validation import make_valid
import multiprocessing as mp

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

    polygons = [
        Polygon(points[face])
        for face in faces
        if len(face) >= 3
    ]

    return MultiPolygon(polygons) if polygons else MultiPolygon()

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

    # round the points to avoid floating point precision issues
    # rounded to 1cm (2 decimal places) because we are working 
    # with EPSG:25830 / EPSG:25829
    shadow_points = np.round(shadow_points, decimals=2)

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

def parallel_unary_union(multi_polygon, n_chunks=None):

    # Check if the input is a MultiPolygon    
    if multi_polygon.is_empty or len(multi_polygon.geoms) == 0:
        return MultiPolygon()    

    # Check the number of chunks
    if n_chunks is None:
        n_chunks = mp.cpu_count()

    # Extracts the polygons from the MultiPolygon
    polygons = list(multi_polygon.geoms)
    
    # Divides the polygons into chunks
    chunk_size = (len(polygons) + n_chunks - 1) // n_chunks
    chunks = [polygons[i:i + chunk_size] for i in range(0, len(polygons), chunk_size)]

    # Process each chunk in parallel
    with mp.Pool(processes=n_chunks) as pool:
        partial_results = pool.map(unary_union, chunks)

    # Filter out None and empty results
    partial_results = [g for g in partial_results if g is not None and not g.is_empty]
    if not partial_results:
        return None

    # Merge the results
    return unary_union(partial_results)

def save_shadows_to_geojson(shadow_mesh, file_path, all_buildings_bases, epsg_source, remove_bases, verbose=True):
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
    verbose : bool, optional
        Whether to print timing information. Defaults to True.
    """
    start_time = time()
    shadow_polygons = polydata_to_shapely(shadow_mesh)
    end_time = time()
    if verbose:
        print(f"\tConversion to Shapely took {end_time - start_time:.2f} seconds.")

    # Use parallel unary union
    start_time = time()
    merged = parallel_unary_union(shadow_polygons)
    end_time = time()
    if verbose:
        print(f"\tParallel unary union took {end_time - start_time:.2f} seconds.")

    if remove_bases and all_buildings_bases:
        start_time = time()
        merged = merged.difference(all_buildings_bases)
        end_time = time()
        if verbose:
            print(f"\tRemoving bases took {end_time - start_time:.2f} seconds.")

    # Sets the reference system to EPSG:4326 and converts to GeoJSON
    start_time = time()
    merged_4326 = convert_multipolygon_coordinates_EPSG_to_4326(merged, epsg_source)
    end_time = time()
    if verbose:
        print(f"\tCoordinate conversion took {end_time - start_time:.2f} seconds.")

    # Convert to GeoJSON format
    start_time = time()
    geojson_dict = mapping(merged_4326)
    end_time = time()
    if verbose:
        print(f"\tMapping to GeoJSON took {end_time - start_time:.2f} seconds.")

    # Save to file
    start_time = time()
    os.makedirs(os.path.dirname(file_path), exist_ok=True) # Ensure the directory exists
    with open(file_path, 'w') as f:
        json.dump(geojson_dict, f)
    end_time = time()
    if verbose:
        print(f"\tSaving GeoJSON took {end_time - start_time:.2f} seconds.")

    return file_path