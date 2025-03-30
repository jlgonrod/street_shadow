import pyvista as pv
from shapely.geometry import Polygon, MultiPolygon, mapping
import numpy as np
import json
from shapely.ops import unary_union

from .coordinates import convert_multipolygon_coordinates_25829_to_4326 

def polydata_to_shapely(poly: pv.PolyData) -> MultiPolygon:
    """
    Convierte un pv.PolyData (2D en z=0) a un Shapely MultiPolygon.
    Se extraen los polígonos sin unirlos.
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

    # Se retornan los polígonos sin unir.
    return MultiPolygon(polygons)

def shapely_to_polydata(shp_geom) -> pv.PolyData:
    """
    Convierte una geometría Shapely (Polygon o MultiPolygon) a un pv.PolyData.
    Se asume z=0 para toda la geometría.
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
    Proyecta 'mesh' sobre el plano z=0 usando la dirección dada.
    No se unifican las caras resultantes de la proyección.

    Parameters
    ----------

    mesh : pv.PolyData
        Malla a proyectar.
    direction : np.ndarray
        Vector de dirección (dx, dy, dz) para la proyección.

    Returns
    -------
    pv.PolyData
        Nueva malla proyectada sobre z=0.
    """
    if direction[2] == 0:
        raise ValueError("El componente z del vector de proyección es cero. No se puede proyectar sobre z=0.")

    original_points = mesh.points.copy()
    shadow_points = []

    for p in original_points:
        # Ecuación paramétrica: p + t * direction y buscamos t tal que p.z+t*dz = 0
        t = -p[2] / direction[2]
        p_proj = p + t * direction
        shadow_points.append(p_proj)

    shadow_points = np.array(shadow_points)
    # Se genera la nueva malla con la misma conectividad sin unificar las caras.
    shadow_mesh = pv.PolyData(shadow_points, mesh.faces)
    return shadow_mesh

def process_shadows(buildings_combined_mesh: pv.PolyData, sunlight_direction: np.ndarray) -> pv.PolyData:
    """
    Calcula las sombras de 'buildings_combined_mesh' proyectándola sobre z=0.
    """
    # Proyecta la malla sobre z=0
    shadow_mesh = project_mesh_onto_z0(buildings_combined_mesh, sunlight_direction)
    shadow_mesh = shadow_mesh.triangulate().clean()

    return shadow_mesh

def save_shadows_to_geojson(shadow_mesh, file_path, all_buildings_bases, remove_bases):
    """
    Guarda las sombras proyectadas en un archivo GeoJSON.
    Une los polígonos solapados preservando los anillos internos (huecos).
    Convierte las coordenadas de EPSG:25829 a EPSG:4326.
    """
    shadow_polygons = polydata_to_shapely(shadow_mesh)

    # Se unen los polígonos solapados preservando los anillos internos (huecos)
    merged = unary_union(shadow_polygons)

    if remove_bases and all_buildings_bases:
        merged = merged.difference(all_buildings_bases)

    # Se establece el sistema de referencia EPSG:4326 y se convierte a GeoJSON
    merged_4326 = convert_multipolygon_coordinates_25829_to_4326(merged)
    geojson_dict = mapping(merged_4326)

    with open(file_path, 'w') as f:
        json.dump(geojson_dict, f)

    return file_path