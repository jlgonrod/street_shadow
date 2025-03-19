import pyvista as pv
from shapely.geometry import Polygon, MultiPolygon, mapping
import numpy as np
import json
from shapely.ops import unary_union
from shapely.validation import make_valid
import geopandas as gpd


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

def process_shadows(buildings_combined_mesh: pv.PolyData, sunlight_direction: np.ndarray, building_footprints: MultiPolygon) -> pv.PolyData:
    """
    Calcula las sombras de 'buildings_combined_mesh' proyectándola sobre z=0 y 
    le quita las bases de los edificios restando la unión de las huellas.
    """
    # Proyecta la malla sobre z=0
    shadow_mesh = project_mesh_onto_z0(buildings_combined_mesh, sunlight_direction)
    shadow_mesh = shadow_mesh.triangulate().clean()

    # Convierte la malla proyectada a Shapely y "limpia" su geometría
    shadow_polygons = polydata_to_shapely(shadow_mesh).buffer(0)

    # Une todas las huellas y aplica un buffer para ampliar las bases
    building_footprints = unary_union(building_footprints).buffer(0)
    
    # Realiza la diferencia y "limpia" el resultado
    shadow_polygons = make_valid(shadow_polygons)
    building_footprints = make_valid(building_footprints)
    shadows_without_bases = shadow_polygons.difference(building_footprints).buffer(0)
    
    return shapely_to_polydata(shadows_without_bases)

def save_shadows_to_geojson(shadow_mesh, file_path):
    """
    Guarda las sombras proyectadas en un archivo GeoJSON.
    Se convierten las coordenadas de EPSG:25829 a EPSG:4326.
    """
    shadow_polygons = polydata_to_shapely(shadow_mesh)
    shadow_polygons = convert_multipolygon_coordinates_25829_to_4326(shadow_polygons)
    geojson_dict = mapping(shadow_polygons)
    with open(file_path, 'w') as f:
        json.dump(geojson_dict, f)
    return file_path