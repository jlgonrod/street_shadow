import xml.etree.ElementTree as ET
from os.path import basename, splitext
import os
import pyvista as pv
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from tqdm import tqdm

from .map_image import draw_base_map
from .sun import get_sulight_vector
from .shadow import process_shadows, save_shadows_to_geojson


def process_buildings(root, ns):
    """
    Procesa los edificios desde el árbol XML y genera las mallas 3D y las huellas 2D.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        Raíz del árbol XML del archivo GML.
    ns : dict
        Diccionario con los espacios de nombres (namespaces) del GML.

    Returns
    -------
    tuple
        buildings_3d : list
            Lista de mallas 3D de los edificios.
        footprints_polygons : list
            Lista de polígonos de las huellas de los edificios.
        all_coords_for_map : list
            Lista de todas las coordenadas para el mapa base.
    """
    buildings_3d = []
    footprints_polygons = []
    all_coords_for_map = []

    print("Processing buildings...")

    for building in tqdm(root.findall(".//bu-ext2d:BuildingPart", ns), desc="Processing buildings"):
        # Obtener el número de pisos sobre el suelo
        floors_above_element = building.find(".//bu-ext2d:numberOfFloorsAboveGround", ns)
        num_floors_above = int(floors_above_element.text) if floors_above_element is not None else 0

        if num_floors_above == 0:
            continue

        total_height = num_floors_above * 3
        base_z = 0

        for poslist in building.findall(".//gml:posList", ns):
            if poslist.text:
                coords_text = poslist.text.strip()
                splitted = coords_text.split()
                coords_2d = np.array([
                    list(map(float, splitted[i:i+2]))
                    for i in range(0, len(splitted), 2)
                ])

                all_coords_for_map.extend(coords_2d)

                if not np.array_equal(coords_2d[0], coords_2d[-1]):
                    coords_2d = np.vstack([coords_2d, coords_2d[0]])

                shp_polygon = Polygon(coords_2d)
                if shp_polygon.is_valid and shp_polygon.area > 0:
                    footprints_polygons.append(shp_polygon)

                coords_3d = np.hstack([
                    coords_2d,
                    np.full((coords_2d.shape[0], 1), base_z)
                ])
                faces = [[len(coords_3d)] + list(range(len(coords_3d)))]
                base_polydata = pv.PolyData(coords_3d, faces)
                base_polydata = base_polydata.triangulate()

                extruded = base_polydata.extrude((0, 0, total_height), capping=True)
                extruded = extruded.extract_surface().triangulate()
                buildings_3d.append(extruded)

    return buildings_3d, footprints_polygons, all_coords_for_map


def combine_and_save_meshes(buildings_3d, output_file):
    """
    Combina las mallas 3D y las guarda en un archivo.

    Parameters
    ----------
    buildings_3d : list
        Lista de mallas 3D de los edificios.
    output_file : str
        Ruta del archivo donde se guardará la malla combinada.
    """
    if not buildings_3d:
        print("No hay mallas para combinar.")
        return

    print("Combining building meshes...")
    combined_mesh = buildings_3d[0].copy()
    for bld in tqdm(buildings_3d[1:], desc="Combining building meshes"):
        combined_mesh = combined_mesh.merge(bld)

    print(f"Saving combined mesh to {output_file}...")
    combined_mesh.save(output_file)


def load_or_process_buildings(gml_file_path, root, ns):
    """
    Carga la malla combinada desde un archivo si existe, o procesa los edificios y guarda la malla combinada.

    Parameters
    ----------
    gml_file_path : str
        Ruta del archivo GML.
    root : xml.etree.ElementTree.Element
        Raíz del árbol XML del archivo GML.
    ns : dict
        Diccionario con los espacios de nombres (namespaces) del GML.

    Returns
    -------
    pv.PolyData
        Malla combinada de los edificios.
    list
        Lista de polígonos de las huellas de los edificios.
    list
        Lista de todas las coordenadas para el mapa base.
    """
    combined_mesh_file = f"./data/vtk/{splitext(basename(gml_file_path))[0]}_combined_mesh.vtk"

    if os.path.exists(combined_mesh_file):
        print(f"Loading existing combined mesh from {combined_mesh_file}...")
        combined_mesh = pv.read(combined_mesh_file)
        return combined_mesh, None, None

    print("Processing buildings...")
    buildings_3d, footprints_polygons, all_coords_for_map = process_buildings(root, ns)
    combine_and_save_meshes(buildings_3d, combined_mesh_file)
    combined_mesh = pv.read(combined_mesh_file)

    return combined_mesh, footprints_polygons, all_coords_for_map


def process_footprints(footprints_polygons):
    """
    Unifica las huellas de los edificios en un MultiPolygon.

    Parameters
    ----------
    footprints_polygons : list
        Lista de polígonos de las huellas de los edificios.

    Returns
    -------
    MultiPolygon
        Huellas unificadas de los edificios.
    """
    if footprints_polygons:
        return unary_union(footprints_polygons)
    return MultiPolygon()


def add_base_plane(plotter, all_coords_for_map, texture_map):
    """
    Genera y agrega el plano base con textura al renderizador.

    Parameters
    ----------
    plotter : pv.Plotter
        Renderizador de PyVista.
    all_coords_for_map : list
        Lista de coordenadas para el mapa base.
    texture_map : bool
        Indica si se debe agregar una textura al plano base.

    Returns
    -------
    tuple
        Coordenadas centrales del plano base.
    """
    if texture_map and all_coords_for_map:
        pad = 30
        base_plane, texture, center_xy = draw_base_map(all_coords_for_map, pad)
        plotter.add_mesh(base_plane, texture=texture)
        return center_xy
    return np.mean(all_coords_for_map, axis=0) if all_coords_for_map else (0, 0)


def calculate_and_add_shadows(plotter, combined_mesh, all_buildings_footprints, center_xy, dt, gml_file_path):
    """
    Calcula las sombras y las agrega al renderizador.

    Parameters
    ----------
    plotter : pv.Plotter
        Renderizador de PyVista.
    combined_mesh : pv.PolyData
        Malla combinada de los edificios.
    all_buildings_footprints : MultiPolygon
        Huellas unificadas de los edificios.
    center_xy : tuple
        Coordenadas centrales para calcular la dirección de la luz solar.
    dt : pd.Timestamp
        Fecha y hora para calcular la posición del sol.
    """
    if dt:
        sunlight_direction = get_sulight_vector(center_xy[0], center_xy[1], dt)
        shadow_mesh_no_bases = process_shadows(combined_mesh, sunlight_direction, all_buildings_footprints)
        save_shadows_to_geojson(
            shadow_mesh_no_bases, 
            f"./data/shadow_geojson/{splitext(basename(gml_file_path))[0]}_{sunlight_direction[0]}_{sunlight_direction[1]}_{sunlight_direction[2]}.geojson"
        )
        plotter.add_mesh(shadow_mesh_no_bases, color="gray", opacity=0.8, show_edges=False, label="Shadows")


def render_scene(plotter, combined_mesh, dt):
    """
    Configura y muestra la escena en el renderizador.

    Parameters
    ----------
    plotter : pv.Plotter
        Renderizador de PyVista.
    combined_mesh : pv.PolyData
        Malla combinada de los edificios.
    dt : pd.Timestamp
        Fecha y hora para mostrar en la escena.
    """
    plotter.add_mesh(combined_mesh, color="lightblue", opacity=1, show_edges=False, label="Unified buildings")
    if dt:
        plotter.add_text(f"Date and Time: {dt}", position='upper_left', font_size=10, color='black')
    plotter.show_grid()
    plotter.view_isometric()
    plotter.show()


def gml_3d_from_file(gml_file_path, dt, texture_map=True):
    """
    Procesa un archivo GML y genera la visualización 3D de los edificios.

    Parameters
    ----------
    gml_file_path : str
        Ruta del archivo GML.
    dt : pd.Timestamp
        Fecha y hora para calcular la posición del sol.
    texture_map : bool, optional
        Indica si se debe agregar una textura al plano base. Por defecto es True.
    """
    try:
        # Cargar el archivo GML
        tree = ET.parse(gml_file_path)
        root = tree.getroot()
        ns = {
            "gml": "http://www.opengis.net/gml/3.2",
            "bu-ext2d": "http://inspire.jrc.ec.europa.eu/schemas/bu-ext2d/2.0"
        }

        # Cargar o procesar edificios
        combined_mesh, footprints_polygons, all_coords_for_map = load_or_process_buildings(gml_file_path, root, ns)

        # Procesar huellas de edificios
        all_buildings_footprints = process_footprints(footprints_polygons)

        # Crear el renderizador
        plotter = pv.Plotter()

        # Agregar el plano base
        center_xy = add_base_plane(plotter, all_coords_for_map, texture_map)

        # Calcular y agregar sombras
        calculate_and_add_shadows(plotter, combined_mesh, all_buildings_footprints, center_xy, dt, gml_file_path)

        # Renderizar la escena
        render_scene(plotter, combined_mesh, dt)

    except ET.ParseError as e:
        print(f"XML parsing error: {e}")
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")