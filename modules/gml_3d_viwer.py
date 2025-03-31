import xml.etree.ElementTree as ET
from os.path import basename, splitext
import os
import pyvista as pv
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from tqdm import tqdm
import pickle as pkl

from .map_image import draw_base_map
from .sun import get_sulight_vector
from .shadow import process_shadows, save_shadows_to_geojson


def process_buildings(root, ns):
    """
    Processes buildings from the XML tree and generates 3D meshes and 2D footprints.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        Root of the XML tree from the GML file.
    ns : dict
        Dictionary with the GML namespaces.

    Returns
    -------
    tuple
        buildings_3d : list
            List of 3D meshes of the buildings.
        footprints_polygons : list
            List of polygons representing the building footprints.
        all_coords_for_map : list
            List of all coordinates for the base map.
    """
    buildings_3d = []
    footprints_polygons = []
    all_coords_for_map = []

    print("Processing buildings...")

    for building in tqdm(root.findall(".//bu-ext2d:BuildingPart", ns), desc="Processing buildings"):
        # Get the number of floors above ground
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


def combine_and_save_processed_files(buildings_3d, buildings_3d_file,
                                     footprints_polygons, footprints_polygons_file,
                                     all_coords_for_map, all_coords_for_map_file):
    """
    Combines 3D meshes and saves them along with footprints and coordinates to files.

    Parameters
    ----------
    buildings_3d : list
        List of 3D meshes of the buildings.
    buildings_3d_file : str
        Path to the file where the combined mesh will be saved.
    footprints_polygons : list
        List of polygons representing the building footprints.
    footprints_polygons_file : str
        Path to the file where the footprints will be saved.
    all_coords_for_map : list
        List of all coordinates for the base map.
    all_coords_for_map_file : str
        Path to the file where the coordinates will be saved.
    """
    if not buildings_3d:
        print("No meshes to combine.")
        return

    print("Combining building meshes...")
    combined_mesh = buildings_3d[0].copy()
    for bld in tqdm(buildings_3d[1:], desc="Combining building meshes"):
        combined_mesh = combined_mesh.merge(bld)

    print(f"Saving combined mesh to {buildings_3d_file}...")
    combined_mesh.save(buildings_3d_file)

    print(f"Saving footprints to {footprints_polygons_file}...")
    with open(footprints_polygons_file, "wb") as f:
        pkl.dump(footprints_polygons, f)

    print(f"Saving coordinates to {all_coords_for_map_file}...")
    with open(all_coords_for_map_file, "wb") as f:
        pkl.dump(all_coords_for_map, f)


def load_or_process_buildings(gml_file_path, root, ns):
    """
    Loads the combined mesh from a file if it exists, or processes the buildings and saves the combined mesh.

    Parameters
    ----------
    gml_file_path : str
        Path to the GML file.
    root : xml.etree.ElementTree.Element
        Root of the XML tree from the GML file.
    ns : dict
        Dictionary with the GML namespaces.

    Returns
    -------
    pv.PolyData
        Combined mesh of the buildings.
    list
        List of polygons representing the building footprints.
    list
        List of all coordinates for the base map.
    """
    combined_mesh_file = f"./data/processed_files/{splitext(basename(gml_file_path))[0]}_combined_mesh.vtk"
    footprints_polygons_file = f"./data/processed_files/{splitext(basename(gml_file_path))[0]}_footprints_polygons.pkl"
    all_coords_for_map_file = f"./data/processed_files/{splitext(basename(gml_file_path))[0]}_all_coords_for_map.pkl"

    if os.path.exists(combined_mesh_file) and os.path.exists(footprints_polygons_file) and os.path.exists(all_coords_for_map_file):
        print(f"Loading existing combined mesh from {combined_mesh_file}...")
        combined_mesh = pv.read(combined_mesh_file)

        footprints_polygons = pkl.load(open(footprints_polygons_file, "rb"))
        all_coords_for_map = pkl.load(open(all_coords_for_map_file, "rb"))
        
        return combined_mesh, footprints_polygons, all_coords_for_map

    print("Processing buildings...")
    buildings_3d, footprints_polygons, all_coords_for_map = process_buildings(root, ns)
    combine_and_save_processed_files(buildings_3d, combined_mesh_file,
                                     footprints_polygons, footprints_polygons_file,
                                     all_coords_for_map, all_coords_for_map_file)
    combined_mesh = pv.read(combined_mesh_file)

    return combined_mesh, footprints_polygons, all_coords_for_map


def process_footprints(footprints_polygons):
    """
    Unifies building footprints into a MultiPolygon.

    Parameters
    ----------
    footprints_polygons : list
        List of polygons representing the building footprints.

    Returns
    -------
    MultiPolygon
        Unified building footprints.
    """
    if footprints_polygons:
        return unary_union(footprints_polygons)
    return MultiPolygon()


def add_base_plane(plotter, all_coords_for_map, texture_map):
    """
    Generates and adds the textured base plane to the renderer.

    Parameters
    ----------
    plotter : pv.Plotter
        PyVista renderer.
    all_coords_for_map : list
        List of coordinates for the base map.
    texture_map : bool
        Indicates whether to add a texture to the base plane.

    Returns
    -------
    tuple
        Central coordinates of the base plane.
    """
    if texture_map and all_coords_for_map:
        pad = 30
        base_plane, texture, center_xy = draw_base_map(all_coords_for_map, pad)
        plotter.add_mesh(base_plane, texture=texture)
        return center_xy
    return np.mean(all_coords_for_map, axis=0) if all_coords_for_map else (0, 0)


def calculate_and_add_shadows(plotter, combined_mesh, all_buildings_footprints, center_xy, dt, gml_file_path, remove_bases=True):
    """
    Calculates shadows and adds them to the renderer.

    Parameters
    ----------
    plotter : pv.Plotter
        PyVista renderer.
    combined_mesh : pv.PolyData
        Combined mesh of the buildings.
    all_buildings_footprints : MultiPolygon
        Unified building footprints.
    center_xy : tuple
        Central coordinates to calculate the sunlight direction.
    dt : pd.Timestamp
        Date and time to calculate the sun's position.
    remove_bases : bool, optional
        Indicates whether to remove the bases of the shadows. Default is True.
    """
    if dt:
        sunlight_direction = get_sulight_vector(center_xy[0], center_xy[1], dt, convert_coords=True)
        shadow_mesh = process_shadows(combined_mesh, sunlight_direction)
        
        save_shadows_to_geojson(
            shadow_mesh,
            f"./data/shadow_geojson/{splitext(basename(gml_file_path))[0]}_{sunlight_direction[0]}_{sunlight_direction[1]}_{sunlight_direction[2]}.geojson",
            all_buildings_footprints,
            remove_bases=remove_bases
        )
        
        plotter.add_mesh(shadow_mesh, color="gray", opacity=1, show_edges=False, label="Shadows")


def render_scene(plotter, combined_mesh, dt):
    """
    Configures and displays the scene in the renderer.

    Parameters
    ----------
    plotter : pv.Plotter
        PyVista renderer.
    combined_mesh : pv.PolyData
        Combined mesh of the buildings.
    dt : pd.Timestamp
        Date and time to display in the scene.
    """
    plotter.add_mesh(combined_mesh, color="lightblue", opacity=1, show_edges=False, label="Unified buildings")
    if dt:
        plotter.add_text(f"Date and Time: {dt}", position='upper_left', font_size=10, color='black')
    plotter.show_grid()
    plotter.view_isometric()
    plotter.show()


def gml_3d_from_file(gml_file_path, dt, texture_map=True):
    """
    Processes a GML file and generates a 3D visualization of the buildings.

    Parameters
    ----------
    gml_file_path : str
        Path to the GML file.
    dt : pd.Timestamp
        Date and time to calculate the sun's position.
    texture_map : bool, optional
        Indicates whether to add a texture to the base plane. Default is True.
    """
    try:
        # Load the GML file
        tree = ET.parse(gml_file_path)
        root = tree.getroot()
        ns = {
            "gml": "http://www.opengis.net/gml/3.2",
            "bu-ext2d": "http://inspire.jrc.ec.europa.eu/schemas/bu-ext2d/2.0"
        }

        # Load or process buildings
        combined_mesh, footprints_polygons, all_coords_for_map = load_or_process_buildings(gml_file_path, root, ns)

        # Process building footprints
        all_buildings_footprints = process_footprints(footprints_polygons)

        # Create the renderer
        plotter = pv.Plotter()

        # Add the base plane
        center_xy = add_base_plane(plotter, all_coords_for_map, texture_map)

        # Calculate and add shadows
        # Shadows coinciding with the building bases are not removed
        # But you can set the parameter remove_bases to True to remove them
        calculate_and_add_shadows(plotter, combined_mesh, all_buildings_footprints, center_xy, dt, gml_file_path, remove_bases=True)

        # Render the scene
        render_scene(plotter, combined_mesh, dt)

    except ET.ParseError as e:
        print(f"XML parsing error: {e}")
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")