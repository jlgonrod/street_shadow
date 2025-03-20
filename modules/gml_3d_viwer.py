import xml.etree.ElementTree as ET
from os.path import basename, splitext
import pyvista as pv
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from tqdm import tqdm

from .map_image import draw_base_map
from .sun import get_sulight_vector
from .shadow import process_shadows, save_shadows_to_geojson


def gml_3d_from_file(gml_file_path, dt, texture_map=True):
    """
    Processes a GML file and generates:
      - 3D geometries of buildings
      - Optionally shadows on z=0 (if dt is not None)
      - A base with texture (optional)

    Removes the shadow projection on the base of ALL buildings,
    so that no shadow is seen on any footprint.

    Parameters
    ----------
    gml_file_path : str
        Path to the GML file.
    dt : pd.Timestamp
        Date/time (with timezone) to calculate the sun position.
        If None, shadows are not projected.
    texture_map : bool, optional
        Indicates if a map texture is added to the base.
        Default is True.
    """
    try:
        # Load the GML
        tree = ET.parse(gml_file_path)
        root = tree.getroot()

        # Typical namespaces
        ns = {
            "gml": "http://www.opengis.net/gml/3.2",
            "bu-ext2d": "http://inspire.jrc.ec.europa.eu/schemas/bu-ext2d/2.0"
        }

        # Lists for 3D meshes and 2D footprints
        buildings_3d = []
        footprints_polygons = []  # to store footprints (Polygons)
        all_coords_for_map = []   # to draw the base texture

        # Process each building
        print("Processing buildings...")

        for building in tqdm(root.findall(".//bu-ext2d:BuildingPart", ns), desc="Processing buildings"):
            # Get number of floors above ground
            floors_above_element = building.find(".//bu-ext2d:numberOfFloorsAboveGround", ns)
            num_floors_above = int(floors_above_element.text) if floors_above_element is not None else 0

            if num_floors_above == 0:
                # If there are no floors above ground, ignore this building
                continue

            # Total height (e.g., 3 m per floor; adjust according to your real model)
            total_height = num_floors_above * 3
            base_z = 0

            # Look for posList (base 2D coordinates of the polygon)
            for poslist in building.findall(".//gml:posList", ns):
                if poslist.text:
                    coords_text = poslist.text.strip()
                    # Convert to array of floats [x,y], in pairs
                    splitted = coords_text.split()
                    coords_2d = np.array([
                        list(map(float, splitted[i:i+2]))
                        for i in range(0, len(splitted), 2)
                    ])

                    # Accumulate for the global footprint (texture bounding box)
                    all_coords_for_map.extend(coords_2d)

                    # Ensure the polygon is closed
                    if not np.array_equal(coords_2d[0], coords_2d[-1]):
                        coords_2d = np.vstack([coords_2d, coords_2d[0]])

                    # Save the footprint in Shapely (to unify later)
                    # Note: there are cases with holes, etc. Here simple:
                    shp_polygon = Polygon(coords_2d)
                    # Ensure it is not empty or invalid
                    if shp_polygon.is_valid and shp_polygon.area > 0:
                        footprints_polygons.append(shp_polygon)

                    # Generate the 3D base in PyVista
                    coords_3d = np.hstack([
                        coords_2d,
                        np.full((coords_2d.shape[0], 1), base_z)
                    ])
                    # Create a base polygon in PyVista
                    faces = [[len(coords_3d)] + list(range(len(coords_3d)))]
                    base_polydata = pv.PolyData(coords_3d, faces)
                    base_polydata = base_polydata.triangulate()

                    # Extrude to create the wall + top cap
                    extruded = base_polydata.extrude((0, 0, total_height), capping=True)
                    # Extract the surface
                    extruded = extruded.extract_surface().triangulate()
                    buildings_3d.append(extruded)

        # If there are no buildings, exit
        if not buildings_3d:
            print("No 3D model generated (no building parts or num_floors_aboveGround=0).")
            return

        print("Combining building meshes...")
        # 1) Unify all meshes into one (for visualization)
        combined_mesh = buildings_3d[0].copy()
        for bld in tqdm(buildings_3d[1:], desc="Combining building meshes"):
            combined_mesh = combined_mesh.merge(bld)

        print("Processing building footprints...")
        # 2) Unify all footprints into a MultiPolygon
        #    to then subtract the shadows
        if footprints_polygons:
            all_buildings_footprints = unary_union(footprints_polygons)
        else:
            all_buildings_footprints = MultiPolygon()

        # Create the plotter
        plotter = pv.Plotter()

        # 3) Add the base plane with texture (optional)
        if texture_map and all_coords_for_map:
            pad = 30
            base_plane, texture, center_xy = draw_base_map(all_coords_for_map, pad)
            # Add the mesh with texture
            plotter.add_mesh(base_plane, texture=texture)
        else:
            # If draw_base_map is not used, we can at least
            # estimate the center for the solar direction:
            center_xy = np.mean(all_coords_for_map, axis=0) if all_coords_for_map else (0, 0)

        print("Adding meshes to the plot...")
        # 4) Add the building mesh
        plotter.add_mesh(
            combined_mesh,
            color="lightblue", opacity=1, show_edges=False,
            label="Unified buildings"
        )

        print("Calculating shadows...")
        # 5) If dt is not None, calculate and project shadows
        if dt:
            # Sunlight direction vector
            sunlight_direction = get_sulight_vector(center_xy[0], center_xy[1], dt)

            # Calculate the shadows
            shadow_mesh_no_bases = process_shadows(combined_mesh, sunlight_direction, all_buildings_footprints)

            # Save the shadows without the bases
            save_shadows_to_geojson(shadow_mesh_no_bases, f"./data/temp/{splitext(basename(gml_file_path))[0]}_{sunlight_direction[0]}_{sunlight_direction[1]}_{sunlight_direction[2]}.geojson")

            # Add the final shadow
            plotter.add_mesh(
                shadow_mesh_no_bases,
                color="gray", opacity=0.8, show_edges=False,
                label="Shadow without bases"
            )

            # Add text with date/time if needed
            plotter.add_text(f"Date and Time: {dt}", position='upper_left', font_size=10, color='black')

        print("Showing the plot...")
        # Optional: grid and view
        plotter.show_grid()
        plotter.view_isometric()
        plotter.show()

    except ET.ParseError as e:
        print(f"XML parsing error: {e}")
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
