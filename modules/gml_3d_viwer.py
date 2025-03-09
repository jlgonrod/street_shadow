import xml.etree.ElementTree as ET
import pyvista as pv
import numpy as np
from .map_image import draw_base_map
from .coordinates import get_square_coords_from_coords
import os

def gml_3d_from_file(gml_file_path, texture_map=True):
    """
    Processes a GML file and generates 3D models of buildings.

    Parameters
    ----------
    gml_file_path : str
        Path to the GML file.
    texture_map : bool, optional
        Whether to add a texture map to the base plane, by default True.
    """
    try:
        # Load the GML file
        tree = ET.parse(gml_file_path)
        root = tree.getroot()

        # Define namespaces
        ns = {
            "gml": "http://www.opengis.net/gml/3.2",
            "bu-ext2d": "http://inspire.jrc.ec.europa.eu/schemas/bu-ext2d/2.0"
        }

        # List to store 3D models of buildings
        buildings_3d = []
        # List to store building coordinates
        all_coords = []

        for building in root.findall(".//bu-ext2d:BuildingPart", ns):
            # Get the number of floors above ground
            floors_above_element = building.find(".//bu-ext2d:numberOfFloorsAboveGround", ns)
            num_floors_above = int(floors_above_element.text) if floors_above_element is not None else 0

            # If there are no floors above ground, ignore this part of the building
            if num_floors_above == 0:
                continue

            # Define height
            height_above = num_floors_above * 3  # 3m per floor
            base_z = 0  # Base at ground level
            total_height = height_above  

            # Get the base geometry (2D) and extrude it to 3D
            for poslist in building.findall(".//gml:posList", ns):
                if poslist.text:
                    coords_text = poslist.text.strip()
                    coords_2d = np.array([list(map(float, coords_text.split()[i:i+2])) for i in range(0, len(coords_text.split()), 2)])

                    # Store coordinates to calculate the bounding box
                    all_coords.extend(coords_2d)

                    # Ensure the polygon is closed
                    if not np.array_equal(coords_2d[0], coords_2d[-1]):
                        coords_2d = np.vstack([coords_2d, coords_2d[0]])
                    
                    # Convert 2D coordinates to 3D by adding the correct base Z
                    coords_3d = np.hstack([coords_2d, np.full((coords_2d.shape[0], 1), base_z*0.1)])
                    
                    # Create a closed mesh object
                    faces = [[len(coords_3d)] + list(range(len(coords_3d)))]
                    polydata = pv.PolyData(coords_3d, faces)
                    
                    # Triangulate the polygon to handle complex polygons
                    polydata = polydata.triangulate()
                    
                    # Extrude the polygon in height ensuring it is a closed solid
                    extruded = polydata.extrude((0, 0, total_height), capping=True)
                    
                    buildings_3d.append(extruded)

        # If models are generated, visualize them
        if buildings_3d:
            plotter = pv.Plotter()
            
            # Add the base map
            if texture_map:
                pad=30
                base_plane, texture = draw_base_map(all_coords, pad)
                plotter.add_mesh(base_plane, texture=texture)
            
            # Generate the 3d extrusion
            for mesh in buildings_3d:
                plotter.add_mesh(mesh, color="lightblue", opacity=1, show_edges=False)

            plotter.show_grid()
            plotter.view_isometric()
            plotter.show()
        else:
            print("No 3D models were generated.")

    except Exception as e:
        print(f"Error processing the GML: {e}")