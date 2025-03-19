import xml.etree.ElementTree as ET
import pyvista as pv
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union

from .map_image import draw_base_map
from .sun import get_sulight_vector
from tqdm import tqdm


def unify_shadow_mesh(shadow_mesh: pv.PolyData) -> pv.PolyData:
    """
    Combines all overlapping faces in shadow_mesh into a single polygon.
    Then triangulates and cleans it.

    Parameters
    ----------
    shadow_mesh : pv.PolyData
        Mesh with possibly overlapping faces.

    Returns
    -------
    pv.PolyData
        Clean, unified mesh without overlaps.
    """
    # Extract polygons from the shadow mesh
    faces = shadow_mesh.faces.reshape(-1, 4)[:, 1:]  # take vertex indices
    points = shadow_mesh.points[:, :2]  # only XY coordinates (flat at z=0)

    polygons = []
    for face in faces:
        polygon_coords = points[face]
        poly = Polygon(polygon_coords)
        if poly.is_valid and poly.area > 0:
            polygons.append(poly)

    # Union of all polygons to eliminate overlaps
    unified_poly = unary_union(polygons)

    # Generate a PolyData from the union
    if unified_poly.is_empty:
        # Case where there is no shadow
        return pv.PolyData()

    if unified_poly.geom_type == 'MultiPolygon':
        all_coords = []
        all_faces = []
        face_offset = 0
        for poly in unified_poly.geoms:
            x, y = poly.exterior.coords.xy
            coords_2d = np.column_stack((x, y))
            coords_3d = np.hstack((coords_2d, np.zeros((len(coords_2d), 1))))
            n_pts = len(coords_3d)
            faces_poly = [n_pts] + list(range(face_offset, face_offset + n_pts))
            all_coords.append(coords_3d)
            all_faces.extend(faces_poly)
            face_offset += n_pts
        all_coords = np.vstack(all_coords)
        all_faces = np.array(all_faces, dtype=np.int32)
    else:
        # Simple polygon case
        x, y = unified_poly.exterior.coords.xy
        coords_2d = np.column_stack((x, y))
        coords_3d = np.hstack((coords_2d, np.zeros((len(coords_2d), 1))))
        n_pts = len(coords_3d)
        all_coords = coords_3d
        all_faces = np.array([n_pts] + list(range(n_pts)), dtype=np.int32)

    # Build the resulting mesh
    unified_shadow_mesh = pv.PolyData(all_coords, all_faces)
    unified_shadow_mesh = unified_shadow_mesh.triangulate().clean()
    return unified_shadow_mesh


def polydata_to_shapely(poly: pv.PolyData) -> MultiPolygon:
    """
    Converts a pv.PolyData (2D at z=0) to a Shapely MultiPolygon.
    Takes the faces and points and unifies them.

    Parameters
    ----------
    poly : pv.PolyData
        2D mesh (at z=0) to be converted.

    Returns
    -------
    MultiPolygon
        Resulting polygon or multipolygon.
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

    return unary_union(polygons)  # This can return Polygon or MultiPolygon


def shapely_to_polydata(shp_geom) -> pv.PolyData:
    """
    Converts a Shapely polygon or multipolygon to a PyVista PolyData.
    Assumes z=0 for all geometry.

    Parameters
    ----------
    shp_geom : shapely.geometry.BaseGeometry
        Shapely polygon or multipolygon.

    Returns
    -------
    pv.PolyData
        Corresponding triangulated and clean mesh.
    """
    if shp_geom.is_empty:
        return pv.PolyData()

    # Ensure it is MultiPolygon
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

        # If the polygon has holes, they could be handled separately
        # by adding additional faces, if necessary.

    all_coords = np.vstack(all_coords)
    all_faces = np.array(all_faces, dtype=np.int32)

    pd = pv.PolyData(all_coords, all_faces)
    return pd.triangulate().clean()


def project_mesh_onto_z0(mesh: pv.PolyData, direction: np.ndarray) -> pv.PolyData:
    """
    Projects the 'mesh' onto the z=0 plane using the given 'direction'.

    Parameters
    ----------
    mesh : pv.PolyData
        Mesh to be projected (surface).
    direction : np.ndarray
        Direction vector (dx, dy, dz) for the projection.

    Returns
    -------
    shadow_mesh : pv.PolyData
        New mesh projected onto z=0 with the same connectivity.

    Raises
    ------
    ValueError
        If direction[2] == 0, as there is no intersection with z=0.
    """
    if direction[2] == 0:
        raise ValueError("The z component of the projection vector is zero. "
                         "Cannot project onto z=0.")

    original_points = mesh.points.copy()
    shadow_points = []

    for p in original_points:
        # Parametric equation: p + t * direction
        # We seek z=0 => p.z + t*dz = 0 => t = -p.z/dz
        t = -p[2] / direction[2]
        p_proj = p + t * direction
        shadow_points.append(p_proj)

    shadow_points = np.array(shadow_points)

    # Build the new mesh with the same connectivity
    shadow_mesh = pv.PolyData(shadow_points, mesh.faces)

    # Unify possible overlapping faces
    shadow_mesh = unify_shadow_mesh(shadow_mesh)

    return shadow_mesh


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

        # 1) Unify all meshes into one (for visualization)
        combined_mesh = buildings_3d[0].copy()
        for bld in buildings_3d[1:]:
            combined_mesh = combined_mesh.merge(bld)

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

        # 4) Add the building mesh
        plotter.add_mesh(
            combined_mesh,
            color="lightblue", opacity=1, show_edges=False,
            label="Unified buildings"
        )

        # 5) If dt is not None, calculate and project shadows
        if dt:
            # Sunlight direction vector
            sunlight_direction = get_sulight_vector(center_xy[0], center_xy[1], dt)

            # Project the entire mesh onto z=0
            shadow_mesh = project_mesh_onto_z0(combined_mesh, sunlight_direction)
            shadow_mesh = shadow_mesh.triangulate().clean()

            # Convert the shadow to Shapely
            shadow_polygons = polydata_to_shapely(shadow_mesh)
            # Subtract ALL building footprints
            final_shadow = shadow_polygons.difference(all_buildings_footprints)

            # Convert the resulting shadow to PolyData
            shadow_mesh_no_bases = shapely_to_polydata(final_shadow)

            # Add the final shadow
            plotter.add_mesh(
                shadow_mesh_no_bases,
                color="gray", opacity=0.8, show_edges=False,
                label="Shadow without bases"
            )

            # Add text with date/time if needed
            plotter.add_text(f"Date and Time: {dt}", position='upper_left', font_size=10, color='black')

        # Optional: grid and view
        plotter.show_grid()
        plotter.view_isometric()
        plotter.show()

    except Exception as e:
        print(f"Error processing the GML: {e}")
