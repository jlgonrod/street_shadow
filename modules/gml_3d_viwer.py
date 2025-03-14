import xml.etree.ElementTree as ET
import pyvista as pv
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union


from .map_image import draw_base_map
from .sun import get_sulight_vector

def unify_shadow_mesh(shadow_mesh: pv.PolyData) -> pv.PolyData:
    """
    Combina todas las caras solapadas en shadow_mesh en un único polígono.

    Parámetros
    ----------
    shadow_mesh : pv.PolyData
        Malla con caras posiblemente solapadas.

    Retorna
    -------
    pv.PolyData
        Malla limpia, unificada y sin solapamientos.
    """

    # Extraer polígonos de la malla de sombra
    faces = shadow_mesh.faces.reshape(-1, 4)[:, 1:]  # Extraer solo índices de vértices
    points = shadow_mesh.points[:, :2]  # Sólo coordenadas XY (planas)

    polygons = []
    for face in faces:
        polygon_coords = points[face]
        poly = Polygon(polygon_coords)
        if poly.is_valid and poly.area > 0:
            polygons.append(poly)

    # Unión de polígonos para eliminar solapamientos
    unified_poly = unary_union(polygons)

    # Puede haber multipolígonos, manejar esta posibilidad:
    if unified_poly.geom_type == 'MultiPolygon':
        all_coords, all_faces = [], []
        face_offset = 0
        for poly in unified_poly.geoms:
            x, y = poly.exterior.coords.xy
            coords_2d = np.column_stack((x, y))
            coords_3d = np.hstack((coords_2d, np.zeros((len(coords_2d), 1))))
            faces_poly = [[len(coords_3d)] + list(range(face_offset, face_offset + len(coords_3d)))]
            face_offset += len(coords_3d)

            all_coords.extend(coords_3d)
            all_faces.extend(faces_poly)
    else:
        x, y = unified_poly.exterior.coords.xy
        coords_2d = np.column_stack((x, y))
        coords_3d = np.hstack((coords_2d, np.zeros((len(coords_2d), 1))))
        all_coords = coords_3d
        all_faces = [[len(coords_3d)] + list(range(len(coords_3d)))]

    # Reconstruir la malla
    unified_shadow_mesh = pv.PolyData(np.array(all_coords), np.hstack(all_faces))

    # Opcionalmente triangula para garantizar simplicidad
    unified_shadow_mesh = unified_shadow_mesh.triangulate().clean()

    return unified_shadow_mesh

def project_mesh_onto_z0(mesh: pv.PolyData, direction: np.ndarray) -> pv.PolyData:
    """
    Proyecta la malla 'mesh' sobre el plano z=0 usando la 'direction' dada.
    
    Parámetros
    ----------
    mesh : pv.PolyData
        Malla a proyectar (superficie).
    direction : np.ndarray
        Vector de dirección (dx, dy, dz) para la proyección.
    
    Retorna
    -------
    shadow_mesh : pv.PolyData
        Nueva malla proyectada en z=0 con la misma conectividad.
    
    Lanza
    -----
    ValueError
        Si direction[2] == 0, pues no hay intersección con z=0.
    """
    if direction[2] == 0:
        raise ValueError("El componente z del vector de proyección es cero. "
                         "No puede proyectarse sobre z=0.")
    
    original_points = mesh.points.copy()
    shadow_points = []

    for p in original_points:
        # Ecuación: p + t * direction
        # Queremos z=0 => p.z + t*dz = 0 => t = -p.z/dz
        t = -p[2] / direction[2]
        p_proj = p + t * direction
        shadow_points.append(p_proj)

    shadow_points = np.array(shadow_points)

    # Construimos la nueva malla con la misma conectividad (mesh.faces)
    shadow_mesh = pv.PolyData(shadow_points, mesh.faces)

    # Combinamos las sombras solapadas
    shadow_mesh = unify_shadow_mesh(shadow_mesh)

    return shadow_mesh


def gml_3d_from_file(gml_file_path, dt, texture_map=True):
    """
    Processes a GML file and generates 3D models of buildings,
    projecting shadows on the z=0 plane if a datetime is provided.
    Optionally, a texture map can be added to the base plane.

    - Extrudes all buildings to a given height (3 m per floor)
    - Merges all buildings into a single 3D mesh
    - (Optional) Projects and unifies the shadow
    - (Optional) Adds a texture map to the base plane

    Parameters
    ----------
    gml_file_path : str
        Path to the GML file.
    dt : pd.Timestamp
        Date and time to calculate the shadows based on the sun position.
        It must have a timezone (e.g., Europe/Madrid).
        If None, no shadows are projected.
    texture_map : bool, opcional
        Whether to add a texture map to the base plane. Por defecto True.
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

        # List to store the 3D meshes of each building
        buildings_3d = []
        # List to store 2D coordinates (to draw the base texture)
        all_coords = []

        for building in root.findall(".//bu-ext2d:BuildingPart", ns):
            # Amount of floors above ground
            floors_above_element = building.find(".//bu-ext2d:numberOfFloorsAboveGround", ns)
            num_floors_above = int(floors_above_element.text) if floors_above_element is not None else 0

            if num_floors_above == 0:
                # No floors above ground, don't extrude it
                continue

            # Total height (3 m per floor, adjust according to your data)
            total_height = num_floors_above * 3
            base_z = 0

            # Get the base geometry (2D) and extrude
            for poslist in building.findall(".//gml:posList", ns):
                if poslist.text:
                    coords_text = poslist.text.strip()
                    coords_2d = np.array([
                        list(map(float, coords_text.split()[i:i+2]))
                        for i in range(0, len(coords_text.split()), 2)
                    ])

                    # Accumulate coords for the base bounding box (texture)
                    all_coords.extend(coords_2d)

                    # Ensure the polygon is closed
                    if not np.array_equal(coords_2d[0], coords_2d[-1]):
                        coords_2d = np.vstack([coords_2d, coords_2d[0]])
                    
                    # Generate 3D coordinates (z = base_z)
                    coords_3d = np.hstack([
                        coords_2d,
                        np.full((coords_2d.shape[0], 1), base_z)
                    ])

                    # Create a base polygon mesh
                    faces = [[len(coords_3d)] + list(range(len(coords_3d)))]
                    polydata = pv.PolyData(coords_3d, faces)

                    # Triangulate the base
                    polydata = polydata.triangulate()

                    # Extrude with capping=True (closes the top cap)
                    extruded = polydata.extrude((0, 0, total_height), capping=True)

                    # Extract complete surface and triangulate
                    # to ensure it is a closed solid in surface mesh
                    extruded = extruded.extract_surface().triangulate()

                    buildings_3d.append(extruded)

        # If there are buildings, unify for projection and connectivity
        if buildings_3d:
            combined_mesh = buildings_3d[0].copy()
            for bld in buildings_3d[1:]:
                combined_mesh = combined_mesh.merge(bld)

            # Create a plotter
            plotter = pv.Plotter()

            # Compute the center of the Bbox for the X and Y axes
            center = np.mean(all_coords, axis=0)

            # Add base plane with texture (optional)
            if texture_map:
                pad = 30
                base_plane, texture, center = draw_base_map(all_coords, pad) # Center is overwritten
                plotter.add_mesh(base_plane, texture=texture)

            # Add the unified mesh of all buildings
            plotter.add_mesh(
                combined_mesh,
                color="lightblue", opacity=1, show_edges=False,
                label="Edificios unificados"
            )

            # If a datetime is provided, calculate the shadow projection
            if dt:
                # Sunlight direction vector is calculated pointing to the origin
                sunlight_direction = get_sulight_vector(center[0], center[1], dt)

                # Shadows are projected on the z=0 plane
                shadow_mesh = project_mesh_onto_z0(combined_mesh, sunlight_direction)

                # Clean the mesh (remove unused points and cells)
                shadow_mesh = shadow_mesh.triangulate().clean()

                # Add the shadow mesh to the plot
                plotter.add_mesh(
                    shadow_mesh,
                    color="gray", opacity=0.8, show_edges=False,
                    label="Sombra proyectada"
                )

            # Add datetime information as text (if provided)
            if dt:
                plotter.add_text(f"Date and Time: {dt}", position='upper_left', font_size=10, color='black')

            # Configure the plotter
            plotter.show_grid()
            plotter.view_isometric()
            plotter.show()

        else:
            print("No 3D models were generated.")

    except Exception as e:
        print(f"Error processing the GML: {e}")
