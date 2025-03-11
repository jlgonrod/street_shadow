import xml.etree.ElementTree as ET
import pyvista as pv
import numpy as np
from .map_image import draw_base_map

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

    return shadow_mesh


def gml_3d_from_file(gml_file_path, texture_map=True, projection_direction=None):
    """
    Processes a GML file and generates 3D models of buildings, opcionalmente
    proyectando una "sombra" de cada edificio sobre el plano z=0.

    - Fusiona todas las edificaciones en una sola malla 3D
    - (Opcional) Proyecta y unifica la sombra para no ver patches sobrepuestos

    Parámetros
    ----------
    gml_file_path : str
        Path to the GML file.
    texture_map : bool, opcional
        Whether to add a texture map to the base plane. Por defecto True.
    projection_direction : np.ndarray, opcional
        Vector de proyección para calcular la sombra sobre z=0.
        Ejemplo: np.array([0,0,-1]).
        Por defecto None (no se proyecta).
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

        # Lista para almacenar las mallas 3D de cada edificio
        buildings_3d = []
        # Lista para almacenar coordenadas 2D (para dibujar la textura base)
        all_coords = []

        for building in root.findall(".//bu-ext2d:BuildingPart", ns):
            # Nº de pisos sobre el suelo
            floors_above_element = building.find(".//bu-ext2d:numberOfFloorsAboveGround", ns)
            num_floors_above = int(floors_above_element.text) if floors_above_element is not None else 0

            if num_floors_above == 0:
                # Sin pisos sobre el suelo, no lo extruimos
                continue

            # Altura total (3 m por piso, ajusta según tus datos)
            total_height = num_floors_above * 3
            base_z = 0

            # Obtener la geometría base (2D) y extruir
            for poslist in building.findall(".//gml:posList", ns):
                if poslist.text:
                    coords_text = poslist.text.strip()
                    coords_2d = np.array([
                        list(map(float, coords_text.split()[i:i+2]))
                        for i in range(0, len(coords_text.split()), 2)
                    ])

                    # Acumular coords para el bounding box de la base (textura)
                    all_coords.extend(coords_2d)

                    # Asegurar que el polígono está cerrado
                    if not np.array_equal(coords_2d[0], coords_2d[-1]):
                        coords_2d = np.vstack([coords_2d, coords_2d[0]])
                    
                    # Generar coords en 3D (z = base_z)
                    coords_3d = np.hstack([
                        coords_2d,
                        np.full((coords_2d.shape[0], 1), base_z)
                    ])

                    # Crear una malla de polígono base
                    faces = [[len(coords_3d)] + list(range(len(coords_3d)))]
                    polydata = pv.PolyData(coords_3d, faces)

                    # Triangular la base
                    polydata = polydata.triangulate()

                    # Extruir con capping=True (cierra tapa superior)
                    extruded = polydata.extrude((0, 0, total_height), capping=True)

                    # Extraer superficie completa y triangular
                    # para que sea un sólido cerrado en malla superficial
                    extruded = extruded.extract_surface().triangulate()

                    buildings_3d.append(extruded)

        # Si hay edificios, unificamos para proyección y conectividad
        if buildings_3d:
            combined_mesh = buildings_3d[0].copy()
            for bld in buildings_3d[1:]:
                combined_mesh = combined_mesh.merge(bld)

            # Visualización
            plotter = pv.Plotter()

            # Añadir plano base con textura (opcional)
            if texture_map:
                pad = 30
                base_plane, texture = draw_base_map(all_coords, pad)
                plotter.add_mesh(base_plane, texture=texture)

            # Añadimos la malla de todos los edificios unificada
            plotter.add_mesh(
                combined_mesh,
                color="lightblue", opacity=1, show_edges=False,
                label="Edificios unificados"
            )

            # Si se proporcionó vector para proyección, generamos sombra
            if projection_direction is not None:
                shadow_mesh = project_mesh_onto_z0(combined_mesh, projection_direction)
                
                # "Limpiamos" para unificar posibles solapes de la malla
                shadow_mesh = shadow_mesh.triangulate().clean()

                # Añadir la sombra
                plotter.add_mesh(
                    shadow_mesh,
                    color="gray", opacity=0.8, show_edges=False,
                    label="Sombra proyectada"
                )

            plotter.show_grid()
            plotter.view_isometric()
            plotter.show()

        else:
            print("No 3D models were generated.")

    except Exception as e:
        print(f"Error processing the GML: {e}")
