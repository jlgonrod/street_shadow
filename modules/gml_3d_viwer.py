import xml.etree.ElementTree as ET
import pyvista as pv
import numpy as np

def gml_3d_from_file(gml_file_path):
    """
    Procesa un archivo GML y genera modelos 3D de los edificios.

    :param gml_file_path: Ruta del archivo GML
    """
    try:
        # Cargar el archivo GML
        tree = ET.parse(gml_file_path)
        root = tree.getroot()

        # Definir namespaces
        ns = {
            "gml": "http://www.opengis.net/gml/3.2",
            "bu-ext2d": "http://inspire.jrc.ec.europa.eu/schemas/bu-ext2d/2.0"
        }

        # Lista para almacenar los modelos 3D de los edificios
        buildings_3d = []

        for building in root.findall(".//bu-ext2d:BuildingPart", ns):
            # Obtener número de plantas sobre rasante
            floors_above_element = building.find(".//bu-ext2d:numberOfFloorsAboveGround", ns)
            num_floors_above = int(floors_above_element.text) if floors_above_element is not None else 0

            # Si no hay pisos sobre rasante, se ignora esta parte del edificio
            if num_floors_above == 0:
                continue

            # Definir altura
            height_above = num_floors_above * 3  # 3m por planta
            base_z = 0  # Base en el suelo
            total_height = height_above  

            # Obtener la geometría base (2D) y extruirla a 3D
            for poslist in building.findall(".//gml:posList", ns):
                if poslist.text:
                    coords_text = poslist.text.strip()
                    coords_2d = np.array([list(map(float, coords_text.split()[i:i+2])) for i in range(0, len(coords_text.split()), 2)])
                    
                    # Asegurar que el polígono está cerrado
                    if not np.array_equal(coords_2d[0], coords_2d[-1]):
                        coords_2d = np.vstack([coords_2d, coords_2d[0]])
                    
                    # Convertir coordenadas 2D a 3D añadiendo la base Z correcta
                    coords_3d = np.hstack([coords_2d, np.full((coords_2d.shape[0], 1), base_z)])
                    
                    # Crear un objeto de malla cerrada
                    faces = [[len(coords_3d)] + list(range(len(coords_3d)))]
                    polydata = pv.PolyData(coords_3d, faces)
                    
                    # Triangular el polígono para manejar polígonos complejos
                    polydata = polydata.triangulate()
                    
                    # Extruir el polígono en altura asegurando que sea un sólido cerrado
                    extruded = polydata.extrude((0, 0, total_height), capping=True)
                    
                    buildings_3d.append(extruded)

        # Si hay modelos generados, visualizarlos
        if buildings_3d:
            plotter = pv.Plotter()
            for mesh in buildings_3d:
                plotter.add_mesh(mesh, color="lightblue", opacity=1, show_edges=False)

            plotter.show_grid()
            plotter.view_isometric()
            plotter.show()
        else:
            print("No se generaron modelos 3D.")

    except Exception as e:
        print(f"Error al procesar el GML: {e}")