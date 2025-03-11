from modules.gml_3d_viwer import gml_3d_from_file
import numpy as np

ref = "3526101TG3432N"
GML_FILE_PATH = f"./data/gml/buildings/Building_{ref}.gml"
gml_3d_from_file(GML_FILE_PATH, texture_map=True, projection_direction=np.array([0.7, 0.5, -1]))