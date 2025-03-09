from modules.gml_3d_viwer import gml_3d_from_file

ref = "3526101TG3432N"
GML_FILE_PATH = f"./data/gml/buildings/Building_{ref}.gml"
gml_3d_from_file(GML_FILE_PATH, texture_map=False)