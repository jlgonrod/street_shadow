from modules.catastro import get_building_gml

# Define the reference of the parcel
ref_cat = "7118005TG3471N"
save_folder = "data/gml/buildings"

get_building_gml(ref_cat, save_folder)