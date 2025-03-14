from modules.catastro import get_building_gml

# Define the reference of the parcel
ref_cat = "3981901PB8238S"
save_folder = "data/gml/buildings"

get_building_gml(ref_cat, save_folder)