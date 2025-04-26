import os
from modules.catastro import get_building_gml

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Define the reference of the parcel
ref_cat = "3981901PB8238S"
save_folder = os.path.join(BASE_DIR, "data", "gml", "buildings")

get_building_gml(ref_cat, save_folder)