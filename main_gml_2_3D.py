from modules.gml_3d_viwer import gml_3d_from_file
from pandas import Timestamp as pd_Timestamp
from zoneinfo import ZoneInfo

# Get the GML file path for a building
#ref = "3981901PB8238S"
#GML_FILE_PATH = f"./data/gml/buildings/Building_{ref}.gml"

# Get the GML file path for a town
name = "candon"
GML_FILE_PATH = f"./data/gml/towns/{name}.gml"

# Set the datetime to calculate the shadows based on the sun position
dt_sample = pd_Timestamp(
    year=2025,
    month=3,
    day=14,
    hour=15,
    minute=0,
    second=0,
    tz=ZoneInfo('Europe/Madrid') # It cares about the change of hour in summer and winter in Spain
)

gml_3d_from_file(GML_FILE_PATH,dt=dt_sample, texture_map=False)