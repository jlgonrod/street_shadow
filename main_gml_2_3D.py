from modules.gml_3d_viwer import gml_3d_from_file
from pandas import Timestamp as pd_Timestamp
from zoneinfo import ZoneInfo

# Get the GML file path for a building
#ref = "3981901PB8238S"
#GML_FILE_PATH = f"./data/gml/buildings/Building_{ref}.gml"

# Get the GML file path for a town
GML_FILE_PATH = "./data/gmL/towns/A.ES.SDGC.BU.21011.buildingpart.gml"

# Set the datetime to calculate the shadows based on the sun position
dt_sample = pd_Timestamp(
    year=2023,
    month=3,
    day=1,
    hour=17,
    minute=45,
    second=00,
    tz=ZoneInfo('Europe/Madrid') # It cares about the change of hour in summer and winter in Spain
)

gml_3d_from_file(GML_FILE_PATH,dt=dt_sample, texture_map=False)