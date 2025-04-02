import pandas as pd
from zoneinfo import ZoneInfo
from pvlib.solarposition import sun_rise_set_transit_ephem
from modules.sun import get_sulight_vector
from modules.coordinates import convert_coordinates_EPSG_to_4326
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import pickle as pkl
import numpy as np
import os

# GLOBAL VARIABLES
YEAR = 2025
ALL_COORDS_PATH = "./data/processed_files/candon_all_coords_for_map.pkl"
epsg_source = "EPSG:25830"  # CHANGE THIS TO THE EPSG CODE OF THE COORDINATES IN THE GML FILE

### To execute this script, it is necessary to have a precomputed `all_coords_for_map.pkl`
### file for the town or building of interest. This file is used to determine the center
### coordinates of the building or observation point for calculating the sun vectors.
### If the file is not available, can be generated executing the script main_gml_2_3D.py

# Load the coordinates from the pickle file and get the center coordinates
if not os.path.exists(ALL_COORDS_PATH):
    raise FileNotFoundError(f"File {ALL_COORDS_PATH} not found. Please generate it first excuting the script main_gml_2_3D.py")

all_coords = pkl.load(open(ALL_COORDS_PATH, "rb")) 

center_xy = np.mean(all_coords, axis=0) if all_coords else (0, 0)
longitude, latitude = convert_coordinates_EPSG_to_4326(center_xy[0], center_xy[1], epsg_source)


## For each day of the year, get the sun rise and sun set times (based on the timezone of the location) ##
# Create a DatetimeIndex for all days of the year
date_range = pd.date_range(
    start=f'{YEAR}-01-01',
    end=f'{YEAR}-12-31',
    freq='D',
    tz=ZoneInfo('Europe/Madrid')  # Localize to Madrid timezone
)

# Calculate sun rise and set times for each day
sun_times = sun_rise_set_transit_ephem(
    date_range,
    latitude,
    longitude,
).drop(columns=['transit'])  # keep only sunrise and sunset times

# Round to second precision
sun_times['sunrise'] = sun_times['sunrise'].dt.ceil('s')
sun_times['sunset'] = sun_times['sunset'].dt.ceil('s')

## Get each instant of time between sunrise and sunset with the desired frequency ##
freq = '1min'  # CHANGE THIS TO THE DESIRED FREQUENCY
sun_times_list = []

for index, row in sun_times.iterrows():
    # Create a date range for each second between sunrise and sunset
    daily_range = pd.date_range(
        start=row['sunrise'],
        end=row['sunset'],
        freq=freq,
        tz=ZoneInfo('Europe/Madrid')  # Localize to Madrid timezone, CHANGE THIS TO THE DESIRED TIMEZONE
    )
    
    # Append the daily range to the list
    sun_times_list.append(daily_range)

## Get the sert of unique sun light vectors for each second of the year ##
sun_vectors = set()

def process_daily_range(daily_range):
    daily_vectors = set()
    for dt in daily_range:
        # Get the light vector for each second
        x, y, z = get_sulight_vector(
            x=longitude,
            y=latitude,
            dt=dt,
            convert_coords=False  # Coordinates already in EPSG:4326
        )
        daily_vectors.add((x, y, z))  # Use a tuple to ensure uniqueness
    return daily_vectors

if __name__ == '__main__':
    # Use multiprocessing to parallelize the computation
    with Pool(cpu_count()) as pool:
        results = list(tqdm(pool.imap(process_daily_range, sun_times_list), 
                            desc="Processing daily ranges", 
                            total=len(sun_times_list)))

    # Combine the sets of results to ensure uniqueness
    for daily_vectors in results:
        sun_vectors.update(daily_vectors)

    # Save the values in a csv file use the year and loc in the name
    sun_vectors_df = pd.DataFrame(list(sun_vectors), columns=['x', 'y', 'z'])
    sun_vectors_df.to_csv(f'./data/sun_vectors/sun_vectors_{YEAR}_{latitude}_{longitude}.csv', index=False)