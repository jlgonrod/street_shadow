import os
import pickle as pkl
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from zoneinfo import ZoneInfo
from tqdm import tqdm

# Third-party packages
from pvlib.solarposition import sun_rise_set_transit_ephem

# Internal modules
from modules.coordinates import convert_coordinates_EPSG_to_4326
from modules.sun import get_sulight_vector

# GLOBAL CONFIGURATION
CITY = "malaga"
EPSG_SOURCE = "EPSG:25830"
FREQ = '1min'
TIMEZONE = 'Europe/Madrid'
OUTPUT_DIR = './data/sun_vectors'
START_DATETIME = "2025-06-01 18:00:00"  # changed from START_DATE
END_DATETIME = "2025-06-01 18:05:00"      # changed from END_DATE

ALL_COORDS_PATH = f"./data/processed_files/{CITY}_all_coords_for_map.pkl"

def load_coordinates(filepath):
    """Carga las coordenadas y calcula el punto central."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"File {filepath} not found. Please generate it first by executing main_gml_2_3D.py"
        )
    all_coords = pkl.load(open(filepath, "rb"))
    center_xy = np.mean(all_coords, axis=0) if all_coords else (0, 0)
    return center_xy

def get_date_range(start_datetime, end_datetime, tz):
    """Genera un rango de fechas (un día por entrada) usando la parte de fecha de los datetime."""
    start_date = pd.to_datetime(start_datetime).date()
    end_date = pd.to_datetime(end_datetime).date()
    return pd.date_range(
        start=start_date,
        end=end_date,
        freq='D',
        tz=ZoneInfo(tz)
    )

def get_sun_times(date_range, latitude, longitude):
    """Calcula y redondea los tiempos de amanecer y atardecer."""
    sun_times = sun_rise_set_transit_ephem(date_range, latitude, longitude).drop(columns=['transit'])
    sun_times['sunrise'] = sun_times['sunrise'].dt.ceil('s')
    sun_times['sunset'] = sun_times['sunset'].dt.ceil('s')
    return sun_times

def generate_daily_ranges(sun_times):
    """Genera intervalos de tiempo entre amanecer y atardecer y retorna además sunrise y sunset."""
    ranges = []
    for _, row in sun_times.iterrows():
        daily_range = pd.date_range(
            start=row['sunrise'],
            end=row['sunset'],
            freq=FREQ,
            tz=ZoneInfo(TIMEZONE)
        )
        ranges.append((daily_range, row['sunrise'], row['sunset']))
    return ranges

def process_daily_range(daily_range, sunrise, sunset, longitude, latitude, global_start, global_end):
    """Calcula los vectores de luz únicos para un rango diario, ignorando instantes no visibles y fuera del rango global."""
    daily_vectors = set()
    for dt in daily_range:
        if dt <= sunrise or dt >= sunset:
            continue
        if dt < global_start or dt > global_end:
            continue
        x, y, z = get_sulight_vector(
            x=longitude,
            y=latitude,
            dt=dt,
            epsg_source=EPSG_SOURCE,
            convert_coords=False
        )
        daily_vectors.add((x, y, z))
    return daily_vectors

def process_daily_range_wrapper(args):
    # Función wrapper para evitar pickling de lambdas.
    daily_range, sunrise, sunset, longitude, latitude, global_start, global_end = args
    return process_daily_range(daily_range, sunrise, sunset, longitude, latitude, global_start, global_end)

if __name__ == '__main__':
    center_xy = load_coordinates(ALL_COORDS_PATH)
    longitude, latitude = convert_coordinates_EPSG_to_4326(center_xy[0], center_xy[1], EPSG_SOURCE)
    
    # Se utiliza el rango de fechas definido por START_DATETIME y END_DATETIME
    date_range = get_date_range(START_DATETIME, END_DATETIME, TIMEZONE)
    sun_times = get_sun_times(date_range, latitude, longitude)
    
    # Guardar CSV con datos de sunrise y sunset en columnas separadas
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    sun_times_output_path = f'{OUTPUT_DIR}/sun_times_{CITY}_{START_DATETIME.split(" ")[0]}_{END_DATETIME.split(" ")[0]}.csv'
    formatted_sun_times = pd.DataFrame({
        'date': sun_times['sunrise'].dt.strftime('%Y-%m-%d'),
        'sunrise': sun_times['sunrise'].dt.strftime('%H:%M:%S'),
        'sunset': sun_times['sunset'].dt.strftime('%H:%M:%S')
    })
    formatted_sun_times.to_csv(sun_times_output_path, index=False)
    print(f"Sun times saved to {sun_times_output_path}")
    
    daily_ranges = generate_daily_ranges(sun_times)
    
    # Definir el rango global de datetimes usando START_DATETIME y END_DATETIME
    global_start = pd.Timestamp(START_DATETIME, tz=ZoneInfo(TIMEZONE))
    global_end = pd.Timestamp(END_DATETIME, tz=ZoneInfo(TIMEZONE))
    
    # Procesamiento en paralelo utilizando todos los núcleos disponibles
    with Pool(cpu_count()) as pool:
        tasks = [(dr, sr, ss, longitude, latitude, global_start, global_end) for (dr, sr, ss) in daily_ranges]
        results = list(tqdm(
            pool.imap(process_daily_range_wrapper, tasks),
            desc="Processing daily ranges",
            total=len(daily_ranges)
        ))

    sun_vectors = set()
    for vectors in results:
        sun_vectors.update(vectors)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path = f'{OUTPUT_DIR}/sun_vectors_{CITY}_{START_DATETIME.replace(" ", "T")}_{END_DATETIME.replace(" ", "T")}.csv'
    pd.DataFrame(list(sun_vectors), columns=['x', 'y', 'z']).to_csv(output_path, index=False)
    print(f"Sunlight vectors saved to {output_path}")