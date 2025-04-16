import os
import numpy as np
import pandas as pd
import pyvista as pv
import xml.etree.ElementTree as ET
from os.path import splitext, basename
from glob import glob

from modules.gml_3d_viwer import load_or_process_buildings, process_footprints
from modules.shadow import process_shadows, save_shadows_to_geojson
from tqdm import tqdm

# CONFIGURE YOUR PATHS HERE
GML_FILE_PATH = "./data/gml/towns/malaga.gml"
SUN_VECTORS_CSV = "./data/sun_vectors/malaga/sun_vectors_malaga_2025-06-01T18:00:00_2025-06-01T19:00:00.csv"
REMOVE_BASES = True
SAVE_PATH = f"/mnt/d/JLGon/Descargas/street_shadow_data/shadow_geojson/{splitext(basename(GML_FILE_PATH))[0]}/" # The files will be saved here.

if __name__ == "__main__":
    # Data source configuration
    tree = ET.parse(GML_FILE_PATH)
    root = tree.getroot()
    ns = {
        "gml": "http://www.opengis.net/gml/3.2",
        "bu-ext2d": "http://inspire.jrc.ec.europa.eu/schemas/bu-ext2d/2.0"
    }
    epsg_source = root.find(".//gml:Surface", ns).attrib.get("srsName", "").split(":")[-1] # Assume all the coordinates are in the same EPSG in the file.
    epsg_source = f"EPSG:{epsg_source}"

    # Load geometry files
    combined_mesh, footprints_polygons, all_coords_for_map = load_or_process_buildings(GML_FILE_PATH, root, ns)
    all_buildings_footprints = process_footprints(footprints_polygons) # Process the footprints

    # Get the center of the combined mesh
    min_x, min_y = np.min(all_coords_for_map, axis=0)
    max_x, max_y = np.max(all_coords_for_map, axis=0)
    coords_min_max = np.array([[min_x, min_y], [max_x, max_y]])
    center_xy = np.mean(coords_min_max, axis=0)
    
    # Load solar vectors
    sun_vectors_df = pd.read_csv(SUN_VECTORS_CSV)

    # Retrieve existing geojson files for the city
    existing_files_pattern = f"{SAVE_PATH}{splitext(basename(GML_FILE_PATH))[0]}_*.geojson"
    existing_files = glob(existing_files_pattern)
    existing_vectors = []
    for file in existing_files:
        bname = os.path.splitext(basename(file))[0]
        parts = bname.split('_')
        if len(parts) >= 4:
            existing_vectors.append({'x': float(parts[1]), 'y': float(parts[2]), 'z': float(parts[3])})
    existing_df = pd.DataFrame(existing_vectors)
    print(f"Found {len(existing_df)} existing vectors in {SAVE_PATH}")
    
    # Remove solar vectors that already exist from the sun_vectors_df dataframe
    if not existing_df.empty:
        merged = sun_vectors_df.merge(existing_df, on=["x", "y", "z"], how='left', indicator=True)
        sun_vectors_df = merged[merged['_merge'] == 'left_only'].drop(columns=['_merge'])

    # For each vector, calculate the shadow and save it as a geojson.
    for index, row in tqdm(sun_vectors_df.iterrows(), total=len(sun_vectors_df), desc="Processing shadows with sun vectors", unit="sun vector"):

        # Convert the vector to an array
        sun_vector = np.array([row["x"], row["y"], row["z"]])

        # Create the pyvista plotter
        plotter = pv.Plotter()

        # Calculate the shadows
        shadow_mesh = process_shadows(combined_mesh, sun_vector)

        # Save as geojson
        save_shadows_to_geojson(
            shadow_mesh,
            f"{SAVE_PATH}{splitext(basename(GML_FILE_PATH))[0]}_{sun_vector[0]}_{sun_vector[1]}_{sun_vector[2]}.geojson",
            all_buildings_footprints,
            epsg_source,
            remove_bases=REMOVE_BASES,
            verbose=True
        )