import os
import re
import pandas as pd
import pickle as pkl
from tqdm import tqdm
import geopandas as gpd

from modules.graphs import (load_graph_pkl,
                            process_graph_using_geojson,
                            save_graph,
                            get_graph_from_osm,
                            get_nodes_edges)

# GLOBAL VARIABLES
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

CITY = "malaga"
FOLDER_SHADOWS = os.path.join(BASE_DIR, "data", "shadow_geojson", CITY)
FOLDER_GRAPHS = os.path.join(BASE_DIR, "data", "osmnx", CITY)
GRAPH_BASE_PATH = os.path.join(BASE_DIR, "data", "osmnx", CITY, "graph_base.pkl")
GRAPH_BASE_EDGES = os.path.join(BASE_DIR, "data", "osmnx", CITY, "_edges_gdfs.pkl")

MESH_PATH = os.path.join(BASE_DIR, "data", "processed_files", f"{CITY}_combined_mesh.vtk")
EPSG_SOURCE = "EPSG:25830"
POLYGON_QUERY_GRAPH_PATH = os.path.join(BASE_DIR, "data", "osmnx", CITY, f"{CITY}_polygon_geometry_to_query_graph.pkl")

# FUNCTIONS
def get_data_basename(basename):
    splitted = re.split(r"_", os.path.splitext(basename)[0])
    x, y, z = splitted[-3:]
    return basename, x, y, z

def load_or_build_base_graph():
    """
    Load the base graph from a file or download it from OSM if it does not exist.
    
    Parameters
    ----------
    graph_base_path : str
        Path to the graph file.
        Example: "/mnt/d/JLGon/Descargas/street_shadow_data/osmnx/malaga/malaga_base.pkl"
        
    Returns
    -------
    graph : networkx.Graph
        Graph object.
    """
    if os.path.exists(GRAPH_BASE_PATH):
        # Load the graph from file
        print("Graph already exists, loading it...")
        graph = load_graph_pkl(GRAPH_BASE_PATH)
    else:
        # Check if there is a polygon file to query the graph
        if os.path.exists(POLYGON_QUERY_GRAPH_PATH):
            with open(POLYGON_QUERY_GRAPH_PATH, 'rb') as f:
                polygon = pkl.load(f)
        else:
            polygon = None

        # Download the graph from OSM using osmnx library
        print("Graph does not exist, downloading it from OSM...")
        graph = get_graph_from_osm(GRAPH_BASE_PATH, polygon, MESH_PATH, EPSG_SOURCE)

    return graph

if __name__ == "__main__":
    
    # Get the list of all geojson files in the directory
    pattern = re.compile(rf"{CITY}_[+-]?\d+\.\d+_[+-]?\d+\.\d+_[+-]?\d+\.\d+\.geojson")
    shadows_files = [f for f in os.listdir(FOLDER_SHADOWS) if pattern.match(f)]

    # Create a DataFrame to store the shadows data of the basename
    shadows_files = [get_data_basename(f) for f in shadows_files]
    shadows_files_df = pd.DataFrame(shadows_files, columns=["basename", "x", "y", "z"])

    # Get the list of all graph (except base) files in the directory
    pattern = re.compile(rf"graph_[+-]?\d+\.\d+_[+-]?\d+\.\d+_[+-]?\d+\.\d+\.pkl")
    graphs_files = [f for f in os.listdir(FOLDER_GRAPHS) if pattern.match(f)]
    
    graphs_files = [get_data_basename(f) for f in graphs_files]
    graphs_files_df = pd.DataFrame(graphs_files, columns=["basename", "x", "y", "z"])

    # conserva solo los geojson que no tengan su graph por coincidencia de x, y, z
    shadows_without_graph = shadows_files_df[~shadows_files_df[["x", "y", "z"]].apply(tuple, 1).isin(graphs_files_df[["x", "y", "z"]].apply(tuple, 1))]
    shadows_without_graph.reset_index(drop=True, inplace=True)

    # Create or load the base_graph and base edges
    G = load_or_build_base_graph()
    nodes, edges = get_nodes_edges(G, GRAPH_BASE_PATH)

    # Iterate over the geojson files without graph
    for geojson_file, x, y, z in tqdm(shadows_without_graph.values, desc="Processing geojson files"):
        
        # Get the geojson path
        geojson_path = os.path.join(FOLDER_SHADOWS, geojson_file)

        # Process the graph
        G_weighted = process_graph_using_geojson(G, edges, geojson_path)

        # Save the graph
        graph_save_path = os.path.join(FOLDER_GRAPHS, f"graph_{x}_{y}_{z}.pkl")
        save_graph(G_weighted, graph_save_path)

