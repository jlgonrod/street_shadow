import os
from pandas import Timestamp as pd_Timestamp
from zoneinfo import ZoneInfo
import pickle as pkl
import numpy as np
import geopandas as gpd
import pandas as pd

from modules.graphs import (get_graph_from_osm,
                            load_graph_pkl,
                            get_nodes_edges,
                            apply_shadow_fractions,
                            get_new_weights,
                            add_weights_to_graph,
                            save_graph,
                            calculate_routes,
                            route_to_list_coordinates,
                            display_all_routes_on_map)
from modules.sun import get_sulight_vector
from modules.coordinates import convert_coordinates_EPSG_to_4326
from tqdm import tqdm

# GLOBAL VARIABLES
CITY = "malaga"
MESH_PATH = f"./data/processed_files/{CITY}_combined_mesh.vtk"
ALL_COORDS_PATH = f"./data/processed_files/{CITY}_all_coords_for_map.pkl"
GEOJSON_PATH = f"/mnt/d/JLGon/Descargas/street_shadow_data/shadow_geojson/{CITY}"
EPSG_SOURCE = "EPSG:25830"
filename = f"graph_base.pkl"
POLYGON_QUERY_GRAPH_PATH = f"/mnt/d/JLGon/Descargas/street_shadow_data/osmnx/{CITY}/{CITY}_polygon_geometry_to_query_graph.pkl"
GRAPH_BASE_PATH = f"/mnt/d/JLGon/Descargas/street_shadow_data/osmnx/{CITY}/{filename}"
DATETIME = "2025-06-01 18:05:00"

# FUNCTIONS
def load_base_graph():
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

    # Get needed data to calculate the sun vector
    dt = pd_Timestamp(DATETIME).tz_localize(ZoneInfo("Europe/Madrid"))
    all_coords = pkl.load(open(ALL_COORDS_PATH, 'rb'))
    x, y = np.mean(all_coords, axis=0) if all_coords else (0, 0)
    # Convert coordinates to EPSG:4326 (lat and lon)
    x, y = convert_coordinates_EPSG_to_4326(x, y, EPSG_SOURCE)

    sun_vector = get_sulight_vector(x, y, dt, "EPSG:4326", convert_coords=False) # Coordinates already in EPSG:4326

    # Check if the weighted graph already exists
    weighted_graph_path = GRAPH_BASE_PATH.replace("_base.pkl", f"_{sun_vector[0]}_{sun_vector[1]}_{sun_vector[2]}.pkl")
    
    if os.path.exists(weighted_graph_path):
        print("Weighted graph already exists, loading it...")
        G_weighted = load_graph_pkl(weighted_graph_path)

        # Get the alpha values from the graph
        edges = G_weighted.edges(data=True)
        edges = pd.DataFrame([attr for _, _, attr in edges])
        
        alpha_values = [col.split("_")[1] for col in edges.columns if col.startswith("weight_")]
        
    else:
        print("Weighted graph does not exist, creating it...")
        # Load the base graph
        G = load_base_graph()

        # Get nodes and edges
        nodes, edges = get_nodes_edges(G, GRAPH_BASE_PATH)

        # Get the geojson
        geojson_path = os.path.join(GEOJSON_PATH, f"{CITY}_{sun_vector[0]}_{sun_vector[1]}_{sun_vector[2]}.geojson")
        geojson_shadows = gpd.read_file(geojson_path)

        # Get the shadows fractions for each edge
        print("Calculating shadow fractions...")
        edges_geometries = edges["geometry"]
        edges["shadow_fraction"] = apply_shadow_fractions(geojson_shadows, edges_geometries)

        # Calculate updated edge weights using the alpha value range
        edges = get_new_weights(edges, 0.1)
        
        # New weights are stored in the graph
        print("Adding new weights to the graph...")
        G_weighted, alpha_values = add_weights_to_graph(G, edges)

        print("Saving graph with new weights...")
        # Save the graph with new weights
        save_graph(G_weighted, weighted_graph_path)


    # Calculate routes using the weighted graph
    print("Calculating routes...")
    origen = (36.711829, -4.431232) 
    destination = (36.726485, -4.421981)

    routes = calculate_routes(origen, destination, G_weighted, alpha_values)

    # Get the list of coordinates for each route
    routes_coords = {}
    for alpha, route in routes.items():
        list_coords = route_to_list_coordinates(route, G_weighted)
        routes_coords[alpha] = list_coords

    # Display all routes on a map
    display_all_routes_on_map(routes_coords, "all_routes_map.html")
    