from flask import Flask, request, render_template, url_for, send_from_directory, Response
import os
import pickle as pkl
import numpy as np
import geopandas as gpd
import pandas as pd
from pandas import Timestamp as pd_Timestamp
from zoneinfo import ZoneInfo
import sys
import folium
import threading
import queue

original_sys_path = sys.path.copy()
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)
from modules.graphs import (
    load_graph_pkl,
    get_nodes_edges,
    calculate_routes,
    remove_repeated_routes,
    route_to_list_coordinates,
    process_routes_distances,
    route_time_from_distances,
    get_all_routes_on_map,
    save_and_format_map_html
)
from modules.sun import get_sulight_vector, get_existing_sun_vector
from modules.coordinates import convert_coordinates_EPSG_to_4326
sys.path = original_sys_path

app = Flask(__name__)

# GLOBAL VARIABLES
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, ".."))
CITY = "malaga"
MESH_PATH = os.path.join(PROJECT_ROOT, "data", "processed_files", f"{CITY}_combined_mesh.vtk")
ALL_COORDS_PATH = os.path.join(PROJECT_ROOT, "data", "processed_files", f"{CITY}_all_coords_for_map.pkl")
GEOJSON_PATH = os.path.join(PROJECT_ROOT, "data", "shadow_geojson", CITY)
EPSG_SOURCE = "EPSG:25830"
GRAPH_BASE_PATH = os.path.join(PROJECT_ROOT, "data", "osmnx", CITY, "graph_base.pkl")
USER_SPEED = 4.7  # km/h

# Create a queue to store messages
messages_queue = queue.Queue()

def log_message(msg):
    """
    Log a message to the console and put it in the messages queue.
    """
    print(msg)
    messages_queue.put(msg)

def load_or_build_base_graph():
    """
    Load the base graph from a file or download it from OSM if it does not exist.
    """
    POLYGON_QUERY_GRAPH_PATH = os.path.join(PROJECT_ROOT, "data", "osmnx", CITY, f"{CITY}_polygon_geometry_to_query_graph.pkl")
    if os.path.exists(GRAPH_BASE_PATH):
        log_message("The base graph exists, loading it...")
        graph = load_graph_pkl(GRAPH_BASE_PATH)
    else:
        polygon = None
        if os.path.exists(POLYGON_QUERY_GRAPH_PATH):
            with open(POLYGON_QUERY_GRAPH_PATH, 'rb') as f:
                import pickle as pkl
                polygon = pkl.load(f)
        from modules.graphs import get_graph_from_osm
        log_message("The base graph does not exist, downloading it from OSM...")
        graph = get_graph_from_osm(GRAPH_BASE_PATH, polygon, MESH_PATH, EPSG_SOURCE)
    return graph

def process_request(origen, destination, date, time):
    """
    Process the request to calculate routes and generate a map.
    
    Parameters
    ----------
    origen : str
        Origin address in as a string.
    destination : str
        Destination address in as a string.
    date : str
        Date in the format YYYY-MM-DD.
    time : str
        Time in the format HH:MM.
    """

    dt_str = f"{date} {time}:00"
    dt = pd_Timestamp(dt_str).tz_localize(ZoneInfo("Europe/Madrid"))
    
    # Load all coordinates from the pickle file to get the center
    with open(ALL_COORDS_PATH, "rb") as f:
        all_coords = pkl.load(f)
    x, y = np.mean(all_coords, axis=0) if all_coords else (0, 0)
    x, y = convert_coordinates_EPSG_to_4326(x, y, EPSG_SOURCE)

    # Get the sun vector for the given coordinates and datetime
    sun_vector = get_sulight_vector(x, y, dt, "EPSG:4326", convert_coords=False)

    # Ensure the sun vector matches an existing one or find a similar precomputed vector
    sun_vector = get_existing_sun_vector(sun_vector, GEOJSON_PATH, max_mse_allowed=0.01)
    
    # Get the name of the weighted graph based on the sun vector and load it
    weighted_graph_path = GRAPH_BASE_PATH.replace("_base.pkl", f"_{sun_vector[0]}_{sun_vector[1]}_{sun_vector[2]}.pkl")
    
    if os.path.exists(weighted_graph_path):
        # Load the weighted graph from the file
        log_message("Loading weighted graph...")
        G_weighted = load_graph_pkl(weighted_graph_path)
        log_message("Extracting alpha values from the graph...")
        edges = G_weighted.edges(data=True, keys=True)
        edges = gpd.GeoDataFrame.from_records(
            [(u, v, key, data) for u, v, key, data in edges],
            columns=["u", "v", "key", "attributes"]
        )
        attributes_df = pd.DataFrame(edges["attributes"].tolist())
        edges = pd.concat([edges.drop(columns=["attributes"]), attributes_df], axis=1)
        edges.set_index(["u", "v", "key"], inplace=True)
        alpha_values = [col.split("_")[1] for col in edges.columns if col.startswith("weight_")]
    
    else:
        # Load the base graph and create the weighted graph
        log_message("Creating weighted graph...")
        G = load_or_build_base_graph()
        nodes, edges = get_nodes_edges(G, GRAPH_BASE_PATH)
        geojson_path = os.path.join(
            GEOJSON_PATH,
            f"{CITY}_{sun_vector[0]}_{sun_vector[1]}_{sun_vector[2]}.geojson"
        )
        geojson_shadows = gpd.read_file(geojson_path)
        log_message("Applying shadow fractions and weights to edges...")
        edges["shadow_fraction"] = __import__('modules.graphs').graphs.apply_shadow_fractions(geojson_shadows, edges["geometry"])
        edges = __import__('modules.graphs').graphs.get_new_weights(edges, 0.1)
        G_weighted, alpha_values = __import__('modules.graphs').graphs.add_weights_and_shadow_fractions_to_graph(G, edges)
        __import__('modules.graphs').graphs.save_graph(G_weighted, weighted_graph_path)
    
    # Generate the routes
    log_message("Generating routes...")
    routes = calculate_routes(origen, destination, G_weighted, alpha_values)

    # Remove repeated routes
    log_message("Removing repeated routes...")
    routes = remove_repeated_routes(routes)
    routes_coords = {}
    
    # Convert routes from a list of OSMID to a list of
    # coordinates adding the real origin and destination points
    log_message("Converting routes to coordinates...")
    for alpha, route in routes.items():
        routes_coords[alpha] = route_to_list_coordinates(origen, destination, route, G_weighted)

    # Compute the routes metrics    
    log_message("Calculating distances...")
    routes_distances = process_routes_distances(routes, edges)
    log_message("Calculating times...")
    routes_times = route_time_from_distances(routes_distances, USER_SPEED)
    
    # Load the shadows as a Geojson file
    log_message("Reading geojson file for shadows...")
    geojson_path = os.path.join(
        GEOJSON_PATH,
        f"{CITY}_{sun_vector[0]}_{sun_vector[1]}_{sun_vector[2]}.geojson"
    )
    shadows = gpd.read_file(geojson_path)

    # Geenrate the map with the routes and shadows
    log_message("Creating map with routes and shadows...")
    map_with_routes = get_all_routes_on_map(routes_coords, shadows, routes_distances)
    
    # Save the map to the static folder
    log_message("Saving map to static folder...")
    static_folder = os.path.join(BASE_DIR, "static")
    os.makedirs(static_folder, exist_ok=True)
    map_html_path = os.path.join(static_folder, "map_with_routes.html")
    save_and_format_map_html(map_with_routes, dt, CITY, origen, destination, routes_coords, routes_times, False, map_html_path)
    
    with open(map_html_path, "r", encoding="utf-8") as f:
        content = f.read()
    modified_content = content.replace('assets/', '/assets/')
    with open(map_html_path, "w", encoding="utf-8") as f:
        f.write(modified_content)
    
    log_message("Process completed.")

@app.route("/progress")
def progress():
    def generate():
        while True:
            msg = messages_queue.get()
            yield f"data: {msg}\n\n"
    return Response(generate(), mimetype="text/event-stream")

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        origen = request.form.get("origen")
        destination = request.form.get("destination")
        date = request.form.get("date")   # format: YYYY-MM-DD
        time = request.form.get("time")   # format: HH:MM
        
        # Start the heavy task in a background thread
        thread = threading.Thread(target=process_request, args=(origen, destination, date, time))
        thread.start()
        # Respond immediately so the form does not block
        return "Process started", 200
    else:
        static_folder = os.path.join(BASE_DIR, "static")
        os.makedirs(static_folder, exist_ok=True)
        map_html_path = os.path.join(static_folder, "map_with_routes.html")
        # Coordinates of Málaga
        malaga_coords = [36.72093, -4.42404]

        # Set the default map
        default_map = folium.Map(location=malaga_coords, zoom_start=12, zoom_control="topright")
        default_map.save(map_html_path)
        return render_template("index.html", map_url=url_for('static', filename="map_with_routes.html"))

@app.route("/assets/<path:filename>")
def assets(filename):
    assets_path = os.path.join(PROJECT_ROOT, "assets")
    return send_from_directory(assets_path, filename)

if __name__ == "__main__":
    app.run(debug=True)