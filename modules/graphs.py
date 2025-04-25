from time import time
from modules.coordinates import get_mesh_bounds_coords
import pyvista as pv
import pickle
import numpy as np
import osmnx as ox
import folium
import numpy as np
import pandas as pd
from shapely.geometry import Polygon,MultiPolygon, box
from shapely.strtree import STRtree
import branca.colormap as color_map
import os
import geopandas as gpd
import json

BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def load_graph_pkl(graph_path):
    """
    Load a graph from a pickle file.
    """
    start_time = time()
    with open(graph_path, 'rb') as f:
        graph = pickle.load(f)
    end_time = time()
    print("\tGraph loading time: ", end_time - start_time)
    return graph

def save_graph(graph, graph_path):
    """
    Save a graph to a pickle file.
    """
    os.makedirs(os.path.dirname(graph_path), exist_ok=True)
    with open(graph_path, 'wb') as f:
        pickle.dump(graph, f)
    print(f"Graph saved to {graph_path}")

def get_graph_from_osm(graph_base_path, polygon=None, mesh_path=None, epsg_mesh_source=None):
    """
    Get a graph from OSM using the mesh bounds and save it to a file.

    Parameters
    ----------
    graph_base_path : str
        Path to save the graph file.
        Example: "data/osmnx/malaga/malaga_base.pkl"
    polygon : shapely.geometry.Polygon, optional
        Polygon to use for graph extraction. If None, mesh bounds will be used.
        Polingons coordinates must be in EPSG:4326.
    mesh_path : str, optional
        Path to the mesh file. Required if polygon is None.
    epsg_mesh_source : str, optional
        EPSG code of the mesh source. Required if polygon is None.
    
    Returns
    -------
    graph : networkx.Graph
        Graph object.
    """

    if polygon is not None:
        print("Polygon provided, using it to get the graph.")
        print("Downloading graph from OSM...")
        # Get the graph from OSM using the polygon
        start_time = time()
        graph = ox.graph_from_polygon(polygon,
                                    network_type='walk',
                                    simplify=False,
                                    retain_all=True,
                                    truncate_by_edge=False
                                    )
        end_time = time()
        print("\tGraph download and generation time: ", end_time - start_time)

    else:
        print("Polygon not provided, using mesh bounds.")
        print("Loading mesh...")
        mesh = pv.read(mesh_path)

        coordinates_bbox = np.round(np.array(get_mesh_bounds_coords(mesh, epsg_mesh_source)).flatten(), 4).tolist()
        print("\tCoordinates BBOX: ", coordinates_bbox)

        print("Downloading graph from OSM...")
        # Get the graph from OSM using the bounds
        start_time = time()
        graph = ox.graph_from_bbox(coordinates_bbox,
                                network_type='walk',
                                simplify=False,
                                retain_all=True,
                                truncate_by_edge=False
                                )
        end_time = time()
        print("\tGraph download and generation time: ", end_time - start_time)

    # Save the graph to a file
    print("Saving graph to file...")
    save_graph(graph, graph_base_path)

    return graph

def get_nodes_edges(G, graph_base_path):
    """
    Get nodes and edges from a graph and save them to files
    if they do not exist.
    
    Parameters
    ----------
    G : networkx.Graph
        Graph to get nodes and edges from.
    graph_base_path : str
        Path to the graph base file.
        Example: "data/osmnx/malaga/malaga_base.pkl"
    
    Returns
    -------
    nodes : GeoDataFrame
        Nodes of the graph.
    edges : GeoDataFrame
        Edges of the graph.
    """

    # Check if nodes and edges already exist
    nodes_path = graph_base_path.replace("graph_base.pkl", "_nodes_gdfs.pkl")
    edges_path = graph_base_path.replace("graph_base.pkl", "_edges_gdfs.pkl")
    exist_nodes = os.path.exists(nodes_path)
    exist_edges = os.path.exists(edges_path)

    # If nodes and edges already exist, load them
    if exist_nodes and exist_edges:
        print("Nodes and edges already exist, loading them...")
        nodes = pickle.load(open(nodes_path, 'rb'))
        edges = pickle.load(open(edges_path, 'rb'))
        return nodes, edges
    # If nodes and edges do not exist, generate them
    else:
        print("Nodes or edges do not exist, generating them...")

        nodes, edges = ox.graph_to_gdfs(G,
                                        nodes=True,
                                        edges=True,
                                        node_geometry=True,
                                        fill_edge_geometry=True
                                        )
        
        # Save nodes and edges to files overwriting if one of them exists
        print("Saving nodes...")
        with open(nodes_path, 'wb') as f:
            pickle.dump(nodes, f)
        print("Saving edges...")
        with open(edges_path, 'wb') as f:
            pickle.dump(edges, f)
        print("Nodes and edges saved to files.")

        return nodes, edges
    
def add_shadows_to_map(map, shadows_geojson):
    """
    Add shadows to a folium map from a geojson file.
    
    Parameters
    ----------
    map : folium.Map
        Folium map to add shadows to.
    shadows_geojson : gpd.GeoDataFrame
        GeoDataFrame with the shadows.
    Returns
    -------
    map : folium.Map
        Folium map with the shadows added. 
    """
    
    # Add the shadows to the map
    folium.GeoJson(
        shadows_geojson,
        name="shadows",
        style_function=lambda x: {
            "fill_color": "black",
            "color": "black",
            "weight": 0.5,
            "fillOpacity": 0.25
        },
    ).add_to(map).get_root().html.add_child(folium.Element('<div data-type="shadow"></div>'))
    
    return map
    
def geojson_to_individual_polygons(shadows_geojson):
    """
    Convert a GeoDataFrame with polygons and multipolygons to a list of individual polygons.

    Parameters
    ----------
    shadows_geojson : GeoDataFrame
        GeoDataFrame with polygons and multipolygons.

    Returns
    -------
    polygons : list
        List of individual polygons.
    """

    polygons = []
    # Iterate over each feature in the GeoJSON
    for geometry in shadows_geojson.geometry:
        
        # Check if the geometry is not empty
        if geometry.is_empty or geometry is None:
            continue

        # If it's a polygon, add it to the list
        if isinstance(geometry, Polygon):
            polygons.append(geometry)

        # If it's a multipolygon, flat it and append to polygons
        elif isinstance(geometry, MultiPolygon):
            polygons.extend([g for g in geometry.geoms if not g.is_empty])

    return polygons

def get_shadow_fractions(geometry, strtree, list_poligons):
    try:
        # Check if the geometry is not empty and is valid
        if geometry.is_empty or not geometry.is_valid:
            return 0.0
        
        canditates = strtree.query(geometry)

        # Chef if there are no candidates
        if canditates.size == 0:
            return 0.0
        
        # Verify if it's candidate index or geometries
        if isinstance(canditates[0], (int, np.integer)):
            # If it's a list of indexes, get the geometries
            canditates = [list_poligons[i] for i in canditates]

        # Intersect the candidates with the geometry and remove empty ones
        inter = [s.intersection(geometry) for s in canditates if s.intersects(geometry)]
        inter = [segment for segment in inter if not segment.is_empty]

        # If there are no intersections, return 0.0
        if not inter:
            return 0.0
        
        # Calculate the shadow longitude
        longitude_shadow = sum([segment.length for segment in inter])

        # Calculate the shadow fraction and return it or 1.0 if it's greater than 1.0
        return min(1.0, longitude_shadow / geometry.length)
        
    except Exception as e:
        print(f"Error calculating shadow fraction: {e}")
        return 0.0

def apply_shadow_fractions(shadows_geojson, edges_geom_series):
    """
    This function calculates the shadow fractions for each edge in the graph
    using the shadows from a geojson file.

    Parameters
    ----------
    shadows_geojson : GeoDataFrame
        GeoDataFrame with the shadows.
    edges_geom_series : pandas.Series
        Series with the edges geometries.

    Returns
    -------
    pandas.Series
        Series with the shadow fractions for each edge.
    """

    list_polygons = geojson_to_individual_polygons(shadows_geojson)
    tree = STRtree(list_polygons)
    
    return edges_geom_series.apply(get_shadow_fractions, args=(tree, list_polygons))
    
def add_shadow_fraction_segments_to_map(map, edges):
    """
    
    Parameters
    ----------
    map : folium.Map
        Folium map to add shadows to.
    edges : GeoDataFrame
        GeoDataFrame with the edges. It must contain
        the geometry and the shadow_fraction columns.

    Returns
    -------
    map : folium.Map
        Folium map with the shadows fracction segments added.
    """
    
    # Define the colormap
    cmap = color_map.LinearColormap(
        colors=["#ff271c", "#1795d4"],
        vmin=0, vmax=1,
        caption='Fracción de sombra')
    
    # Convert the edges GeoDataFrame to GeoJSON
    edges_geojson = edges[["geometry", "shadow_fraction"]]
    edges_geojson = edges_geojson.to_json()

    # Add a geojson layer to the map with a style function
    def style_function(feature):
        shadow_fraction = feature['properties']['shadow_fraction']
        color = cmap(shadow_fraction)
        return {
            'color': color,
            'weight': 3,
            'opacity': 0.7
        }
    
    # Add the GeoJSON layer to the map
    folium.GeoJson(
        edges_geojson,
        style_function=style_function,
        name='Shadow Fractions'
    ).add_to(map)
    
    return map

def compute_new_weight(base_weight, shadow_fraction, alpha):
    """
    Compute a new weight for the edges based on the shadow fraction
    and a base weight (length). The new weight considers both the 
    total distance and how much of it is under sunlight or shadow.

    Parameters
    ----------
    base_weight : float
        The original weight of the edge, typically its length.
    shadow_fraction : float
        The fraction of the edge's length that is in shadow (0.0 to 1.0).
    alpha : float
        A parameter (0.0 to 1.0) that defines the preference for shadow or distance:
        - alpha = 1.0: Prioritize shorter distance, even if under sunlight.
        - alpha = 0.5: Balance between distance and shadow.
        - alpha = 0.0: Prioritize shadow, even if it means walking a longer distance.

    Returns
    -------
    float
        The new weight for the edge, adjusted based on the shadow fraction and alpha.
    """

    return base_weight * (alpha + (1 - alpha) * (1 - shadow_fraction))

def get_new_weights(edges ,alpha_res=0.1):
    """
    This function calculates the new weights for each edge in the graph
    based on the shadow fraction and a base weight (length). The new weights
    are computed using the compute_new_weight function, which considers both
    the total distance and how much of it is under sunlight or shadow.

    The new weights are stored in the edges GeoDataFrame with the format
    new_weight_{alpha}, where alpha is the parameter that defines the preference
    for shadow or distance. The alpha values are defined in the range [0, 1]
    with a resolution of alpha_res.

    Parameters
    ----------
    edges : GeoDataFrame
        GeoDataFrame with the edges. It must contain
        the geometry and the shadow_fraction columns.
        It must also contain the length column. (base
        weight)
    alpha_res : float
        The resolution of the alpha values. The default is 0.1.
        The alpha values are defined in the range [0, 1]
        with a resolution of alpha_res.

    Returns
    -------
    edges : GeoDataFrame
        GeoDataFrame with the new weights for each edge.
        The new weights are stored in the edges GeoDataFrame
        with the format new_weight_{alpha}, where alpha is
        the parameter that defines the preference for shadow
        or distance.
    """

    # Define all the possible alpha values that are available to choose by the user
    alpha_values = np.arange(0, 1 + alpha_res, alpha_res)[:-1]  # Remove the last value (1.0) the same weight as the base weight

    for alpha in alpha_values:
        rounded_alpha = round(alpha, 2)  # Round alpha to 2 decimal places avoid floating point issues
        # Compute the new weights for each edge
        edges[f'weight_{rounded_alpha}'] = edges.apply(
            lambda row: compute_new_weight(row['length'], row['shadow_fraction'], alpha), axis=1
        )

    # Rename "length" column to "weight_1.0" to keep the same format
    edges.rename(columns={"length": "weight_1.0"}, inplace=True)

    # Move the "weight_1.0" column to the end of the DataFrame
    cols = list(edges.columns)
    cols.append(cols.pop(cols.index("weight_1.0")))
    edges = edges[cols]

    return edges

def add_weights_and_shadow_fractions_to_graph(G, edges_weight):
    """
    This function adds the new weights and shadow fractions to the graph. 
    The new weights are stored in the edges GeoDataFrame with the format weight_{alpha}.
    The shadow fractions are stored in the shadow_fraction column.
    Edges weights must contain the u, v, key as columns or index. And at least
    one column with the format weight_{alpha}.

    Parameters
    ----------
    G : osmnx.graph
        Graph to add the weights and shadow fractions to.
    edges_weight : GeoDataFrame
        GeoDataFrame with the edges. It must contain
        the geometry and the shadow_fraction columns.
        It must also contain the length column. (base
        weight)
        It must also contain the u, v, key columns or index.

    Returns
    -------
    G : osmnx.graph
        Graph with the new weights and shadow fractions added.
        The weights are stored in the columns weight_{alpha} for each alpha value.
        The shadow fractions are stored in the shadow_fraction column.
    alpha_values : list
        List with the alpha values used to calculate the new weights.
        The alpha values are defined in the range [0, 1] with a resolution
        of alpha_res.
    """

    # Get all available alpha values
    alpha_values = [col.split("_")[1] for col in edges_weight.columns if col.startswith('weight')]

    # Check alpha values is not empty
    if not alpha_values:
        raise ValueError("No alpha values found in the edges GeoDataFrame.")

    # New weights and shadow fractions are stored in the graph
    edges_weight.reset_index(inplace=True)  # u, v, key are now columns

    # Create a DataFrame from the graph edges to facilitate the merge
    graph_edges = pd.DataFrame(
        [(u, v, k) for u, v, k in G.edges(keys=True)],
        columns=["u", "v", "key"]
    )

    # Perform a merge between the DataFrame of the graph edges and the `edges_weight` DataFrame
    merged_edges = graph_edges.merge(edges_weight, on=["u", "v", "key"], how="left")

    # Update the weights and shadow fractions directly in the graph
    for _, row in merged_edges.iterrows():
        u, v, k = row["u"], row["v"], row["key"]
        data = G[u][v][k]
        # Add shadow_fraction
        if pd.notna(row["shadow_fraction"]):  # Verify that the value is not NaN
            data["shadow_fraction"] = row["shadow_fraction"]
        # Add weights for each alpha
        for alpha in alpha_values:
            if pd.notna(row[f"weight_{alpha}"]):  # Verify that the value is not NaN
                data[f"weight_{alpha}"] = row[f"weight_{alpha}"]

    return G, alpha_values

def calculate_routes(origen, destination, G, alpha_list):
    """
    
    Parameters
    ----------
    orgine : str
        Address of the origin point.
        Example: "Calle Rios Rosas 1, Malaga, España"
    destination : str
        Address of the destination point.
        Example: "Calle Purificación 4, Malaga, España"
    G : osmnx.graph
        Graph to calculate the routes from.
    alpha_list : list
        List with the alpha values used to calculate the new weights.
        The alpha values are defined in the range [0, 1].
    """

    # Get the coordinates of the origin and destination
    origen_coords = ox.geocoder.geocode(origen)
    destination_coords = ox.geocoder.geocode(destination)

    # Get the nearest nodes in the graph
    org_node = ox.distance.nearest_nodes(G, origen_coords[1], origen_coords[0])
    dest_node = ox.distance.nearest_nodes(G, destination_coords[1], destination_coords[0])

    # Calculate the routes for each alpha value
    routes = {}
    for alpha in alpha_list:
        # Calculate the route
        route = ox.shortest_path(G, org_node, dest_node, weight=f'weight_{alpha}', cpus=6)
        routes[alpha] = route

    return routes

def route_to_list_coordinates(origen, destination, route, G):
    """
    This function takes a list of nodes as the list of their identifiers
    and returns a list of tuples with the coordinates of each node. It
    is important to note that all the nodes identified by their id must
    be in the graph. If not, the function will raise an error.    
    
    For each node, a tuple with the coordinates (lat, lon) is returned.
    
    Parameters
    ----------
    origen : str
        Address of the origin point.
        Example: "Calle Rios Rosas 1, Malaga, España"
    destination : str
        Address of the destination point.
        Example: "Calle Purificación 4, Malaga, España"
    route : list
        List of nodes identifiers.
    G : osmnx.graph
        Graph to get the coordinates from.

    Returns
    -------
    list_coordinates
        List of tuples with the coordinates (lat, lon) of each node.
        The coordinates are in EPSG:4326.
    """

    org_coords = ox.geocoder.geocode(origen)
    dest_coords = ox.geocoder.geocode(destination)

    coordinates_list = [[G.nodes[node]['y'], G.nodes[node]['x']] for node in route]

    # If the first node is not the origin, add the origin coordinates
    if coordinates_list[0] != [org_coords[0], org_coords[1]]:
        coordinates_list.insert(0, [org_coords[0], org_coords[1]])

    # If the last node is not the destination, add the destination coordinates
    if coordinates_list[-1] != [dest_coords[0], dest_coords[1]]:
        coordinates_list.append([dest_coords[0], dest_coords[1]])

    return coordinates_list

def remove_repeated_routes(routes_list):
    """
    This function takes a dictionary with alpha values as keys
    and routes (list of nodes) as values. It removes the repeated
    routes from the dictionary. The routes are considered repeated
    if they have the same nodes in the same order. The function
    returns a dictionary with the unique routes. The keys are the
    alpha values and the values are the unique routes.

    Parameters
    ----------
    routes_list : dict
        Dictionary with the alpha values as keys and the routes
        as values. The routes are lists of nodes.

    Returns
    -------
    unique_routes : dict
        Dictionary with the unique routes. The keys are the alpha
        values and the values are the unique routes.
    """

    unique_routes = {}

    # Iterate over the routes
    for alpha, route in routes_list.items():
        # Check if the route is already in the unique_routes dictionary
        if route in unique_routes.values():
            continue
        else:
            # If the route is not in the unique_routes dictionary, add it
            unique_routes[alpha] = route

    return unique_routes

def max_min_center_coords_routes(routes_coords, pad=0):
    """
    This function takes a dictionary with the routes coordinates
    (one for each alpha value) and returns the max and min latitude
    and longitude of the routes. It also returns the center of the
    coordinates. The coordinates are in EPSG:4326.

    A padding is added to the max and min coordinates to avoid
    clipping the routes. The padding is in degrees.

    Parameters
    ----------
    routes_coords : dict
        Dictionary with the routes coordinates. The keys are the alpha
        values and the values are lists of tuples with the coordinates
        (lat, lon) of each node.
    pad : float
        Padding to add to the max and min coordinates. The padding is
        in degrees. The default is 0.

    Returns
    -------
    center_coords : tuple
        Tuple with the center coordinates (lat, lon).
    max_lat : float
        Maximum latitude of the routes.
    min_lat : float
        Minimum latitude of the routes.
    max_lon : float
        Maximum longitude of the routes.
    min_lon : float
        Minimum longitude of the routes.
    """
    # Get a plain list with the tuples of coordinates
    coords = [coord for route in routes_coords.values() for coord in route]
    coords = np.array(coords)

    # Get the max and min lat and lon and add a padding
    max_lat = coords[:, 0].max() + pad
    min_lat = coords[:, 0].min() - pad
    max_lon = coords[:, 1].max() + pad
    min_lon = coords[:, 1].min() - pad

    # Get the center of the coordinates
    center_lat = (max_lat + min_lat) / 2
    center_lon = (max_lon + min_lon) / 2

    return (center_lat, center_lon), max_lat, min_lat, max_lon, min_lon

def add_route_to_map(route, alpha, dist_dict, color, map):
    """
    This function takes a list of coordinates and adds
    a polyline to the folium map representing the route.

    Parameters
    ----------
    route : list
        List of coordinates (lat, lon) of the route.
    aplha : str
        Alpha value of the route. It is used to identify
        the route in the map.
    dist_dict : dict
        Dictionary with the distances of the route. The keys
        are the alpha values and the values are nested dictionaries
        with the distances of shadow and sun.
    color : str
        Color of the polyline. It can be a hex color code
        or a color name. Example: "#ff271c" or "red".
    map : folium.Map
        Folium map to add the polyline to. The map must
        be already created and centered on the coordinates
        of the routes. The map must be in EPSG:4326.

    Returns
    -------
    None
        The function does not return anything. It adds the
        polyline to the map.
    """
    # Create a feature group for the route
    route_layer = folium.FeatureGroup(name=alpha)

    # Create a tooltip with the distances
    shadow_distance = dist_dict["distance_shadow"]
    sun_distance = dist_dict["distance_sun"]
    dist_info = f"""
    <div style="text-align: left;">
        <b>Shadow distance:</b>&emsp;{int(shadow_distance)} m<br>
        <b>Sun distance:</b>&emsp;&nbsp;&nbsp;&nbsp;{int(sun_distance)} m<br>
        <b>Total distance:</b>&emsp;&nbsp;&nbsp;&nbsp;{int(shadow_distance + sun_distance)} m
    </div>
    """

    # Add the route to the layer
    folium.PolyLine(
        locations=route,
        color=color,
        weight=5,
        opacity=0.7,
        tooltip=dist_info
    ).add_to(route_layer)

    # Add the layer to the map
    route_layer.add_to(map)

def get_all_routes_on_map(routes_coords, shadows, route_distances):
    """
    This function takes a dictionary with the routes coordinates
    (one for each alpha value) and displays them on a folium map.
    The map is saved to a HTML file.

    Parameters
    ----------
    routes_coords : dict
        Dictionary with the routes coordinates. The keys are the alpha
        values and the values are lists of tuples with the coordinates
        (lat, lon) of each node. Coordinates are in EPSG:4326.
    shadows : GeoDataFrame
        GeoDataFrame with the shadows. It must contain the geometry
        column. The geometry must be in EPSG:4326. If None, the shadows
        will not be added to the map.
    route_distances : dict
        Dictionary with the distances of each route. The keys are the
        alpha values and the values are nested dictionaries with the
        distances of shadow and sun. The distances are in meters.

    Returns
    -------
    None
        The function does not return anything. It saves the map to
        a HTML file.
    """
    
    # Get the needed coordinates to load the map
    center_coords, max_lat, min_lat, max_lon, min_lon = max_min_center_coords_routes(routes_coords)
    
    # Create a folium map centered on the coordinates
    map = folium.Map(location=center_coords, zoom_start=14)
    map.fit_bounds([[min_lat, min_lon], [max_lat, max_lon]])

    # Generate a list of colors for the routes between two values interpolating len(routes_coords) colors
    colors = color_map.LinearColormap(
        colors=["#1795d4", "#ff271c"],
        vmin=0, vmax=len(routes_coords)-1,
    )
    list_colors = [colors(i) for i in range(len(routes_coords))]

    # Ensure the routes are sorted in ascending order (most shadow first)
    routes_coords = dict(sorted(routes_coords.items(), key=lambda item: item[0], reverse=False)) 

    # Add the routes to the map
    for i, (alpha, route) in enumerate(routes_coords.items()):
        distances = route_distances[alpha]
        add_route_to_map(route, alpha, distances, list_colors[i], map)

    if shadows is not None:
        # Add the shadows to the map but cutting using the min and max lat and lon
        pad = 0.001 # Pad in degrees not in meters
        bounding_box = box(min_lon - pad, min_lat - pad, max_lon + pad, max_lat + pad)
        shadows = shadows.clip(bounding_box)
        shadows = shadows.reset_index(drop=True)

        folium.GeoJson(
            shadows,
            name="shadows",
            style_function=lambda x: {
                "fill_color": "black",
                "color": "black",
                "weight": 0.5,
                "fillOpacity": 0.25
            },
        ).add_to(map)

    return map

def load_assets(filepath):
    """
    This function loads the assets from the assets folder.
    The assets folder is in the same directory as this file.

    Parameters
    ----------
    filepath : str
        Path to the asset file. The path is relative to the
        assets folder. Example: "/templates/panel.html"

    Returns
    -------
    str
        The content of the asset file.
    """
    with open(filepath, 'r') as f:
        return f.read()

def save_and_format_map_html(map, datetime, city, origen, destination, routes_coords, dist_times, info_panel, map_path_html):

    # Format the datetime into a readable string
    datetime_str = datetime.strftime("%d/%m/%Y %H:%M:%S")

    # Add a panel to the map with the datetime and city
    if info_panel:
        panel_path = os.path.join(BASE_DIR, "assets", "templates", "panel.html")
        panel_template = load_assets(panel_path)
        panel_html = panel_template.format(
            datetime_str=datetime_str,
            city=city.capitalize(),
            origen=origen,
            destination=destination
        )
        map.get_root().html.add_child(folium.Element(panel_html))

    # Add markers for the origin and destination
    folium.Marker(
        location=ox.geocoder.geocode(origen),
        popup=origen,
        icon=folium.Icon(color="gray", icon="play"),
    ).add_to(map).get_root().html.add_child(folium.Element('<div data-type="marker"></div>'))

    folium.Marker(
        location=ox.geocoder.geocode(destination),
        popup=destination,
        icon=folium.Icon(color="gray", icon="stop"),
    ).add_to(map).get_root().html.add_child(folium.Element('<div data-type="marker"></div>'))

    # Add JavaScript for slider and checkbox
    slider_js_path = os.path.join(BASE_DIR, "assets", "templates", "slider.js")
    slider_js = load_assets(slider_js_path)
    map.get_root().html.add_child(folium.Element(f"<script>{slider_js}</script>"))

    # Add slider and checkbox to the map
    slider_path = os.path.join(BASE_DIR, "assets", "templates", "slider.html")
    slider_template = load_assets(slider_path)
    slider_html = slider_template.format(
        len_routes=len(routes_coords) - 1
        )
    map.get_root().html.add_child(folium.Element(slider_html))

    # Display the times for each route
    
    # Add the times box (when only one route is selected)
    times_info_path = os.path.join(BASE_DIR, "assets", "templates", "times_info.html")
    times_info_html = load_assets(times_info_path)
    map.get_root().html.add_child(folium.Element(times_info_html))

    # Order times by route alpha value (asceding)
    sorted_keys = sorted(dist_times, key=lambda k: float(k))
    route_times_array = [dist_times[k] for k in sorted_keys]

    # Add the times to the map
    script = f"<script>var routeTimes = {json.dumps(route_times_array)};</script>"
    map.get_root().html.add_child(folium.Element(script))

    # Load and add the JavaScript for the times box
    times_info_js_path = os.path.join(BASE_DIR, "assets", "templates", "times_info.js")
    times_info_js = load_assets(times_info_js_path)
    map.get_root().html.add_child(folium.Element(f"<script>{times_info_js}</script>"))

    # Save the map to a HTML file
    map.save(map_path_html)

def process_graph_using_geojson(G, edges, geojson_path):
    """
    This function processes a graph using a geojson file with shadows.
    It calculates the shadow fractions for each edge in the graph
    and updates the edge weights based on the shadow fractions.
    The new weights are stored in the graph with the format
    weight_{alpha}, where alpha is the parameter that defines the
    preference for shadow or distance. The alpha values are defined
    in the range [0, 1] with a resolution of alpha_res.

    Parameters
    ----------
    G : osmnx.graph
        Graph base to process.
    edges : GeoDataFrame
        Edges in the graph. It must contain the geometry column.
    geojson_path : str
        Path to the geojson file with the shadows.
        Example: "data/shadow_geojson/malaga/malaga_36.711829_-4.431232_0.0.geojson"

    Returns
    -------
    G_weighted : osmnx.graph
        Graph with the new weights added stored in the columns
        weight_{alpha} for each alpha value.
    """

    # Copy the edges to avoid modifying the original GeoDataFrame
    edges_copy = edges.copy()

    # Load the geojson file
    shadows_geojson = gpd.read_file(geojson_path)

    # Get the shadow fractions
    edges_copy["shadow_fraction"] = apply_shadow_fractions(shadows_geojson, edges_copy["geometry"])

    # Calculate updated edge weights using the alpha value range
    edges_copy = get_new_weights(edges_copy, 0.1) # Resolution of 0.1

    # Store the new weights in the graph
    G_weighted, *_ = add_weights_and_shadow_fractions_to_graph(G, edges_copy)

    return G_weighted

def get_route_distance(route, lengths_edges):
    """
    This function calculates the distance of a route in shadow and sun.
    The route is a list of nodes and the lengths_edges is a DataFrame
    with the lengths of the edges. The lengths_edges DataFrame must
    contain the columns weight_1.0 and shadow_fraction. The weight_1.0
    column is the length of the edge in meters. The shadow_fraction
    column is the fraction of the edge that is in shadow. The indexs
    of the lengths_edges DataFrame must be the u and v columns of the
    edges GeoDataFrame.

    Parameters
    ----------
    route : list
        List of nodes in the route. The nodes are identified by their
        identifiers.
    lengths_edges : DataFrame
        DataFrame with the lengths of the edges.
        
    Returns
    -------
    distance_shadow : int
        Distance of the route in shadow in meters.
    distance_sun : int
        Distance of the route in sun in meters.

    """
    # Convert route to a DataFrame of consecutive node pairs
    route_df = pd.DataFrame({'u': route[:-1], 'v': route[1:]})

    # Merge the route DataFrame with the lengths_edges DataFrame
    merged_edges = route_df.merge(lengths_edges, on=['u', 'v'], how='left')

    # Calculate the distances in shadow and sun
    distance_shadow = (merged_edges["weight_1.0"] * merged_edges["shadow_fraction"]).sum()
    distance_sun = (merged_edges["weight_1.0"] * (1 - merged_edges["shadow_fraction"])).sum()

    # Round the values with no decimals
    distance_shadow = int(round(distance_shadow, 0))
    distance_sun = int(round(distance_sun, 0))

    return distance_shadow, distance_sun

def process_routes_distances(routes, edges):
    """
    This function calculates the distances of the routes in shadow and sun.
    The routes are a dictionary with the alpha values as keys and the
    routes as values.

    Parameters
    ----------
    routes : dict
        Dictionary with the alpha values as keys and the routes
        as values. The routes are lists of nodes.
    edges : GeoDataFrame
        GeoDataFrame with the edges. It must contain the u, v
        as index (origen and destination nodes) to identify the edges.
        It also must contain the weight_1.0 column, which is the length
        of the edge in meters and the shadow_fraction column, which is the
        fraction of the edge that is in shadow.

    Returns
    -------
    distances : dict
        Dictionary with the alpha values as keys and an nested
        dictionary with the distances as values. The nested dictionary
        contains the keys distance_shadow and distance_sun, which are
        the distances of the route in shadow and sun, and the values
        are the distances in meters.
    """
    # Get the length of each edge in the graph
    edges_length = edges[["weight_1.0", "shadow_fraction"]] # Length in meters is the base weight

    # Create a dictionary to store the distances
    distances = {}

    # Iterate over the routes
    for alpha, route in routes.items():
        # Get the length of each edge in the route
        dist_shadow, dist_sun = get_route_distance(route, edges_length)

        # Store the distances in the dictionary
        distances[alpha] = {
            "distance_shadow": dist_shadow,
            "distance_sun": dist_sun
        }

    return distances
        
def route_time_from_distances(routes_distances_dic, speed_kmh=5):
    """
    This function calculates the time of each route in shadow and sun.
    The routes are a dictionary with the alpha as keys and the values
    are a nested dictionary with the distances at sun and shadow.
    
    Parameters
    ----------
    routes_distances_dic : dict
        Dictionary with the alpha values as keys and a nested
        dictionary with the distances as values.
    speed_kmh : int or float
        Speed in km/h. The speed is used to calculate the time
        of the route. The default is 5 km/h.

    Returns
    -------
    routes_times : dict
        Dictionary with the alpha values as keys and a nested
        dictionary with the times as values. The nested dictionary
        contains the keys time_shadow and time_sun, which are the
        times of the route in shadow and sun, and the values are
        the times in seconds.
    """
    # Convert speed from km/h to m/s
    speed_mps = speed_kmh * 1000 / 3600

    # Create a dictionary to store the times
    routes_times = {}

    # Iterate over the routes distances
    for alpha, distances in routes_distances_dic.items():
        # Calculate the time in shadow and sun
        time_shadow = distances["distance_shadow"] / speed_mps
        time_sun = distances["distance_sun"] / speed_mps

        # Get the time in minutes
        time_shadow = time_shadow / 60
        time_sun = time_sun / 60

        # Round the times to integers
        time_shadow = int(round(time_shadow, 0))
        time_sun = int(round(time_sun, 0))

        # Store the times in the dictionary
        routes_times[alpha] = {
            "time_shadow": time_shadow,
            "time_sun": time_sun
        }

    return routes_times
