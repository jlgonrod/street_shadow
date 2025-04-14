from time import time
from modules.coordinates import get_mesh_bounds_coords
import pyvista as pv
import pickle
import numpy as np
import osmnx as ox
import folium
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon
from shapely.strtree import STRtree
import branca.colormap as color_map
import os

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
        Example: "/mnt/d/JLGon/Descargas/street_shadow_data/osmnx/malaga/malaga_base.pkl"
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
        Example: "/mnt/d/JLGon/Descargas/street_shadow_data/osmnx/malaga/malaga_base.pkl"
    
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
    ).add_to(map)
    
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

def add_weights_to_graph(G, edges_weight):
    """
    This function adds the new weights to the graph. The new weights
    are stored in the edges GeoDataFrame with the format weight_{alpha}.
    Edges weigths must contain the u, v, key as columns or index. And at least
    one column with the format weight_{alpha}.

    Parameters
    ----------
    G : osmnx.graph
        Graph to add the weights to.
    edges_weight : GeoDataFrame
        GeoDataFrame with the edges. It must contain
        the geometry and the shadow_fraction columns.
        It must also contain the length column. (base
        weight)
        It must also contain the u, v, key columns or index.

    Returns
    -------
    G : osmnx.graph
        Graph with the new weights added stored in the columns
        weight_{alpha} for each alpha value.
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

    # New weights are stored in the graph
    edges_weight.reset_index(inplace=True) # u, v, key are now columns

    # Create a DataFrame from the graph edges to facilitate the merge
    graph_edges = pd.DataFrame(
        [(u, v, k) for u, v, k in G.edges(keys=True)],
        columns=["u", "v", "key"]
    )

    # Perform a merge between the DataFrame of the graph edges and the `edges_weight` DataFrame
    merged_edges = graph_edges.merge(edges_weight, on=["u", "v", "key"], how="left")

    # Update the weights directly in the graph
    for _, row in merged_edges.iterrows():
        u, v, k = row["u"], row["v"], row["key"]
        data = G[u][v][k]
        for alpha in alpha_values:
            if pd.notna(row[f"weight_{alpha}"]):  # Verificar que el valor no sea NaN
                data[f"weight_{alpha}"] = row[f"weight_{alpha}"]

    return G, alpha_values

def save_custom_graph(weighted_graph_path, graph_with_weights):
    """
    Save a graph with custom weights to a file.

    Parameters
    ----------
    weighted_graph_path : str
        Path to save the graph with weights.
        Example: "/mnt/d/JLGon/Descargas/street_shadow_data/osmnx/malaga/malaga_weighted.pkl"
    graph_with_weights : osmnx.graph
        Graph object with custom weights added as attributes
        (e.g., weight_{alpha} for each alpha value).

    Returns
    -------
    None
    """
    
    # Save the graph to a file
    os.makedirs(os.path.dirname(weighted_graph_path), exist_ok=True)
    with open(weighted_graph_path, 'wb') as f:
        pickle.dump(graph_with_weights, f)

def calculate_routes(orgine, destination, G, alpha_list):
    """
    
    Parameters
    ----------
    orgine : tuple
        Coordinates of the origin point (latitude, longitude).
    destination : tuple
        Coordinates of the destination point (latitude, longitude).
    G : osmnx.graph
        Graph to calculate the routes from.
    alpha_list : list
        List with the alpha values used to calculate the new weights.
        The alpha values are defined in the range [0, 1].
    """

    # Get the nearest nodes in the graph
    org_node = ox.distance.nearest_nodes(G, orgine[1], orgine[0])
    dest_node = ox.distance.nearest_nodes(G, destination[1], destination[0])

    # Calculate the routes for each alpha value
    routes = {}
    for alpha in alpha_list:
        # Calculate the route
        route = ox.shortest_path(G, org_node, dest_node, weight=f'weight_{alpha}')
        routes[alpha] = route

    return routes

def route_to_list_coordinates(route, G):
    """
    This function takes a list of nodes as the list of their identifiers
    and returns a list of tuples with the coordinates of each node. It
    is important to note that all the nodes identified by their id must
    be in the graph. If not, the function will raise an error.    
    
    For each node, a tuple with the coordinates (lat, lon) is returned.
    
    Parameters
    ----------
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

    return [[G.nodes[node]['y'], G.nodes[node]['x']] for node in route]

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

def add_route_to_map(route, alias, color, map):
    """
    This function takes a list of coordinates and adds
    a polyline to the folium map representing the route.

    Parameters
    ----------
    route : list
        List of coordinates (lat, lon) of the route.
    alias : str
        Alias of the route. It will be used as the tooltip
        of the polyline. Example: "Route 1".
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
    route_layer = folium.FeatureGroup(name=alias)

    # Add the route to the layer
    folium.PolyLine(
        locations=route,
        color=color,
        weight=5,
        opacity=0.7,
        tooltip=alias
    ).add_to(route_layer)

    # Add the layer to the map
    route_layer.add_to(map)

def display_all_routes_on_map(routes_coords, save_html_path):
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
    save_html_path : str
        Path to save the HTML file with the map.

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
        add_route_to_map(route, f"Route {alpha}", list_colors[i], map)

    # Add a layer control to the map
    folium.LayerControl().add_to(map)

    # Save the map to an HTML file
    map.save(save_html_path)
    



