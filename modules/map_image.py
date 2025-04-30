import io
import os
import re
import requests
import numpy as np
import pyvista as pv
from PIL import Image

from .coordinates import convert_coordinates_arrays_EPSG_to_4326, get_square_coords_from_coords

BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def get_min_max_coords(lon1, lat1, lon2, lat2):
    """
    This function returns the minimum and maximum coordinates of a bounding box
    defined by two points (4 values) that are opposite corners of the box but 
    not necessarily the maximum and minimum coordinates.

    Parameters
    ----------
    lon1 : float
        Longitude of the first point
    lat1 : float
        Latitude of the first point
    lon2 : float
        Longitude of the second point
    lat2 : float
        Latitude of the second point

    Returns
    -------
    lon_min : float
        The minimum longitude of the bounding box
    lat_min : float
        The minimum latitude of the bounding box
    lon_max : float
        The maximum longitude of the bounding box
    lat_max : float
        The maximum latitude of the bounding box
    """
    lon_min = min(lon1, lon2)
    lon_max = max(lon1, lon2)
    lat_min = min(lat1, lat2)
    lat_max = max(lat1, lat2)

    return lon_min, lat_min, lon_max, lat_max

def get_query_static_image(lon1, lat1, lon2, lat2):
    """
    This function returns a query to get a static image from Mapbox API.
    The query is based on a bounding box defined by two points that are opposite.

    Parameters
    ----------
    lon1 : float
        Longitude of the first point
    lat1 : float
        Latitude of the first point
    lon2 : float
        Longitude of the second point
    lat2 : float
        Latitude of the second point

    Returns
    -------
    query : str
        The query to get the static image
    """
    
    # Get the minimum and maximum coordinates of the bounding box
    lon_min, lat_min, lon_max, lat_max = get_min_max_coords(lon1, lat1, lon2, lat2)

    username_style = "mapbox" # The username that owns the style
    token = "<YOUR_MAPBOX_ACCESS_TOKEN>" # The Mapbox access token (https://console.mapbox.com/account/access-tokens/)
    style_id = "streets-v12"
    overlay = ""
    bbox = f"[{lon_min},{lat_min},{lon_max},{lat_max}]"
    width = 400
    height = 400

    query=f"https://api.mapbox.com/styles/v1/{username_style}/{style_id}/static/{overlay}/{bbox}/{width}x{height}?access_token={token}"

    # Remove double slashes excluding the one after the protocol
    query = re.sub(r'(?<!:)//', '/', query)

    return query

def get_image_map(lat1, lon1, lat2, lon2):
    """
    This function returns a static image from Mapbox API based on two points that are opposite corners of a bounding box.

    Parameters
    ----------
    lat1 : float
        Latitude of the first point
    lon1 : float
        Longitude of the first point
    lat2 : float
        Latitude of the second point
    lon2 : float
        Longitude of the second point

    Returns
    -------
    image : bytes
        The static image
    """
    query = get_query_static_image(lon1, lat1, lon2, lat2)
    response = requests.get(query)
    image = response.content

    return image

def parse_image_map(image):
    """
    This function parses the image from Mapbox API and returns the image as a PIL Image.

    Parameters
    ----------
    image : bytes
        The image to parse

    Returns
    -------
    img : PIL Image
        The image as a PIL Image
    """
    img = Image.open(io.BytesIO(image))
    return img

def get_image_from_coords(coords_gml, epsg_source):
    """
    This function plots a static image from Mapbox API based on a list of coordinates
    obtained from a GML file and saves the image in a file.

    Parameters
    ----------
    coords_gml : list
        A list of coordinates obtained from a GML file in  EPSG:25829
        The list should be in the form [[lon1, lat1], [lon2, lat2], ...]

    Returns
    -------
    save_path : str
        The path to the saved image
    """
    # Convert the coordinates to an array and get the minimum and maximum coordinates
    coords_array = np.array(coords_gml)

    min_coords = np.min(coords_array, axis=0)
    max_coords = np.max(coords_array, axis=0)
    
    coords_min_max = np.array([min_coords, max_coords])

    # Convert the coordinates to EPSG:4326
    coords_min_max = convert_coordinates_arrays_EPSG_to_4326(coords_min_max, epsg_source)

    # Reorder from [[lon, lat]] to [[lat, lon]]
    coords_min_max = np.array([[coords_min_max[0][1], coords_min_max[0][0]], [coords_min_max[1][1], coords_min_max[1][0]]])

    # Get the image from the Mapbox API
    map_img = get_image_map(coords_min_max[0][0], coords_min_max[0][1], coords_min_max[1][0], coords_min_max[1][1])
    img = parse_image_map(map_img)

    # Save the image in a temporary file
    save_path = os.path.join(BASE_DIR, "data", "temp", "map_image.png")
    if not os.path.exists(os.path.dirname(save_path)):
        os.makedirs(os.path.dirname(save_path))
    img.save(save_path)

    return save_path

def draw_base_map(list_coords, epsg_source, pad=0):
    """
    This function creates a base plane that will be used to contain the image
    and loads the image as a texture.

    Parameters
    ----------
    list_coords : list
        A list of coordinates that define a polygon in EPSG:25829
        with the form [[lon1, lat1], [lon2, lat2], ...]
    epsg_source (str): The EPSG code of the coordinate system to convert to.
        Example: "EPSG:25829"
    pad : float, optional
        Padding to add to the bounding box, by default 0.

    Returns
    -------
    base_plane : pyvista.PolyData
        The base plane
    texture : pyvista.Texture
        The texture of the map
    center : np.array
        The center of the bounding box
    """
    
    # Get the coordinates from the bounding box
    min_x, min_y = np.min(list_coords, axis=0)
    max_x, max_y = np.max(list_coords, axis=0)

    coords_min_max = np.array([[min_x, min_y], [max_x, max_y]])

    # Add padding to the bounding box
    pad = float(pad) # Convert to float to avoid errors
    coords_min_max[0] -= pad
    coords_min_max[1] += pad

    # The coordinates need to define a square because the image provided 
    # by Mapbox API is a square
    coords_min_max = get_square_coords_from_coords(coords_min_max)

    # Get the center of the bounding box
    center = np.mean(coords_min_max, axis=0)

    # Create the base plane that will be used to contain the image
    base_plane = pv.Plane(
        center=(center[0], center[1], -0.2), #-0.1 to avoid shadows overlapping
        i_size=abs(coords_min_max[0, 0] - coords_min_max[1, 0]),
        j_size=abs(coords_min_max[0, 1] - coords_min_max[1, 1]),
    )
    
    base_plane.texture_map_to_plane(inplace=True)

    # Get the image from the Mapbox API and save it in the temporary path
    map_img_path = get_image_from_coords(coords_min_max, epsg_source)

    # Load the image as a texture
    texture = pv.read_texture(map_img_path)

    # Remove the temporary file
    os.remove(map_img_path)

    return base_plane, texture, center