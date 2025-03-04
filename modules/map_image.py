import requests
import re
from PIL import Image
import io
import numpy as np
from .coordinates import convert_coordinates_arrays

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
    token = "pk.eyJ1Ijoiamxnb25yb2QiLCJhIjoiY203c3hjMmhvMWNvdjJqc2Rudm90OWhpOCJ9.-4kBagshD5htWXB7xadG2A"
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

def get_square_coords_from_coords(coords_rectangle):
    """
    This function returns the coordinates of a square that contains a rectangle defined by two points.
    """

    x1 = coords_rectangle[0][0]
    y1 = coords_rectangle[0][1]

    x2 = coords_rectangle[1][0]
    y2 = coords_rectangle[1][1]

    # Get que center of the square
    x_center = (x1 + x2) / 2
    y_center = (y1 + y2) / 2

    # Get the maximun width or height
    long = max(abs(x1 - x2), abs(y1 - y2))

    # Get the new coordinates
    x1 = x_center - long / 2
    y1 = y_center - long / 2

    x2 = x_center + long / 2
    y2 = y_center + long / 2

    return np.array([[x1, y1], [x2, y2]])

def plot_image_from_coords(coords_gml, pad=0):
    """
    This function plots a static image from Mapbox API based on a list of coordinates
    obtained from a GML file.

    Parameters
    ----------
    coords_gml : list
        A list of coordinates obtained from a GML file in  EPSG:25829
        The list should be in the form [[lon1, lat1], [lon2, lat2], ...]
    pad : float
        The padding to add to the bounding box in meters
    """
    # Convert the coordinates to an array and get the minimum and maximum coordinates
    coords_array = np.array(coords_gml)

    min_coords = np.min(coords_array, axis=0)
    max_coords = np.max(coords_array, axis=0)

    # Add the padding
    min_coords -= pad
    max_coords += pad
    
    coords_min_max = np.array([min_coords, max_coords])

    # Convert the coordinates to EPSG:4326
    coords_min_max = convert_coordinates_arrays(coords_min_max)

    # Reorder from [[lon, lat]] to [[lat, lon]]
    coords_min_max = np.array([[coords_min_max[0][1], coords_min_max[0][0]], [coords_min_max[1][1], coords_min_max[1][0]]])

    # Get the coordinates of a square that contains the building
    coords_min_max = get_square_coords_from_coords(coords_min_max)

    # Get the image from the Mapbox API
    map_img = get_image_map(coords_min_max[0][0], coords_min_max[0][1], coords_min_max[1][0], coords_min_max[1][1])
    img = parse_image_map(map_img)

    # Show the image
    img.show()