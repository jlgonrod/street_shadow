from pyproj import Transformer
import numpy as np
from shapely.geometry import Polygon, MultiPolygon

def convert_coordinates_EPSG_to_4326(x, y, epsg_source):
    """
    This function converts a pair of coordinates from EPSG:***** to EPSG:4326.

    Parameters:
    x (float): The X coordinate in EPSG:25829.
    y (float): The Y coordinate in EPSG:25829.
    epsg (str): The EPSG code of the coordinate system to convert to.
                Example: "EPSG:25829"

    Returns:
    tuple: A tuple with the converted coordinates in EPSG
    """

    # Set up the transformer from EPSG:***** (source) to EPSG:4326 (destination)
    transformer = Transformer.from_crs(epsg_source, "EPSG:4326", always_xy=True)

    # Convert the coordinates
    x_out, y_out = transformer.transform(x, y)

    return x_out, y_out

def convert_coordinates_arrays_EPSG_to_4326(coords_arrays, epsg_source):
    """
    This function converts a list of coordinates from EPSG:***** to EPSG:4326.
    The shaped array should be (2, N) where N is the number of coordinates.

    Parameters:
    coords_arrays (np.array): A numpy array with the coordinates to convert.
    epsg_source (str): The EPSG code of the coordinate system to convert to.
                Example: "EPSG:25829"
    
    Returns:
    np.array: A numpy array with the converted coordinates in EPSG:4326.
    """

    # Set up the transformer from EPSG:25829 (source) to EPSG:4326 (destination)
    transformer = Transformer.from_crs(epsg_source, "EPSG:4326", always_xy=True)

    # Convert the coordinates in the array in the same shape
    coords_out = [transformer.transform(x, y) for x, y in coords_arrays]

    return np.array(coords_out)

def get_square_coords_from_coords(coords_rectangle):
    """
    This function returns the coordinates of a square that contains a rectangle defined by two points.

    Parameters:
    ------------
    coords_rectangle (np.array): A numpy array with the two points that define the rectangle.

    Returns:
    ------------
    np.array: A numpy array with the two points that define the square.
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

def convert_multipolygon_coordinates_EPSG_to_4326(multipolygon, epsg_source):
    """
    Convert the coordinates of a MultiPolygon from EPSG:***** to EPSG:4326.

    Parameters
    ----------
    multipolygon : shapely.geometry.MultiPolygon
        The MultiPolygon to convert.
    epsg_source : str
        The EPSG code of the coordinate system to convert from.
        Example: "EPSG:25829"

    Returns
    -------
    shapely.geometry.MultiPolygon
        The converted MultiPolygon.
    """
    # Set up the transformer from EPSG:25829 (source) to EPSG:4326 (destination)
    transformer = Transformer.from_crs(epsg_source, "EPSG:4326", always_xy=True)
    
    converted_polygons = []
    for polygon in multipolygon.geoms:
        exterior_coords = [transformer.transform(x, y) for x, y in polygon.exterior.coords]
        interiors_coords = [[transformer.transform(x, y) for x, y in interior.coords] for interior in polygon.interiors]
        converted_polygons.append(Polygon(exterior_coords, interiors_coords))
    return MultiPolygon(converted_polygons)