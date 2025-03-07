from pyproj import Transformer
import numpy as np

def convert_coordinates(x, y):
    """
    This function converts a pair of coordinates from EPSG:25829 to EPSG:4326.

    Parameters:
    x (float): The X coordinate in EPSG:25829.
    y (float): The Y coordinate in EPSG:25829.

    Returns:
    tuple: A tuple with the converted coordinates in EPSG
    """
    # Set up the transformer from EPSG:25829 (source) to EPSG:4326 (destination)
    transformer = Transformer.from_crs("EPSG:25829", "EPSG:4326", always_xy=True)

    # Convert the coordinates
    x_out, y_out = transformer.transform(x, y)

    return x_out, y_out

def convert_coordinates_arrays(coords_arrays):
    """
    This function converts a list of coordinates from EPSG:25829 to EPSG:4326.
    The shaped array should be (2, N) where N is the number of coordinates.
    """

    # Set up the transformer from EPSG:25829 (source) to EPSG:4326 (destination)
    transformer = Transformer.from_crs("EPSG:25829", "EPSG:4326", always_xy=True)

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