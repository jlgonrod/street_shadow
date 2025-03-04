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