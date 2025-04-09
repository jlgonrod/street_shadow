from pvlib.solarposition import get_solarposition
import pandas as pd
from zoneinfo import ZoneInfo
import numpy as np
from .coordinates import convert_coordinates_EPSG_to_4326

def get_zenith_azimuth(lat, lon, datetime):
    """
    Get the zenith and azimuth angles of the sun given the latitude,
    longitude and datetime.
    
    Parameters
    ----------
    lat : float
        Latitude of the location.
    lon : float
        Longitude of the location.
    datetime : pandas.Timestamp
        Datetime of the calculation.
        Indicate the timezone with tz=ZoneInfo({timezone})."""
    solarposition = get_solarposition(datetime, lat, lon, altitude=0)
    return solarposition['zenith'].values, solarposition['azimuth'].values

def light_angles2vec(azimuth_ang, zenith_angle):
    """
    Calculate the unit vector of light pointing to the origin of the system
    given the azimuth and zenith angles in degrees.

    Parameters
    ----------
    azimuth_ang : float
        Azimuth angle in degrees.
    zenith_angle : float
        Zenith angle in degrees.

    Returns
    -------
    numpy.ndarray
        Unit vector of light pointing to the origin of the system.
    """
    azimuth = np.radians(azimuth_ang)
    zenith = np.radians(zenith_angle)
    sin_zenith = np.sin(zenith)
    cos_zenith = np.cos(zenith)
    sin_azimuth = np.sin(azimuth)
    cos_azimuth = np.cos(azimuth)
    
    x = -sin_zenith * sin_azimuth
    y = -sin_zenith * cos_azimuth
    z = -cos_zenith
    
    return np.round(np.array([x, y, z]), 3)

def get_sulight_vector(x, y, dt, epsg_source, convert_coords):
    """
    Get the light vector pointing to the origin of the system given the
    coordinates in EPSG:25829 or other (set convert_coords= True)
    If the coordinates are in EPSG:4326 (set convert_coords= False).
    The location, and the datetime are also needed. Assumes z=0.

    Parameters
    ----------
    x : float
        Longitude of the location. 
    y : float
        Latitude of the location.
    dt : pandas.Timestamp
        Datetime of the calculation.
        Indicate the timezone with tz=ZoneInfo({timezone}).
    epsg_source : str
        EPSG code of the source coordinate system.
        It should be 'EPSG:25829' or 'EPSG:4326'.
    convert_coords : bool
        If True, convert coordinates from EPSG:xxxxx to EPSG:4326.
        If False, use the coordinates as they are in EPSG:4326.

    Returns
    -------
    numpy.ndarray
        Unit vector of light pointing to the origin of the system.
    """
    
    # Convert the coordinates from EPSG:25829 to EPSG:4326
    if convert_coords:
        # Convert coordinates from EPSG:25829 to EPSG:4326
        lon, lat = convert_coordinates_EPSG_to_4326(x, y, epsg_source)
    else:
        # Use the coordinates as they are in EPSG:4326
        lon = x
        lat = y

    # Get the zenith and azimuth angles (degrees)
    zenith, azimuth = get_zenith_azimuth(lat, lon, dt)

    # Calculate the light vector
    return light_angles2vec(azimuth, zenith).squeeze()


if __name__ == '__main__':
    
    lat = 36.71964
    lon = -4.41993
    dt = pd.Timestamp(
        year=2025,
        month=3,
        day=11,
        hour=19,
        minute=21,
        second=51,
        tz=ZoneInfo('Europe/Madrid') # It cares about the change of hour in summer and winter in Spain
    )

    print(f'Latitude: {lat}\nLongitude: {lon}\nDate: {dt}')
    zenith, azimuth = get_zenith_azimuth(lat, lon, dt)
  
    zenith = zenith.values[0]
    azimuth = azimuth.values[0]
    altitude = 90 - zenith
  
    print(f'\nAzimuth: {round(azimuth, 2)}')
    print(f'Zenith: {round(zenith, 2)}')
    print(f'Altitude: {round(altitude, 2)}')

    light_vector = light_angles2vec(azimuth, zenith)
    print(f'\nLight vector: {light_vector}')