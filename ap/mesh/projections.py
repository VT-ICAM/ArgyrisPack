#! /usr/bin/env python
"""Map projections for use with geophysical data."""
import numpy as np
from numpy import pi


def lambert_azimuthal(coordinate_triples, longitude_offset=pi/8,
                      latitude_offset=pi/8):
    """
    An area-preserving projection from geophysical coordinates to the
    plane.

    Required Arguments
    ------------------
    * coordinate_triples : 3-column, N row triples of sphere
    (cartesian) coordinates.

    Optional Arguments
    ------------------
    * latitude_offset : Offset in radians from the equator. Negative
      moves north, positive moves south.

    * longitude_offset : Offset in radians from Greenwich, England.
      Negative moves towards New York, positive moves towards Moscow.

    Output
    ------
    * Lambert Azimuthal projection (not scaled) of coordinate_triples.

    References
    ----------
    See Wolfram Mathworld

        http://140.177.205.23/LambertAzimuthalEqual-AreaProjection.html

    for details.
    """
    latitudes, longitudes = cartesian_to_geographical(coordinate_triples)
    k = np.sqrt(2/(1 + np.cos(latitudes - latitude_offset)
                   *np.cos(longitudes - longitude_offset)))
    x_projected = (k*np.cos(latitudes - latitude_offset)
                   *np.sin(longitudes - longitude_offset))
    y_projected = k*np.sin(latitudes - latitude_offset)
    return np.array([x_projected, y_projected])


def miller_cylindrical(coordinate_triples, longitude_offset=0,
                       latitude_offset=0):
    """
    A modified Mercador projection from geophysical coordinates to the
    plane.

    Required Arguments
    -------------------
    * coordinate_triples : 3-column, N row triples of sphere
      (cartesian) coordinates.

    Optional Arguments
    ------------------
    * latitude_offset : Offset in radians from the equator. Negative
      moves north, positive moves south.

    * longitude_offset : Offset in radians from Greenwich, England.
      Negative moves towards New York, positive moves towards Moscow.

    Output
    ------
    * Miller cylindrical projection (not scaled) of coordinate_triples.
    """
    latitudes, longitudes = cartesian_to_geographical(coordinate_triples)
    x_projected = longitudes - longitude_offset
    y_projected = 1.25*np.log(np.tan(pi/4 + 0.4*(latitudes - latitude_offset)))
    return np.array([x_projected, y_projected])


def cartesian_to_geographical(coordinate_triples):
    """
    Convert Cartesian coordinates to geophysical (latitude and longitude)
    coordinates.

    Required Arguments
    -------------------
    * coordinate_triples : 3-column, N row triples of sphere
    (cartesian) coordinates.

    Output
    ------
    * The tuple (latitudes, longitudes).
    """
    if len(coordinate_triples.shape) == 1:
        x = coordinate_triples[0]
        y = coordinate_triples[1]
        z = coordinate_triples[2]
    elif len(coordinate_triples.shape) == 2:
        assert coordinate_triples.shape[1] == 3
        x = coordinate_triples[:, 0]
        y = coordinate_triples[:, 1]
        z = coordinate_triples[:, 2]
    radius = np.sqrt(x**2 + y**2 + z**2)
    longitudes = np.arctan2(y, x)
    latitudes = np.arcsin(z/radius)
    return (latitudes, longitudes)
