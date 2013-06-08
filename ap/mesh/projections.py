#! /usr/bin/env python
import numpy as np
from math import pi

def lambert_azimuthal(coordinate_tripples, longitude_offset = pi/8,
                      latitude_offset = pi/8):
    """
    An area-preserving projection from geophysical coordinates to the
    plane.

    Required Arguments
    -------------------
    * coordinate_tripples : 3-column, N row tripples of sphere
      (cartesian) coordinates.

    Optional Arguments
    -------------------
    * longitude_offset : offset in radians from Greenwich, England.
      Negative moves towards New York, positive moves towards Moscow.

    Output
    ------
    * the Lambert Azimuthal projection (not scaled) of the
    coordinate_tripples.
    """
    x = coordinate_tripples[0]
    y = coordinate_tripples[1]
    z = coordinate_tripples[2]
    # convert to latitude / longitude coordinates.
    rho = np.sqrt(x**2) + np.sqrt(y**2) + np.sqrt(z**2)
    azimuth = np.arctan2(y,x)                  # same as longitude
    inclination = -1*(np.arccos(z/rho) - pi/2) # same as latitude

    # project to 2D.
    k = np.sqrt(2/(1 + np.cos(inclination - latitude_offset)
                   *np.cos(azimuth - longitude_offset)))

    x_projected = (k*np.cos(inclination - latitude_offset)
                   *np.sin(azimuth - longitude_offset))
    y_projected = k*np.sin(inclination)

    return np.array([x_projected, y_projected])
