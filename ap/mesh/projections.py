#! /usr/bin/env python
"""Map projections for use with geophysical data."""
import numpy as np

def lambert_azimuthal(coordinate_tripples, longitude_offset = np.pi/8,
                      latitude_offset = np.pi/8):
    """
    An area-preserving projection from geophysical coordinates to the
    plane.

    Required Arguments
    ------------------
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

    References
    ----------
    See Wolfram Mathworld

        http://140.177.205.23/LambertAzimuthalEqual-AreaProjection.html

    for details.
    """
    x = coordinate_tripples[0]
    y = coordinate_tripples[1]
    z = coordinate_tripples[2]
    # convert to latitude / longitude coordinates.
    rho = np.sqrt(x**2 + y**2 + z**2)
    azimuth = np.arctan2(y, x)                    # same as longitude
    inclination = -1*(np.arccos(z/rho) - np.pi/2) # same as latitude

    # project to 2D.
    k = np.sqrt(2/(1 + np.cos(inclination - latitude_offset)
                   *np.cos(azimuth - longitude_offset)))

    x_projected = (k*np.cos(inclination - latitude_offset)
                   *np.sin(azimuth - longitude_offset))
    y_projected = k*np.sin(inclination)

    return np.array([x_projected, y_projected])
