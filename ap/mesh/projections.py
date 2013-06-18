#! /usr/bin/env python
import numpy as np
from math import pi

def lambert_lonal(coordinate_triples, longitude_offset = pi/8.,
                      latitude_offset = pi/8.):
    """
    An area-preserving projection from geophysical coordinates to the
    plane.

    Required Arguments
    -------------------
    * coordinate_triples : 3-column, N row triples of sphere
      (cartesian) coordinates.

    Optional Arguments
    -------------------
    * longitude_offset : offset in radians from Greenwich, England.
      Negative moves towards New York, positive moves towards Moscow.

    Output
    ------
    * the Lambert Azimuthal projection (not scaled) of the
    coordinate_triples.
    """
    lat, lon = convert_to_latlon(coordinate_triples)

    # project to 2D.
    k = np.sqrt(2./(1. + np.cos(lat - latitude_offset)
                   *np.cos(lon - longitude_offset)))

    x_projected = (k*np.cos(lat - latitude_offset)
                   *np.sin(lon - longitude_offset))
    y_projected = k*np.sin(lat)

    return np.array([x_projected, y_projected])

def miller_cylindrical(coordinate_triples):
    """
    A modified mercador projection from geophysical coordinates to the
    plane.

    Required Arguments
    -------------------
    * coordinate_triples : 3-column, N row triples of sphere
      (cartesian) coordinates.

    Output
    ------
    * the Miller cylindrical projection (not scaled) of the
    coordinate_triples.
    """

    lat, lon = convert_to_latlon(coordinate_triples) 
    lon_0 = 0.                           # center meridian
    #now use miller cylindrical projection
    x_projected = lon - lon_0
    y_projected = 5./4.*np.log(np.tan(pi/4.+2./5.*lat))

    return np.array([x_projected, y_projected])

def convert_to_latlon(coordinate_triples): 
    """
    converts cartesian coordinates (x,y,z) to lat,lon

    Required Arguments
    -------------------
    * coordinate_triples : 3-column, N row triples of sphere
      (cartesian) coordinates.

    Output
    ------
    * the corresponding latitude and longitude for a given point x,y,z
    """
    x = coordinate_triples[0]
    y = coordinate_triples[1]
    z = coordinate_triples[2]

    # convert to latitude / longitude coordinates.
    R = np.sqrt(x**2 + y**2 + z**2)           # Radius of the Earth
    lon = np.arctan2(y,x)                     # longitude
    lat = np.arcsin(z/R)                      # latitude

    return np.array([lat, lon])
