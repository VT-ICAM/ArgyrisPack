#! /usr/bin/env python

from MeshParser import ParseMESHFormat
import numpy as np
from math import pi

def main():
    # mesh_file = "earth.mesh"
    # t = ParseMESHFormat(mesh_file, projection = lambert_azimuthal_projection,
    #                     special_borders = {'open' : (306, 374)})
    # t.save_QGE_outfiles()

    mesh_file = "unitsquare.mesh"
    t = ParseMESHFormat(mesh_file, special_borders = {'ocean' : (1,2)})
    t.save_outfiles()

def lambert_azimuthal_projection(coordinate_tripples, longitude_offset = -pi/4):
    """
    An area-preserving projection from geophysical coordinates to the plane.

    Required Arguements
    -------------------
    * coordinate_tripples : 3-column, N row tripples of sphere coordinates.

    Optional Arguements
    -------------------
    * longitude_offset : offset in Radians from Greenwich, England. Negative
      moves towards New York, positive moves towards Moscow.

    Output
    ------
    * the Lambert Azimuthal projection (not scaled) of the coordinate_tripples.
    """
    x = coordinate_tripples[:,0]
    y = coordinate_tripples[:,1]
    z = coordinate_tripples[:,2]
    # convert to latitude / longitude coordinates.
    rho = np.sqrt(x**2) + np.sqrt(y**2) + np.sqrt(z**2)
    azimuth = np.arctan2(y,x)                    # same as longitude
    inclination = -1 * (np.arccos(z/rho) - pi/2) # same as latitude

    # project to 2D.
    k = np.sqrt(2 / (1 + np.cos(inclination) * np.cos(azimuth - longitude_offset)))

    x_projected = (k * np.cos(inclination) * np.sin(azimuth - longitude_offset))
    y_projected = (k * np.sin(inclination))

    return np.hstack((x_projected.reshape((len(x),1)),
                      y_projected.reshape((len(x),1))))

if __name__ == '__main__':
    main()
