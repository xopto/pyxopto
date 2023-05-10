# -*- coding: utf-8 -*-
################################ Begin license #################################
# Copyright (C) Laboratory of Imaging technologies,
#               Faculty of Electrical Engineering,
#               University of Ljubljana.
#
# This file is part of PyXOpto.
#
# PyXOpto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyXOpto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
################################# End license ##################################

from typing import Tuple, List, Type

import numpy as np

def detector_direction(incidence: float, tilt: float, design_angle: float,
                       verbose: bool = False) -> Tuple[float, float, float]:
    '''
    Determines the direction vector of the detector by solving a quadratic
    equation that arises from the following system of equations:

    .. math::

        x\\sin(\\Theta_i) + z\\cos(\\Theta_i) = -\\cos(\\Theta_0)

        y\\cos(\\Theta_t) - z\\sin(\\Theta_t) = 0

        x^2 + y^2 + z^2 = 1

    In the above equations :math:`\Theta_i` is the angle of incidence,
    :math:`\\Theta_t` the detector tilt angle and :math:`\\Theta_0` the
    design angle, i.e. the angle between the source and detector.

    The first equation is the dot product of the detector direction vector
    and the direction of the incident beam. The angle between the two is
    the design angle :math:`180 - \\Theta_0`. The second equation is
    the dot product of the detector direction vector and the normal
    of the tilted plane that contains the detector direction vector. Note
    that the tilted plane is the x-z plane rotated around the x axis.

    The solution with negative z component that has a polar angle
    in the x-z plane smaller than the polar angle of the incident source
    is taken (detector always rotates from source towards the x axis).

    Parameters
    ----------
    incidence: float
        Source incidence angle (radians) in the x-z plane, measured from
        the z axis.
    tilt: float
        Detector tilt angle (radians) measured as the angle between the
        tilted plane (x-z plane rotated around the x axis) and the z axis. 
    design_angle: float
        Design angle between the source and detector (radians).
    verbose: bool
        Enables verbose report.

    Returns
    -------
    dir: np.ndarray
        Direction vector of the detector.
    '''
    # special simplified case for perpendicular incidence
    if incidence == 0:
        z = -np.cos(design_angle)
        y = z*np.sin(tilt)/np.cos(tilt)
        x = np.sqrt(1.0 - y*y - z*z)
        return np.array([x, y, z])

    src_direction = np.array([np.sin(incidence), 0.0, np.cos(incidence)])

    a = 1.0/np.tan(incidence)**2 + np.tan(tilt)**2 + 1.0
    b = 2.0*np.cos(design_angle)*np.cos(incidence)/np.sin(incidence)**2
    c = np.cos(design_angle)**2/np.sin(incidence)**2 - 1.0
    d = (b**2 - 4.0*a*c)**0.5

    z1, z2 = (-b + d)/(2.0*a), (-b - d)/(2.0*a)
    y1, y2 = z1*np.tan(tilt), z2*np.tan(tilt)
    x1 = (-np.cos(design_angle) - z1*np.cos(incidence))/np.sin(incidence)
    x2 = (-np.cos(design_angle) - z2*np.cos(incidence))/np.sin(incidence)

    if verbose:
        print('dir 1             :', x1, y1, z1,)
        print('  src-detector (°):', np.rad2deg(
            np.arccos(-np.dot([x1, y1, z1], src_direction))))
        print('  polar in x-z (°):', np.rad2deg(np.arccos(x1)))
        print('dir 2             :', x2, y2, z2)
        print('  src-detector (°):', np.rad2deg(
            np.arccos(-np.dot([x2, y2, z2], src_direction))))
        print('  polar in x-z (°):', np.rad2deg(np.arccos(x2)))

    if z1 >= 0.0 and z2 >= 0.0:
        raise ValueError('Detector direction could not be resolved!')

    if x1 > np.cos(np.pi*0.5 + incidence):
        direction = np.array([float(x1), float(y1), float(z1)])
    else:
        direction = np.array([float(x2), float(y2), float(z2)])

    return direction
