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

import numpy as np

def cos_critical(n1: float, n2: float) -> float:
    '''
    Cosine of the critical angle of incidence beyond which the incident
    beam is reflected at the boundary n1 => n2.

    Parameters
    ----------
    n1: float
        Refractive index of the material on the incident side of the 
        boundary.
    n2: float
        Refractive index of the material across the boundary.
    '''
    cc = 0.0
    if n1 > n2:
        cc = (1.0 - (n2/n1)**2)**0.5
    return cc

def reflectance(
            n1: float, n2: float, costheta: float=1.0) -> float:
    '''
    Computes reflectance of unpolarized light at the specified boundary.

    Parameters
    ----------
    n1: float
        Refractive index of the material on the side of the incident beam.
    n2: float
        Refractive index of the material across the boundary.
    costheta: float
        Cosine of the angle of incidence -
        perpendicular incidence by default (90 deg, cosine is 1).

    Returns
    -------
    Returns the reflectance of unpolarized light.
    '''
    if costheta < 0.0:
        raise ValueError('The incidence angle cosine must not be negative!')

    n1, n2 = float(n1), float(n2)
    sintheta = (1.0 - costheta**2)**0.5
    sincritical = n2/n1

    if n1 > n2 and sintheta >= sincritical:
        Rs = Rp = 1.0
    else:
        a1, a2 = n1*costheta, n2*(1.0 - (n1/n2*sintheta)**2)**0.5
        Rs = np.abs((a1 - a2)/(a1 + a2))**2

        b1, b2 = n1*(1.0 - (n1/n2*sintheta)**2)**0.5, n2*costheta
        Rp = np.abs((b1 - b2)/(b1 + b2))**2

    return (Rs + Rp)*0.5

def refract(direction: np.ndarray, normal: np.ndarray, n1: float, n2: float) \
        -> np.ndarray:
    '''
    Refract the beam across the given boundary with refractive indices
    n1 and n2.

    Parameters
    ----------
    direction: np.ndarray
        Propagation direction of the incident beam.
    normal: np.ndarray
        Boundary surface normal (pointing inwards or outwards).
    n1: float
        Refractive index on the incident side of the medium.
    n2: float
        Refractive index across the boundary.

    Returns
    -------
    direction: np.ndarray
        Propagation direction of the refracted beam.
    '''
    direction = np.asarray(direction, dtype=np.float)
    direction_len = np.linalg.norm(direction)

    normal = np.asarray(normal, dtype=np.float)
    normal_len = np.linalg.norm(normal)

    # For outwards pointing normal, cos1 is negative.
    cos1 = np.dot(direction, normal)/(direction_len*normal_len)
    if abs(cos1) < cos_critical(n1, n2):
        raise ValueError('Cannot refract, the incidence angle '
                         'exceeds the critical angle!') 
    n1_d_n2 = n1/n2
    sin2_squared = n1_d_n2*n1_d_n2*(1.0 - cos1*cos1)
    k = n1_d_n2*abs(cos1) - np.sqrt(1.0 - sin2_squared)

    return n1_d_n2*direction/direction_len - np.sign(cos1)*k*normal/normal_len

def reflect(direction: np.ndarray, normal: np.ndarray) -> np.ndarray:
    '''
    Reflect the beam direction from a boundary with the given normal.

    Parameters
    ----------
    direction: np.ndarray
        Propagation direction of the incident beam.
    normal: np.ndarray
        Boundary surface normal (pointing inwards or outwards).

    Returns
    -------
    reflected_dir: np.ndarray
        Propagation direction of the reflected beam.
    '''
    direction = np.asarray(direction, dtype=np.float)
    direction_len = np.linalg.norm(direction)
    if direction_len > 0.0:
        direction = direction/direction_len

    normal = np.asarray(normal, dtype=np.float)
    normal_len = np.linalg.norm(normal)
    if normal_len > 0.0:
        normal = normal/normal_len

    k = 2.0*np.dot(direction, normal)

    return direction - k*normal
