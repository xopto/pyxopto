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

from typing import Tuple, List
import numpy as np


Vector3 = Tuple[float, float, float]
SquareMatrix3 = Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]

Vector2 = Tuple[float, float]
SquareMatrix2 = Tuple[Tuple[float, float], Tuple[float, float]]

def rotation_matrix_2d(a: Vector2, b:Vector2):
    '''
    Computes a rotation matrix that rotates a unit vector a
    onto a unit vector b.

    Parameters
    ----------
    a: (float, float)
        A unit vector from which the rotation starts.
    b: (float, float)
        A unit vector onto which the vector a is rotated.

    Returns
    -------
    T: np.ndarray of shape (2, 2)
        Rotation matrix that rotates the unit vector a onto the
        unit vector b (b = T*a)
    '''
    a = np.asarray(a, dtype=np.float)
    b = np.asarray(b, dtype=np.float)
    a = a/np.linalg.norm(a)
    b = b/np.linalg.norm(b)

    cos_theta = np.dot(a, b)
    sin_theta = np.sqrt(1.0 - cos_theta**2)

    return np.array([[cos_theta, -sin_theta],
                     [sin_theta, cos_theta]])


def rotation_matrix(a: Vector3, b: Vector3) -> SquareMatrix3:
    '''
    Computes a rotation matrix that rotates a unit vector a
    onto a unit vector b.
    Parameters
    ----------
    a: (float, float, float)
        A unit vector from which the rotation starts.
    b: (float, float, float)
        A unit vector onto which the vector a is rotated.

    Returns
    -------
    T: np.ndarray of shape (3, 3)
        Rotation matrix that rotates the unit vector a onto the
        unit vector b (b = T*a)

    Note
    ----
    The general solution fails if a == -b (rotation is not uniquely defined).
    '''
    a = np.asarray(a, dtype=np.float)
    b = np.asarray(b, dtype=np.float)
    a = a/np.linalg.norm(a)
    b = b/np.linalg.norm(b)
    v = np.cross(a, b)
    # s = np.linalg.norm(v)
    c = np.dot(a.flat, b.flat)

    vx = np.array([[0.0,  -v[2],  v[1]],
                   [v[2],  0.0,  -v[0]],
                   [-v[1], v[0],  0.0]], dtype=np.float)

    R = np.identity(3) + vx + np.dot(vx, vx)*(1.0/(1.0 + c))

    return R

def transform_base(vfrom: Vector3, vto: Vector3) -> SquareMatrix3:
    '''
    Create a matrix that transforms the orthogonal base defining
    the components of vectors vfrom and vto by rotating vector vfrom onto
    vector vto.
    Use the transformation matrix to transform coordinates defined in the
    transformed orthogonal base back to the initial orthogonal base that
    was used to define the components of vectors vfrom and vto. 

    Parameters
    ----------
    from_v: (float, float, float)
        Initial vector defined with the source orthogonal base.
    to_v: (float, float, float)
        Target vector defined with the source orthogonal base.

    Returns
    -------
    transform: np.ndarray of shape (3, 3)
        Transformation matrix.

    Examples
    --------
    Create a transformation from coordinate system with a standard
    z axis [0, 0, 1] to a coordinate system with a z axis [sin(10), 0, cos(10)].
    This requires additional step that inverts the matrix returned by
    transform_base.

        >>> Tinv = np.linalg.inv(transform_base([0, 0, 1], [np.sin(np.pi/18), 0, np.cos(np.pi/18)]))
        >>> np.dot(Tinv, [[0], [0], [1]])
    '''
    t_from_to = rotation_matrix(vfrom, vto)

    return t_from_to
