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

from scipy.ndimage import map_coordinates
import numpy as np


def interp1(xi: np.ndarray, x: np.ndarray, f: np.ndarray, order: int = 1,
            **kwargs) -> np.ndarray:
    '''
    Interpolates 1D functions that are defined on a uniform grid.

    Parameters
    ----------
    xi: np.ndarray
        Points at which to interpolate f.
    x: np.ndarray
        A vector of uniformly distributed points at which the values
        of f are defined.
    f: np.ndarray
        Function values at x. The shape of f should be (x.size,).
    order: int
        Interpolation order. 0 - nearest neighbour, 1 - linear, ...
    kwargs: dict
        Keyword arguments passed to map_coordinates.

    Returns
    -------
    fxi : np.ndarray
        Interpolated values of f at xi.

    Example
    -------
    >>> import numpy as np
    >>> from matplotlib import pyplot as pp
    >>>
    >>> x = np.linspace(0, np.pi, 5)
    >>> x_ref = np.linspace(x[0], x[-1], 1000)
    >>> f = np.cos(x)
    >>>
    >>> xi = np.linspace(0, np.pi, 50)
    >>> fi_lin = interp1(xi, x, f)
    >>> fi_quad = interp1(xi, x, f, 2)
    >>>
    >>> pp.figure()
    >>> pp.plot(x_ref, np.cos(x_ref), '-k', label='cos(x)')
    >>> pp.plot(x, f, 'or', label='sample points', markersize=6)
    >>> pp.plot(xi, fi_lin, 'xg', label='linear', markersize=6)
    >>> pp.plot(xi, fi_quad, 'xb', label='quadratic', markersize=6)
    >>> pp.legend()
    '''
    if x is not None:
        x = np.asarray(x).flatten()
        if x.size != f.size:
            raise IndexError('Length of vector x and f must be the same.')
        indx = (xi - x[0])*((x.size - 1)/(x[-1] - x[0]))
    else:
        indx = xi

    return map_coordinates(
        np.asarray(f), np.asarray([indx]), order=order, **kwargs)

def interp2(xi: np.ndarray, yi: np.ndarray,
            x: np.ndarray, y: np.ndarray, f: np.ndarray, order: int = 1,
            **kwargs) -> np.ndarray:
    '''
    Interpolates 2D functions that are defined on a uniform grid.

    Parameters
    ---------
    xi: np.ndarray
        X coordinates of the points at which to interpolate f.
    yi: np.ndarray
        Y coordinates of the points at which to interpolate f.
    x: np.ndarray
        A vector of uniformly distributed points along the x axis at which
        the values of f are defined.
    y: np.ndarray
        A vector of uniformly distributed points along the y axis at which
        the values of f are defined.
    f: np.ndarray
        Function values at points (x, y).
        The shape of f should be (y.size, x.size).
    order: int
        Interpolation order. 0 - nearest neighbour, 1 - linear, ...
    kwargs: dict
        Keyword arguments passed to map_coordinates.

    Returns
    -------
    fxyi: np.ndarray
        Interpolated values of f at (xi, yi).

    Example
    -------
    >>> import numpy as np
    >>> from matplotlib import pyplot as pp
    >>> from mpl_toolkits.mplot3d import Axes3D
    >>>
    >>> x = np.linspace(-1, 1, 10)
    >>> y = np.linspace(-1, 1, 10)
    >>> Y, X = np.meshgrid(y, x, indexing='ij')
    >>> f = 1.0/(X**2 + Y**2 + 1)
    >>>
    >>> xi = np.linspace(0, 1, 30)
    >>> yi = np.linspace(0, 1, 30)
    >>> Yi, Xi = np.meshgrid(yi, xi, indexing='ij')
    >>> fi = interp2(Xi, Yi, x, y, f)
    >>>
    >>> fig = pp.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> ax.plot_wireframe(X, Y, f, color='r', label='sample points')
    >>> ax.plot_wireframe(Xi, Yi, fi, color='g', label='interpolated values')
    >>> ax.legend()
    '''

    f = np.asarray(f)
    xi = np.array(xi, ndmin=1, copy=False)
    yi = np.array(yi, ndmin=1, copy=False)

    if x is not None:
        x = np.asarray(x)
        if x.size != f.shape[1]:
            raise IndexError('Length of vector x must match ' \
                'the number of columns of f.')
        indx = (xi - x[0])*((x.size - 1)/(x[-1] - x[0]))
    else:
        indx = xi

    if y is not None:
        y = np.asarray(y)
        if y.size != f.shape[0]:
            raise IndexError('Length of vector y must match ' \
                'the number of rows of f.')
        indy = (yi - y[0])*((y.size - 1)/(y[-1] - y[0]))
    else:
        indy = yi

    return map_coordinates(f, np.asarray([indy, indx]), order=order)

def interp3(xi: np.ndarray, yi: np.ndarray, zi: np.ndarray,
            x: np.ndarray, y: np.ndarray, z: np.ndarray, f: np.ndarray,
            order: int = 1, **kwargs) -> np.ndarray:
    '''
    Interpolates 3D functions that are defined on a uniform grid.

    Parameters
    ---------
    xi: np.ndarray
        X coordinates of the points at which to interpolate f.
    yi: np.ndarray
        Y coordinates of the points at which to interpolate f.
    zi: np.ndarray
        Z coordinates of the points at which to interpolate f.
    x: np.ndarray
        A vector of uniformly distributed points along the x axis at which
        the values of f are defined.
    y: np.ndarray
        A vector of uniformly distributed points along the y axis at which
        the values of f are defined.
    z: np.ndarray
        A vector of uniformly distributed points along the z axis at which
        the values of f are defined.
    f: np.ndarray
        Function values at points (x, y, z).
        The shape of f should be (z.size, y.size, x.size).
    order: int
        Interpolation order. 0 - nearest neighbour, 1 - linear, ...
    kwargs: dict
        Keyword arguments passed to map_coordinates.

    Returns
    -------
    fxyzi: np.ndarray
        Interpolated values of f at (xi, yi, zi).

    Example
    -------
    >>> import numpy as np
    >>>>
    >>> x = np.linspace(-1, 1, 10)
    >>> y = np.linspace(-1, 1, 10)
    >>> z = np.linspace(-1, 1, 10)
    >>> Z, Y, X = np.meshgrid(z, y, x, indexing='ij')
    >>> f = 1.0/(X**2 + Y**2 + Z**2 + 1)
    >>>
    >>> xi = np.linspace(0, 1, 30)
    >>> yi = np.linspace(0, 1, 30)
    >>> zi = np.linspace(0, 1, 30)
    >>> Zi, Yi, Xi = np.meshgrid(zi, yi, xi, indexing='ij')
    >>>
    >>> fi = interp3(Xi, Yi, Zi, x, y, z, f)
    '''
    f = np.asarray(f)
    xi, yi = np.asarray(xi), np.asarray(yi)

    if x is not None:
        x = np.asarray(x)
        if x.size != f.shape[2]:
            raise IndexError('Length of vector x must match ' \
                'the number of columns in f.')
        indx = (xi - x[0])*((x.size - 1)/(x[-1] - x[0]))
    else:
        indx = xi

    if y is not None:
        y = np.asarray(y)
        if y.size != f.shape[1]:
            raise IndexError('Length of vector y must match ' \
                'the number of rows in f.')
        indy = (yi - y[0])*((y.size - 1)/(y[-1] - y[0]))
    else:
        indy = yi

    if z is not None:
        z = np.asarray(z)
        if z.size != f.shape[0]:
            raise IndexError('Length of vector z must match ' \
                'the number of slices in f.')
        indz = (zi - z[0])*((z.size - 1)/(z[-1] - z[0]))
    else:
        indy = zi

    return map_coordinates(
        f, np.asarray([indz, indy, indx]), order=order, **kwargs)

def interpn(ti: list, t: list, f: np.ndarray, order: int = 1,
            **kwargs) -> np.ndarray:
    '''
    Interpolates ND functions that are defined on a uniform grid.

    Parameters
    ---------
    ti: list
        A list of coordinates at which to interpolate f ([zi, yi, xi, ...]).
    t: list
        List of uniformly distributed points along each axis of f.
        (z, y, x, ..)
    f: np.ndarray
        Function values at points (x, y, z, ...).
        The shape of f should be (z.size, y.size, x.size, ...).
    order: int
        Interpolation order. 0 - nearest neighbour, 1 - linear, ...
    kwargs: dict
        Keyword arguments passed to map_coordinates.

    Returns
    -------
    fti: np.ndarray
        Interpolated values of f at f[t[0], t[1], ..., t[-1]].

    Example
    -------
    >>> import numpy as np
    >>> from matplotlib import pyplot as pp
    >>> from mpl_toolkits.mplot3d import Axes3D
    >>>
    >>> x = np.linspace(-1, 1, 15)
    >>> y = np.linspace(-1, 1, 10)
    >>> Y, X = np.meshgrid(y, x, indexing='ij')
    >>> f = 1.0/(X**2 + Y**2 + 1)
    >>>
    >>> xi = np.linspace(0, 1, 30)
    >>> yi = np.linspace(0, 1, 30)
    >>> Yi, Xi = np.meshgrid(yi, xi, indexing='ij')
    >>> fi = interpn([Yi, Xi], [y, x], f)
    >>>
    >>> fig = pp.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> ax.plot_wireframe(X, Y, f, color='r', label='sample points')
    >>> ax.plot_wireframe(Xi, Yi, fi, color='g', label='interpolated values')
    >>> ax.legend()
    '''

    f = np.asarray(f)

    if t is not None:
        N = len(t) # space dimensionality
        if N != len(ti):
            raise IndexError('Dimensions of the coordinates must agree!')

        tind = []
        for i in range(N):
            if t[i] is not None:
                tmp = t[i].flatten()
                tind.append((ti[i] - tmp[0])*
                            ((tmp.size - 1)/(tmp[-1] - tmp[0])))
            else:
                tind.append(ti[i])
    else:
        tind = ti
    return map_coordinates(f, np.asarray(tind), order=order, **kwargs)
