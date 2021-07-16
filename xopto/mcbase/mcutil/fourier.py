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

# Import the Symmetric Fourier Transform class
from scipy.interpolate import interp1d
from scipy.integrate import quad, simps

def _uneven(array):
    tmp = np.diff(array)
    return np.any(np.abs(tmp - tmp[0]) > np.finfo(np.float).eps)

def discreteSimpson(frequency: list or tuple or np.ndarray,
                    xpts: np.ndarray, fpts: np.ndarray, uneven: bool = False) \
                        -> np.ndarray:
    '''
    Computes Fourier transform of a 1D function defined on a grid of evenly or
    unevenly spaced points. To compute transforms of multiple sets (functions),
    the fpts array shape must be (num_xfuns, xpts.size).

    Parameters
    ----------
    frequency: list, tuple, ndarray vector
        A list of frequencies at which to compute the Fourier transform.
    xpts: np.ndarray vector
        A vector of evenly or unevenly spaced points at which the function
        values in fpts are defined.
    fpts: np.ndarray vector of 2D array
        A vector or array of function values defined at points xpts. To compute
        transforms of multiple sets (functions), the fpts array shape must be
        (num_xfuns, xpts.size).
    uneven: bool
        If True, the method assumes unevenly spaced values in xpts. Default is
        False. If set to None, the value of uneven flag is derived from
        the values in the xpts array.

    g(q) = 2*pi*int_0^inf(f(r)*J0(2*pi*q*r)*r*dr)
    Returns
    -------
    F: np.ndarray vector
        The Fourier transfor of xfun at the given frequencies. If the fpts array
        ia a vector (points of one function only) then F is a vector of size
        len(frequencies). If fpts is a 2D array of shape (N, xpts.size) then
        F is a 2D array of shape (N, len(frequencies)).
    '''
    np_freqs = np.asarray(frequency)

    if fpts.ndim > 1:
        out = np.empty((fpts.shape[0], np_freqs.size,), dtype=np.complex)
        xpts = np.reshape(xpts, (1, xpts.size))
    else:
        out = np.empty((np_freqs.size,), dtype=np.complex)

    if uneven is None:
        uneven = _uneven(xpts)

    x = dx = None
    if uneven:
        x = xpts
        if out.ndim > 1:
            x = np.reshape(x, (1, x.size))
    else:
        dx = xpts.flat[1] - xpts.flat[0]

    for index in range(np_freqs.size):
        f = fpts*np.exp(-2.0*np.pi*1j*xpts*np_freqs[index])
        if out.ndim > 1:
            out[:, index] = simps(f, x, dx=dx)
        else:
            out[index] = simps(f, x, dx=dx)

    return out
