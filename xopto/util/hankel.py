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

from scipy.interpolate import interp1d
from scipy.integrate import quad, simps
from scipy.special import j0


def _is_uneven(array: np.ndarray) -> bool:
    '''
    Returns True, if the points in the array are unevenly spaced.
    '''
    tmp = np.diff(array)
    return np.any(np.abs(tmp - tmp[0]) > np.finfo(np.float).eps)

def continuous(frequency: np.ndarray, rfun: np.ndarray, rstop: float,
               out: np.ndarray = None, **kwargs) -> np.ndarray:
    '''
    Computes Hankel transform of a continuous radially symmetric function:

    .. math::

        g(q) &= 2 \\pi \\int_0^{\\infty}f(r)J_0(2 \\pi q r) r dr

    Parameters
    ----------
    frequency: np.ndarray
        A list of frequencies (1/m) at which to compute the Hankel transform.
    rfun: callable
        Callable with one parameter (radius) representing a radially symmetric
        function.
    rstop: float
        Range of numerical integration as [0, rstop].
    out: np.ndarray
        Optional output array for the computed frequencies.
    kwargs: dict
        Optional keyword arguments passed to the
        :py:func:`scipy.integrate.quad` function.

    Returns
    -------
    F: np.ndarray vector
        The Hankel transfor of rfun at the given frequencies.
    '''
    np_freq = np.asarray(frequency)

    if out is None:
        out = np.empty([np_freq.size], dtype=np.float)

    for index in range(np_freq.size):
        out[index] = 2*np.pi*quad(
            lambda r, index=index: rfun(r)*j0(2*np.pi*r*np_freq[index])*r,
            0, rstop, **kwargs)[0]
    return out

def discrete(frequency: np.ndarray, rpts: np.ndarray, fpts: np.ndarray,
             logscale: bool = True, **kwargs) -> np.ndarray:
    '''
    Computes Hankel transform of a discrete function. The discrete function is
    first made continuous by means of interpolation in linear or log scale.
    Finally, the transform is computed by the :py:func:`continuous` function
    using quad.

    .. math::

        g(q) &= 2 \\pi \\int_0^{\\infty}f(r)J_0(2 \\pi q r) r dr

    Parameters
    ----------
    frequency: np.ndarray
        A list of frequencies (1/m) at which to compute the Hankel transform.
    rpts: np.ndarray
        Points at which the discrete radially symmetric function is defined.
    fpts: np.ndarray
        Value of the radially symmetric function at points rpts.
    kwargs: dict
        Optional keyword arguments passed to the
        :py:func:`scipy.integrate.quad` function.

    Returns
    -------
    F: np.ndarray
        The Hankel transfor of rfun at the given frequencies.
    '''
    if logscale:
        fpts[fpts <= 0.0] = np.finfo(np.float).eps
        fr = lambda r: np.exp(
            interp1d(rpts, np.log(fpts), assume_sorted=True,
                     bounds_error=False, fill_value='extrapolate')(r)
            )
    else:
        fr = interp1d(rpts, fpts, assume_sorted=True,
                      bounds_error=False, fill_value='extrapolate')

    return continuous(frequency, fr, rpts[-1], **kwargs)

def discrete_simpson(frequency: np.ndarray, rpts: np.ndarray, fpts: np.ndarray,
                     uneven: bool = None) -> np.ndarray:
    '''
    Computes Hankel transform of a radially symmetric function defined
    on a grid of evenly or unevenly spaced points. To compute transforms of
    multiple sets (functions), the fpts array shape must be
    (num_sets, rpts.size).

    .. math::

        g(q) &= 2 \\pi \\int_0^{\\infty}f(r)J_0(2 \\pi q r) r dr

    Parameters
    ----------
    frequency: np.ndarray
        A list of frequencies (1/m) at which to compute the Hankel transform.
    rpts: np.ndarray
        A vector of evenly or unevenly spaced points at which the function
        values in fpts are defined.
    fpts: np.ndarray
        A vector or array of function values defined at points rpts. To compute
        transforms of multiple sets (functions), the fpts array shape must be
        (num_sets, rpts.size).
    uneven: bool
        If True, the method assumes unevenly spaced values in rpts. Default is
        False. If set to None, the value of uneven flag is derived from
        the values in the rpts array.

    Returns
    -------
    F: np.ndarray
        The Hankel transfor of rfun at the given frequencies. If the fpts array
        ia a vector (points of one function only) then F is a vector of size
        len(frequencies). If fpts is a 2D array of shape (N, rpts.size) then
        F is a 2D array of shape (N, len(frequencies)).
    '''
    np_freqs = np.asarray(frequency)

    if fpts.ndim > 1:
        out = np.empty((fpts.shape[0], np_freqs.size,), dtype=np.float)
        rpts = np.reshape(rpts, (1, rpts.size))
    else:
        out = np.empty((np_freqs.size,), dtype=np.float)

    if uneven is None:
        uneven = _is_uneven(rpts)

    r = dr = None
    if uneven:
        r = rpts
        if out.ndim > 1:
            r = np.reshape(r, (1, r.size))
    else:
        dr = rpts.flat[1] - rpts.flat[0]

    for index in range(np_freqs.size):
        f = fpts*j0(2*np.pi*np_freqs[index]*rpts)*rpts
        if out.ndim > 1:
            out[:, index] = 2*np.pi*simps(f, r, dx=dr)
        else:
            out[index] = 2*np.pi*simps(f, r, dx=dr)

    return out
