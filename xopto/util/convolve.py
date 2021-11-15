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
from scipy.integrate import simps

def fiber_reflectance(r: np.ndarray, reflectance: np.ndarray,
                      sds: float or np.ndarray, dcore: float,
                      nsimps: int = 1000, loginterp: bool = False,
                      uneven: bool = False, **kwargs) -> np.ndarray:
    '''
    Computes reflectance detected through optical fibers with the given
    geometry. The reflectance can be computed by integrating:

    .. math::

        & 2 \\int_{\\max(sds - r_c, 0)}^{\\max(sds + r_c, 0)} \\arccos\\left(\\frac{r^2 + sds^2 - r_c^2}{2 r_c sds}\\right) \\text{reflectance}(r) r dr,

        & \\text{where}\\; r_c = dcore/2 

    Parameters
    ----------
    r: np.ndarray vector
        Vector of points at which the radially symmetric reflectance is defined.
        The points must be sorted in ascending order and include the range
        :math:`[sds - dcore/2, sds + dcore/2]`.
    reflectance: np.ndarray vector or 2D array
        Radially symmetric reflectance at points r.
        To compute the reflectance through fibers for multiple reflectances,
        use a 2D array with individual reflectances in the rows (n, r.size).
    sds: float or list/tuple of float
        Distance to the fiber center.
    dcore: float
        Fiber diameter.
    nsimps: int
        Number of integration points used by the simpson method.
    loginterp: bool
        If True, the interpolation of reflectance is
        computed in log space. This can in some cases improve the accuracy of
        integration.
    uneven: bool
        Set to true if unevenly spaced points can be found in r.
    kwargs: dict
        Keyword arguments are passed to scipy.interpolate.interp1d.

    Returns
    -------
    fiber_reflectance: np.ndarray vector or 2D array
        Reflectance collected through the fibers. The size of the vector
        equals the number of elements in sds. If reflectance is a 2D array
        then the results are stored in rows (n, len(sds)).
    '''
    if isinstance(sds, float):
        sds = (sds,)

    rfib = dcore*0.5
    nsimps = max(nsimps//2*2 + 1, 3)

    reflectance = np.asarray(reflectance)
    if reflectance.ndim > 1:
        num_reflectances = reflectance.shape[0]
    else:
        num_reflectances = 1

    # TODO should we do interpolation in log scale?
    if loginterp:
        f = interp1d(r, np.log(reflectance),
                     assume_sorted=True, bounds_error=False,
                     fill_value='extrapolate', **kwargs)
    else:
        f = interp1d(r, reflectance,
                     assume_sorted=True, bounds_error=False,
                     fill_value='extrapolate', **kwargs)

    fiber_reflectance = np.zeros((num_reflectances, len(sds)))

    for index, d in enumerate(sds):

        rsimps = np.linspace(max(d - rfib, 0.0), max(d + rfib, 0.0), nsimps)
        if uneven:
            simps_dr = None
            simps_r = rsimps
        else:
            simps_dr = rsimps[1] - rsimps[0]
            simps_r = None

        if loginterp:
            fsimps = np.exp(f(rsimps))
        else:
            fsimps = f(rsimps)


        # Take care of special case when sds = 0.0.
        if d == 0.0:
            fiber_reflectance[:, index] = 2.0*np.pi*simps(
                fsimps*rsimps, x=simps_r, dx=simps_dr
            )

        else:
            # First element of rsimps might be zero ... precision!
            d_times_rsimps = rsimps*d
            d_times_rsimps[d_times_rsimps == 0.0] = np.finfo(np.float).tiny
            fiber_reflectance[:, index] = 2.0*simps(
                np.arccos(
                    np.clip(
                        (rsimps**2 + d**2 - rfib**2)/(2.0*d_times_rsimps),
                        -1.0, 1.0
                    )
                )*fsimps*rsimps, x=simps_r, dx=simps_dr
            )

    return fiber_reflectance
