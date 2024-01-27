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

from .pfbase import PfBase
from .hg import Hg


class Hga(PfBase):
    def __init__(self, g: float or np.ndarray):
        '''
        Anisotropic Henyey-Greenstein scattering phase function constructor.

        Parameters
        ----------
        g: float or np.ndarray
            Anisotropy tensor as a 3x3 array. If a vector of 3 elements,
            the values are used to populate the diagonal of the tensor. If
            a scalar value, all elements of the diagonal are assigned the
            value (equals to :py:class:`~xopto.util.pf.hg.Hg`).
        '''

        super().__init__()
        eps = np.finfo(np.float64).eps
        self._g = np.zeros((3, 3))
        if isinstance(g, (float, int)):
            g = max(-1.0 + eps, min(1.0 -eps, g))
            self._g[0, 0] = g
            self._g[1, 1] = g
            self._g[2, 2] = g
        else:
            g = np.asarray(g, dtype=float)
            if g.size == 3:
                g = np.clip(g, -1.0 + eps, 1.0 - eps)
                self._g[0, 0] = g[0]
                self._g[1, 1] = g[1]
                self._g[2, 2] = g[2]
            else:
                g = np.clip(g, -1.0 + eps, 1.0 - eps)
                self._g[:] = g

    def tensor(self) -> np.ndarray:
        '''
        Returns
        -------
        g: np.ndarray
            Returns the anisotropy tensor.
        '''
        return self._g

    def _project_g_tensor(self, dir: np.ndarray) -> float:
        '''
        Computes/projects the g tensor along the given direction.

        Parameters
        ----------
        dir: np.ndarray
            Direction vector of length 3.
            Does not need to be normalized to unit length.
        '''
        dir = np.asarray(dir)/np.linalg.norm(dir) 

        g = np.dot(
            np.reshape(dir, (1, 3)),
            np.dot(self._g, np.reshape(dir, (3, 1)))
        )

        return float(g)

    def __call__(self, dir: np.ndarray, costheta: float or np.ndarray) \
            -> float or np.ndarray:
        '''
        Call method of the Henyey-Greenstein scattering phase function object.

        Parameters
        ----------
        dir: np.ndarray
            Direction vector of length 3 along which to evaluate the scattering
            phase function. Does not need to be normalized to unit length.
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function
            is evaluated.

        Returns
        -------
        f: float or np.ndarray
            Scattering phase function at the specified scattering angle cosines.

        '''
        g = self._project_g_tensor(dir)
        return 0.5*(1.0 - g*g)/(1 + g*g - 2*g*costheta)**1.5

    def g(self, dir: np.ndarray, n) -> float:
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with 
        an analytical solution.

        Parameters
        ----------
        dir: np.ndarray
            Direction vector of length 3 along which to evaluate the scattering
            phase function. Does not need to be normalized to unit length.
        n: int
            Legendre moment to compute.

        Returns
        -------
        gn: float
            Returns the n-th Legendre moment. 
        '''
        g = self._project_g_tensor(dir)
        return g**n

    def gs(self, dir: np.ndarray, last: int) -> np.ndarray:
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.

        Parameters
        ----------
        dir: np.ndarray
            Direction vector of length 3 along which to evaluate the scattering
            phase function. Does not need to be normalized to unit length.
        last: int
            Computes Legendre moments from 0 to last.

        Returns
        -------
        gs: np.ndarray
            A vector of computed Legendre moments. 
        '''
        g = self._project_g_tensor(dir)
        return g**np.arange(last + 1, dtype=np.float64)

    def fastg(self, dir: np.ndarray, n: int, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with
        an analytical solution.

        Parameters
        ----------
        dir: np.ndarray
            Direction vector of length 3 along which to evaluate the scattering
            phase function. Does not need to be normalized to unit length.
        n: int
            Legendre moment to compute.

        Returns
        -------
        gn: float
            Returns the n-th Legendre moment. 
        '''
        g = self._project_g_tensor(dir)
        return g**n

    def fastgs(self, dir: np.ndarray, last: int, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.

        Parameters
        ----------
        dir: np.ndarray
            Direction vector of length 3 along which to evaluate the scattering
            phase function. Does not need to be normalized to unit length.
        last: int
            Computes Legendre moments from 0 to last.

        Returns
        -------
        gs: np.ndarray
            A vector of computed Legendre moments. 
        '''
        g = self._project_g_tensor(dir)
        return g**np.arange(last + 1, dtype=np.float64)

    def cdf(self, dir: np.ndarray, costheta: np.ndarray, meth: str = 'simpson',
            npts: int = 10000, **kwargs) -> np.ndarray:
        '''
        Cumulative probability density function calculated at the specified
        deflection angle cosines.

        Parameters
        ---------
        dir: np.ndarray
            Direction vector of length 3 along which to evaluate the scattering
            phase function. Does not need to be normalized to unit length.
        costheta: np.ndarray
            A vector of scattering angle cosines at which the cumulative
            distribution function is to be computed.
        meth: str
            Numerical integration method:
                - 'trapez'  - trapezoidal
                - 'simpson' - default
                - 'quad'    - adaptive step; accurate but slower

        npts: int
            Number of control points used by the trapezoidal and simpson
            numerical integration methods.
        kwargs: dict
            Optional keyword arguments passed to the scipy.integrate.quad
            function.

        Returns
        -------
        pfcum: np.ndarray
            Returns the cumulative probability density values at the given
            scattering angle cosines.

        Note
        ----
        Computation method 'quad' might be slower on some platforms, however
        produces accurate results. Use 'quad' whenever possible. Otherwise
        check the accuracy of results obtained by 'simpson' or 'trapez' and
        adjust parameter npts if required (increase to improve accuracy,
        decrease for faster computation and increased computational error).
        '''
        g = self._project_g_tensor(dir)
        hg = Hg(g)
        return hg.cdf(costheta, meth, npts, **kwargs)

    def mclut(self, *args, **kwargs):
        raise RuntimeError('Not implemented.')

    def mcluterr(self, *args, **kwargs):
        raise RuntimeError('Not implemented.')

    def lut(self, *args, **kwargs):
        raise RuntimeError('Not implemented.')

    def __repr__(self):
        return 'Hga({})'.format(self._g)
