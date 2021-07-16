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

import warnings

import numpy as np

from .base import RefractiveIndex


class Nootz(RefractiveIndex):
    material = 'seawater'

    def __init__(self, t=293):
        '''
        Refractive index of seawater as given in:
            Gero Nootz V1.0 07/08/2015
            "Xiaohong Quan and Edward S. Fry, Empirical equation for the index
            of refraction of seawater", Applied Optics, Vol. 34, No. 18 (1995)"
            http://www.mathworks.com/matlabcentral/fileexchange/52030-refractive-index-of-seawater

        Parameters
        ----------
        t: float
            Temperature of the medium.
        '''
        self._salinity_range = (0.0, 0.035)

        super().__init__(t=t, trange=(273, 303), wrange=(400e-9, 700e-9))

    def is_valid_salinity(self, salinity: float or np.ndarray) -> bool:
        '''
        Check if the salinity is within the valid range.

        Parameters
        ----------
        salinity: float or np.ndarray
            Water salinity (kg/kg).

        Returns
        -------
        valid: bool
            True if the salinity is within the valid range.
        '''
        if self._salinity_range is not None:
            salinity = np.asarray(salinity)
            res = np.logical_and(salinity >= self._salinity_range[0],
                                 salinity <= self._salinity_range[1])
            return np.all(res)

        return True

    def check_salinity(self, salinity: float or np.ndarray):
        '''
        Check if the salinity is within valid range and display a warning if
        not.

        Parameters
        ----------
        salinity: float
            Salinity in (kg/kg).
        '''
        if not self.is_valid_salinity(salinity):
            warnings.warn(
                'Salinity is out of valid range '\
                '[{:.1f}, {:.1f}] %!'.format(
                    1e2*self._salinity_range[0], 1e2*self._salinity_range[1])
            )

    def __call__(self, wavelength: float or np.ndarray,
                 salinity: float = 0.0, temperature: float or None = None) \
            -> float or np.ndarray:
        '''
        Calculate the refractive index of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)
        salinity: float
            Seawater salinity from in kg/kg from 0 to 0.035 (3.5 %).
        temperature: float or None
            Use this temperature instead of the default temperature.

        Returns
        -------
        ri: float or np.ndarray
            Refractive index at the given wavelength(s).
        '''
        w = np.asarray(wavelength, dtype=np.float)
        self.check_wavelengths(w)

        if temperature is None:
            temperature = self.t
        t = float(temperature)
        self.check_temperature(t)

        salinity = float(salinity)
        self.check_salinity(salinity)

        if temperature is None:
            temperature = self.t
        t = float(temperature)
        self.check_temperature(t)

        # temperature must be in C
        t = t - 273

        # the wavelength must be in nm
        w = w*1e9

        n0 = 1.31405
        n1 = 1.779e-4
        n2 = -1.05e-6
        n3 = 1.6e-8
        n4 = -2.02e-6
        n5 = 15.868
        n6 = 0.01155
        n7 = -0.00423
        n8 = -4382
        n9 = 1.1455e6

        n = n0 + \
            (n1 + n2*t + n3*t**2)*salinity + \
            n4*t**2 + (n5 + n6*salinity + n7*t)/w + \
            n8/w**2 + n9/w**3

        if isinstance(wavelength, float):
            n = float(n)

        return n

nootz = Nootz()

default = nootz
