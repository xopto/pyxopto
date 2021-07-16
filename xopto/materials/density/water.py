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

from .base import Density


class McCutcheon(Density):
    material = 'water'

    def __init__(self, t: float = 293):
        '''
        Density of water as given in:
        Water density as function of temperature and salt concentration
        McCutcheon, S.C., Martin, J.L, Barnwell, T.O. Jr. 1993. Water
        Quality in Maidment, D.R. (Editor). Handbook of Hydrology,
        McGraw-Hill, New York, NY (p. 11.3).

        The density of water :math:`\\rho` follows the folowing temperature T
        (C :sup:`o`) dependence.
        :math:`\\rho = 1000(1 - (T + 288.9414)/(508929.2*(T + 68.12963))*(T - 3.9863)^2)`

        Water density :math:`\\rho' as a function of temperature and salinity
        :math:`s` kg/m :sup:`3` as a function of temperature and
        salinity :math:`s` in g/kg.

        .. math::

            \\rho &= \\rho + As + Bs^{3/2} + Cs^2

            A &= 8.24493 10^{-1} - 4.0899 10^{-3}T + 7.6438 10^{-5}T^2 - 8.2467 10^{-7}T^3 + 5.3675 10^{-9}T^4

            B &= -5.724 10^{-3} + 1.0227 10^{-4}T - 1.6546 10^{-6}T^2

            C &= 4.8314 10^{-4}

        Parameters
        ----------
        t: float
            Temperature of the medium (K).
        '''
        self._salinity_range = None

        super().__init__(t=t, trange=None)

    def is_valid_salinity(self, salinity: float or np.ndarray) -> bool:
        '''
        Check if the salinity is within the valid range.

        Parameters
        ----------
        salinity: float or np.ndarray
            Water salinity from 0.0 to 1.0.

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
            Salinity from 0.0 to 0.035.
        '''
        if not self.is_valid_salinity(salinity):
            warnings.warn(
                'Salinity is out of valid range '\
                '[{:.1f}, {:.1f}] %!'.format(
                    1e2*self._salinity_range[0], 1e2*self._salinity_range[1])
            )

    def __call__(self, temperature: float or None = None,
                 salinity: float = 0.0) -> float or np.ndarray:
        '''
        Calculate the density of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)
        salinity: float
            Seawater salinity from in kg/kg.
        temperature: float or None
            Use this temperature instead of the default temperature.

        Returns
        -------
        ri: float or np.ndarray
            Refractive index at the given wavelength(s).
        '''
        if temperature is None:
            temperature = self.t
        t = float(temperature)
        self.check_temperature(t)

        salinity = float(salinity)
        self.check_salinity(salinity)

        T = temperature - 273.0 # temperature must be in C
        T2 = T*T
        T3 = T2*T
        T4 = T2*T2
        rho = 1000*(1 - (T + 288.9414)/(508929.2*(T + 68.12963))*(T - 3.9863)**2)
        S = salinity*1e3
        A = 8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T2 -8.2467e-7*T3 + 5.3675e-9*T4
        B = -5.724e-3 + 1.0227e-4*T - 1.6546e-6*T2
        C = 4.8314e-4

        rho = rho + A*S + B*S**(3/2) + C*S**2

        if isinstance(temperature, float) and isinstance(salinity, float):
            rho = float(rho)

        return rho

mccutcheon = McCutcheon()

default = mccutcheon
