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

from typing import Tuple
import warnings

import numpy as np

class RefractiveIndex:
    def __init__(self, t=293, trange: Tuple[float, float]=None,
                 wrange: Tuple[float, float]=None):
        '''
        Base class of refractive index.

        Parameters
        ----------
        t: float
            Default medium temperature.
        trange: Tuple[float, float] or None
            Valid temperature range in K as (min, max).
        wrange: Tuple[float, float] or None
            Valid wavelength range in nm as (min, max).
        '''
        self._temperature = float(t)

        if isinstance(trange, (float, int)):
            trange = (float(trange), float(trange))
        elif trange is not None:
            trange = (float(trange[0]), float(trange[1]))
        self._temperature_range = trange

        if isinstance(wrange, (float, int)):
            wrange = (float(wrange), float(wrange))
        elif wrange is not None:
            wrange = (float(wrange[0]), float(wrange[1]))
        self._wavelength_range = wrange

    def is_valid_temperature(self, temperature: float or np.ndarray) -> bool:
        '''
        Checks if the temperature of the medium is within the valid range.

        Parameters
        ----------
        temperature: float or np.ndarray
            Temperature to check.

        Returns
        -------
        valid: bool
            Returns True if the temperature is within the valid range.
            If the temperature input argument is a numpy array, the call
            returns True only if all the temperatures are within the valid
            range.
        '''
        if self._temperature_range is not None:
            t = np.asarray(temperature)
            res = np.logical_and(t >= self._temperature_range[0],
                                 t <= self._temperature_range[1])
            return np.all(res)

        return True

    def is_valid_wavelength(self, wavelength: float or np.ndarray) -> bool:
        '''
        hecks if the wavelength of light is within the valid range.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength(s) to check.

        Returns
        -------
        valid: bool
            Returns True if the wavelength is within the valid range.
            If the wavelength input argument is a numpy array, the call
            returns True only if all the wavelengths are within the valid
            range.
        '''
        if self._wavelength_range is not None:
            w = np.asarray(wavelength)
            res = np.logical_and(w >= self._wavelength_range[0],
                                 w <= self._wavelength_range[1])
            return np.all(res)

        return True

    def plot(self, temperature=None, **kwargs):
        '''
        Plot the wavelength dependence of the refractive index at the given
        temperature.

        Parameters
        ----------
        temperature: float or None
            Temperature of the medium or None. Default medium temperature is
            used if None.
        kwargs: dict
            Keyword arguments passed to the __call__ method of the instance. 
        '''
        import matplotlib.pyplot as pp

        if temperature is None:
            temperature = self.t
        t = float(temperature)

        wavelengths = np.linspace(
            self._wavelength_range[0], self._wavelength_range[1], 100
        )
        pp.plot(wavelengths, self(wavelengths*1e6, temperature=t, **kwargs))
        pp.xlabel('Wavelength (um)')
        pp.title('{:s}.{:s}'.format(self.material, self.__class__.__name__))
        pp.show()

    def check_wavelengths(self, w: float or np.ndarray):
        '''
        Checks if the wavelength of light is within the valid range.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength(s) to check.

        Note
        ----
        Reports a warning if the wavelength of light is not within the
        valid range. If the wavelength input argument is a numpy array,
        the call reports a warning if any of the wavelengths are not
        within the valid range.
        '''
        if not self.is_valid_wavelength(w):
            warnings.warn(
                'One or more wavelengths are out of valid range '\
                '[{:.1f}, {:.1f}] nm!'.format(
                    self._wavelength_range[0]*1e9,
                    self._wavelength_range[1]*1e9)
            )

    def check_temperature(self, t: float or np.ndarray):
        '''
        Checks if the temperature of medium is within the valid range.

        Parameters
        ----------
        temperature: float or np.ndarray
            Temperature(s) to check.

        Note
        ----
        Reports a warning if the temperature of the medium is not within the
        valid range. If the temperature input argument is a numpy array,
        the call reports a warning if any of the temperatures are not
        within the valid range.
        '''
        if not self.is_valid_temperature(t):
            warnings.warn(
                'Temperature is out of valid range '\
                '[{:.1f}, {:.1f}] K!'.format(*self._temperature_range)
            )

    def _get_temperature(self) -> float:
        return self._temperature

    def _set_temperature(self, t: float):
        t = float(t)
        if not self.is_valid_temperature(t):
            warnings.warn('Default temperature {:.1f} K is out of valid range '\
                          '[{:.1f}, {:.1f}] K!'.format(*self._temperature_range))
        self._temperature = t

    temperature = property(_get_temperature, _set_temperature, None,
                           'Default temperature of the material in K.')

    t = property(_get_temperature, _set_temperature, None,
                 'Default temperature of the material in K.')

    def _get_temperature_range(self) -> Tuple[float, float] or None:
        return self._temperature_range

    trange = property(_get_temperature_range, None, None,
                      'Valid temperature range in K.')

    def _get_wavelength_range(self) -> Tuple[float, float] or None:
        return self._wavelength_range

    wrange = property(_get_wavelength_range, None, None,
                      'Valid wavelength range in K.')

    def __str__(self):
        return '{:s}.{:s}(t={})'.format(
            self.material, self.__class__.__name__, self._temperature)

    def __repr__(self):
        return '{} # id {}'.format(self.__str__(), id(self))
