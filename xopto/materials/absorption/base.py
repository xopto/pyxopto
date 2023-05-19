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
from scipy.interpolate import interp1d


DEFAULT_WAVELENGTHS = np.linspace(400e-9, 1000e-9, 601)
''' Default wavelengths for plots. '''

class Interpolator:
    @classmethod
    def fromfile(cls, filename: str, *args, **kwargs) -> 'Interpolator':
        '''
        Create an interpolator from a numpy data file. The independent variable
        (wavelengths) should be stored in the first column, the dependent
        variable (absorption coefficient) in the second column

        Parameters
        ----------
        filename: str
            Source file.
        args: tuple
            Positional arguments passed to 
            `:py:func:scipy.interpolate.interp1d`.
        kwargs: dict
            Keyword arguments passed to 
            `:py:func:scipy.interpolate.interp1d`.

        Returns
        -------
        interpolator: Interpolator
            Interpolator instance.
        '''
        data = np.load(filename)
        return cls(data[:, 0], data[:, 1], *args, **kwargs)

    def __init__(self, wavelength: np.ndarray, absorption: np.ndarray,
                 *args, **kwargs):
        '''
        Prepare an interpolator for the given x and y data.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths (m) of light at which the absorption coefficient is
            defined.
        absorption: np.ndarray
            Absorption coefficient (1/m) at the wavelengths of light.
            Points that will define the interpolation function.
        args, kwargs: tuple, dict
            Optional positional and keyword arguments passed to the interp1d
            function of scipy.interpolate module.
        '''
        self._wavelength = np.array(wavelength)
        self._wavelength_range = (float(wavelength.min()), float(wavelength.max()))
        self._absorption = np.array(absorption)
        self._interpolator = interp1d(
            self._wavelength, self._absorption, *args, **kwargs)

    def __call__(self, wavelength: np.ndarray or float) -> np.ndarray or float:
        '''
        Parameters
        ----------
        wavelength: float, np.ndarray
            The wavelengths of light (m) at which to estimate the absorption
            coefficient.

        Returns
        -------
        result: float, np.ndarray
            The estimated values of the absorption coefficient.
        '''
        res = self._interpolator(wavelength)

        if isinstance(wavelength, (float, int)):
            res = float(res)

        return res

    def is_valid_range(self, wavelength: float or np.ndarray) -> bool:
        '''
        Check the wavelengths of light (m) for valid range.

        Parameters
        ----------
        wavelength: float, np.ndarray
            Wavelengths of light to check for valid range.

        Returns
        -------
        ok: bool
            Returns True if all the wavelengths are within the valid range.
        '''
        return not (np.all(wavelength >= self._wavelength_range[0]) and \
                    np.all(wavelength <= self._wavelength_range[1]))

    def check_range(self, wavelength: float or np.ndarray):
        '''
        Checks if the wavelength of light is within the valid range.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Temperature(s) to check.

        Note
        ----
        Reports a warning if the wavelength of light is not within the
        valid range. If the wavelength input argument is a numpy array,
        the call reports a warning if any of the wavelengths are not
        within the valid range.
        '''
        if not self.is_valid_range(wavelength):
            warnings.warn(
                'Wavelength is out of valid range '\
                '[{:.1f}, {:.1f}] nm!'.format(
                    self._wavelength_range[0]*1e9,
                    self._wavelength_range[1]*1e9)
            )

    def _get_wavelength(self):
        return self._wavelength
    wavelength = property(_get_wavelength, None, None,
                          'Wavelengths of light (m) at which the absorption '
                          'coefficient was measured.')

    def _get_absorption(self):
        return self._absorption
    absorption = property(_get_wavelength, None, None,
                          'Measured values of the absorption coefficient (1/m).')

    def range(self):
        '''
        Return the valid range of wavelengths.

        Returns
        -------
        range: (float, float)
            Valid range of wavelengths as a tuple (low, high).
        '''
        return self._range

    def plot(self, wavelength: np.ndarray = None, label=None, show=False):
        '''
        Plots the wavelength dependence of the absorption coefficient.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths of light.
        label: str
            Plot label.
        show: bool
            Show the plot window if True.
        '''
        import matplotlib.pyplot as pp

        if wavelength is None:
            wavelength = self.wavelength
        pp.semilogy(wavelength*1e9, self(wavelength), label=label)
        pp.xlabel('Wavelength (nm)')
        pp.ylabel('Absorption coefficient (1/m)')

        if show:
            pp.show()


class Absorption:
    def __init__(self, t=293, trange: Tuple[float, float]=None):
        '''
        Base class of material density.

        Parameters
        ----------
        t: float
            Default medium temperature.
        trange: Tuple[float, float] or None
            Valid temperature range in K as (min, max).
        '''
        self._temperature = float(t)

        if isinstance(trange, (float, int)):
            trange = (float(trange), float(trange))
        elif trange is not None:
            trange = (float(trange[0]), float(trange[1]))
        self._temperature_range = trange

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

    def plot(self, **kwargs):
        '''
        Plot the temperature dependence of the density.

        Parameters
        ----------
        kwargs: dict
            Keyword arguments passed to the :py:meth:`Absorption.__call__`
            method of the instance. 
        '''
        import matplotlib.pyplot as pp

        temperature = np.linspace(
            self._temperature_range[0], self._temperature_range[1], 100
        )
        pp.plot(temperature, self(temperature, **kwargs))
        pp.xlabel('Temperature (K)')
        pp.ylabel('Density (kg/m^3)')
        pp.title('{:s}.{:s}'.format(self.material, self.__class__.__name__))
        pp.show()


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

    def __str__(self):
        return '{:s}.{:s}(t={})'.format(
            self.material, self.__class__.__name__, self._temperature)

    def __repr__(self):
        return '{} # id {}'.format(self.__str__(), id(self))
