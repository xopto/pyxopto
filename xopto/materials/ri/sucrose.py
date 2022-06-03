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

import os.path
import warnings

import numpy as np
from .base import RefractiveIndex
from xopto import DATA_PATH
from xopto.materials.ri.util.model import Sellmeier_3


class Naglic(RefractiveIndex):
    material = 'sucrose'

    def __init__(self, t: float = 293.15):
        '''
        Refractive index of sucrose as given/used in:
            "P. Naglič, Y. Zelinskyi, B. Likar, and B. Miran. Determination
            of refractive index, size, and solid content of monodisperse
            polystyrene microsphere suspensions for the characterization
            of optical phantoms, BOE, Vol. 11, Issue 4,
            1901-18 (2020)"

        Parameters
        ----------
        t: float
            Default temperature (K) of the medium. Defaults to 293.15 K (20 C).
        '''
        super().__init__(t=t, trange=(293.15, 298.15), wrange=(400e-9, 800e-9))

        filename_20c = os.path.join(
            DATA_PATH, 'materials', 'ri', 'sucrose_20c_sellmeier_ri.npz')
        filename_25c = os.path.join(
            DATA_PATH, 'materials', 'ri', 'sucrose_25c_sellmeier_ri.npz')

        wavelength_pp = lambda wavelength: wavelength*1e6

        data_20c = np.load(filename_20c)
        self._brix_20c = data_20c['brix'] 
        self._ri_sellmeir_20c = \
            [Sellmeier_3(item, wavelength_pp) for item in data_20c['coeffs']]

        data_25c = np.load(filename_25c)
        self._brix_25c = data_25c['brix'] 
        self._ri_sellmeir_25c = \
            [Sellmeier_3(item, wavelength_pp) for item in data_25c['coeffs']]

        self._brix_range = (0.0, 60.0)

    def is_valid_brix(self, brix: float or np.ndarray) -> bool:
        '''
        Check if the Brix value is within the valid range.

        Parameters
        ----------
        brix: float or np.ndarray
            Degree Brix value (1 Brix equals 1 g of sucrose in 100 g
            of solution). Defaults to 0 Brix, i.e. pure water.

        Returns
        -------
        valid: bool
            True if the Brix value is within the valid range.
        '''
        if self._brix_range is not None:
            brix = np.asarray(brix)
            res = np.logical_and(brix >= self._brix_range[0],
                                 brix <= self._brix_range[1])
            return np.all(res)

        return True

    def check_brix(self, brix: float or np.ndarray):
        '''
        Check if the Brix value is within valid range and display a warning if
        not.

        Parameters
        ----------
        brix: float or np.ndarray
            Degree Brix value (1 Brix equals 1 g of sucrose in 100 g
            of solution). Defaults to 0 Brix, i.e. pure water.
        '''
        if not self.is_valid_brix(brix):
            warnings.warn(
                'Brix value out of valid range '\
                '[{:.1f}, {:.1f}] Brix!'.format(
                    self._brix_range[0], self._brix_range[1])
            )

    def __call__(self, wavelength: float or np.ndarray, brix: float = 0.0,
                 temperature: float or None = None) -> float or np.ndarray:
        '''
        Calculate the refractive index of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m). Valid range is approx. from 400 to 800 nm.
        brix: float
            Degree Brix value (1 Brix equals 1 g of sucrose in 100 g
            of solution). Defaults to 0 Brix, i.e. pure water.
        temperature: float or None
            Temperature (K) of the medium. Valid range is
            from 293 K (20 C) to 298 K (25 C). Defaults to 293.15 K.

        Returns
        -------
        ri: float or np.ndarray
            Refractive index at the given wavelength(s) and Brix value.
        '''
        if temperature is None:
            temperature = self.t
        self.check_temperature(temperature)

        self.check_brix(brix)
        brix = brix/100.0 # using relative Brix values internally

        w = np.asarray(wavelength, dtype=np.float64)
        self.check_wavelengths(w)

        ri_20c = [ri(w) for ri in self._ri_sellmeir_20c]
        ri_25c = [ri(w) for ri in self._ri_sellmeir_25c]

        # model uses C
        t = temperature - 273.15

        p_20c = np.polyfit(self._brix_20c, ri_20c, deg=4)
        p_25c = np.polyfit(self._brix_25c, ri_25c, deg=4)

        ri_20c = np.polyval(p_20c, brix)
        ri_25c = np.polyval(p_25c, brix)

        result = (ri_25c - ri_20c)/(25.0 - 20.0) * (t - 20.0) + ri_20c

        if isinstance(wavelength, float):
            result = float(result)

        return result

# Global instances for immediate use.
naglic = Naglic()

default = Naglic()
