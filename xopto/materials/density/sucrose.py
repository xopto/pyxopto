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

from xopto import DATA_PATH
from .base import Density


class Lit(Density):
    material = 'sucrose'

    def __init__(self, t: float = 297.45):
        '''
        Density of sucrose as measured by:
            LST: calibrated density bottle
            (ISOLAB Laborgeräte GmbH, Eschau, Germany)

        Parameters
        ----------
        t: float
            Default temperature (K) of the medium.
            Defaults to 297.45 K (24.3 C).
        '''
        super().__init__(t=t, trange=None)

        filename = os.path.join(
            DATA_PATH, 'materials', 'density', 'sucrose_24p3c_density.npz')
        data = np.load(filename)
        self._poly = data['polyval_coefficients']

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

    def __call__(self, brix: float = 0.0, temperature: float or None = None) \
            -> float or np.ndarray:
        '''
        Calculate the density of the medium.

        Parameters
        ----------
        brix: float
            Degree Brix value (1 Brix equals 1 g of sucrose in 100 g
            of solution). Defaults to 0 Brix, i.e. pure water.
        temperature: float or None
            Use this temperature instead of the default temperature.
            Not used by this implementation.

        Returns
        -------
        rho: float or np.ndarray
            Density in kg/m^3.
        '''
        self.check_brix(brix)
        brix = brix/100.0 # relative Brix values used internally

        return np.polyval(self._poly, brix)

lit = Lit()

default = lit
