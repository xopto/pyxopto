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


class Lit(Density):
    material = 'polystyrene'

    def __init__(self, t=293):
        '''
        Density of water as given in:
            LST

        Parameters
        ----------
        t: float
            Temperature of the medium.
        '''

        super().__init__(t=t, trange=None)

    def __call__(self, temperature: float or None = None) \
            -> float or np.ndarray:
        '''
        Calculate the density of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)
        temperature: float or None
            Use this temperature instead of the default temperature.

        Returns
        -------
        rho: float or np.ndarray
            Density in kg/m^3.
        '''
        return 1.05e3
        
lit = Lit()

default = lit
