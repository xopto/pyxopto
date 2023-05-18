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

from .base import Absorption, DEFAULT_WAVELENGTHS
import numpy as np


class Andreia(Absorption):
    material = 'epidermis baseline'

    def __init__(self):
        '''
        Baseline absorption coefficient of epidermis:

        - Skin Optics, Oregon Medical Laser Center News, Jan 1998. Steven L. Jacques
          https://omlc.org/news/jan98/skinoptics.html

        - Andreia V. Moço, Scientific REPOrtS | (2018) 8:8501 |
          DOI:10.1038/s41598-018-26068-2
        '''
        super().__init__()

    def __call__(self, wavelength: float or np.ndarray,
                 gamma: float or np.ndarray, t: float = None) \
                     -> np.ndarray or float:
        '''
        Computes the baseline absorption coefficient (1/m) of epidermis.

        Parameters
        ----------
        wavelength: float, np.ndarray
            Wavelength of light (m) at which to compute the absorption
            coefficient.

        gamma: float
            Baseline absorption parameter from 0 to 1.

        t: float
            Temperature not supported by this model.

        Returns
        -------
        mua: float, np.ndarray
            Absorption coefficient (1/m) at the specified wavelength (m).
        '''
        return 1e2*gamma*(0.244 + 85.3*np.exp(-(1e9*wavelength - 154)/66.2))

    def plot(self, wavelength: np.ndarray, gamma : float = 0.5,
             show: bool = True):
        '''
        Plot the absorption coefficient at the given wavelengths.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths of light. If None, use the reference values.
        gamma: float
            Baseline absorption parameter.
        show: bool
            Show the plot window if True.
        '''
        import matplotlib.pyplot as pp

        if wavelength is None:
            wavelength = DEFAULT_WAVELENGTHS
        pp.semilogy(wavelength, self(wavelength, gamma), label=self.material)


andreia = Andreia()

default = andreia
