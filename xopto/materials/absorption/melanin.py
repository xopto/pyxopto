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

class Jacques(Absorption):
    material = 'melanin'

    def __init__(self):
        '''
        Absorption coefficient (1/m) of melanin as a function of
        wavelength of light (m).

            - Jacques SL1, McAuliffe DJ, hotochem Photobiol. 1991 Jun;53(6):769-75.
            - Skin Optics, Oregon Medical Laser Center News, Jan 1998. Steven L. Jacques
                https://omlc.org/news/jan98/skinoptics.html
        '''
        super().__init__()

    def __call__(self, wavelength: float or np.ndarray, t: float = 293.15) \
            -> np.ndarray or float:
        '''
        Absorption coefficient (1/m) of melanin as a function of
        wavelength of light (m).

        Parameters
        ----------
        wavelength: float, np.ndarray
            Wavelength of light (m) at which to compute the melanin absorption
            coefficient.

        t: float
            Temperature not supported by this model.

        Returns
        -------
        mua: float, np.ndarray
            Absorption coefficient (1/m) at the specified wavelength (m).
        '''
        return 6.6*1e13*(1e9*wavelength)**(-3.33)

    def plot(self, wavelength: np.ndarray, gamma : float = 0.5,
             show: bool = True):
        '''
        Plot the absorption coefficient at the given wavelengths.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths of light. If None, use the reference values.
        show: bool
            Show the plot window if True.
        '''
        import matplotlib.pyplot as pp

        if wavelength is None:
            wavelength = DEFAULT_WAVELENGTHS
        pp.semilogy(wavelength, self(wavelength), label=self.material)
        if show:
            pp.show()

jacques = Jacques()
default = jacques
