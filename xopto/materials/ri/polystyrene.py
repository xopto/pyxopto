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

from .base import RefractiveIndex

class Ma(RefractiveIndex):
    material = 'polystyrene'

    def __init__(self, t: float = 293.15):
        '''
        Refractive index of water as given in:
            Ma, Lu, Brock: Phys. Med. Biol. 48 (2003) 4165-4172.

        .. math::

            n &= 1.5725 + 0.0031080/(\\lambda*10^{6})^2 + 0.00034779/(\\lambda*10^{6})^4

        Parameters
        ----------
        t: float
            Defau.t temperature (K) of the medium.
        '''
        super().__init__(t=t, trange=None, wrange=(370e-9, 1610e-9))

    def __call__(self, wavelength: float or np.ndarray,
                 temperature: float or None = None) -> float or np.ndarray:
        '''
        Calculate the refractive index of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)
        temperature: float or None
            Use this temperature instead of the default temperature.

        Returns
        -------
        ri: float or np.ndarray
            Refractive index at the given wavelength(s).
        '''
        w = np.asarray(wavelength, dtype=np.float64)
        self.check_wavelengths(w)

        w2 = (w*1e6)**2 # wavelength must be in units of um
        n = 1.5725 + 0.0031080/w2 + 0.00034779/(w2**2)

        if isinstance(wavelength, float):
            n = float(n)

        return n


class Nikolov(RefractiveIndex):
    material = 'polystyrene'

    def __init__(self, t: float = 293.15):
        '''
        Refractive index of water as given in:
            Nikolov, Ivanov: Appl. Opt. 39 (2000) 2067-2070

        Refractive indices
        :math:`n = [1.617, 1.606, 1.596, 1.592, 1.587, 1.586, 1.582, 1.579, 1.578, 1.577, 1.576, 1.572]`
        measured at :math:`\\lambda = [436, 486, 546, 588, 633, 656, 703, 752, 804, 833, 879, 1052]` nm.

        Parameters
        ----------
        t: float
            Default temperature (K) of the medium.

        Note
        ----
        The computation of refractive index is based on the one-term Sellmeier
        approximation:
        .. math::

            n = \\sqrt(1.0 + a/(1.0 - b/\\lambda^2))

        '''
        super().__init__(t=t, trange=None, wrange=(442e-9, 1060e-9))

    def __call__(self, wavelength: float or np.ndarray,
                 temperature: float or None = None) -> float or np.ndarray:
        '''
        Calculate the refractive index of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)
        temperature: float or None
            Use this temperature instead of the default temperature.

        Returns
        -------
        ri: float or np.ndarray
            Refractive index at the given wavelength(s).
        '''
        w = np.asarray(wavelength, dtype=np.float64)
        self.check_wavelengths(w)

        n = np.sqrt(
            1.0 + 1.4431234515349096/(1.0 - 2.0219196759502886e-14/(w**2))
        )

        if isinstance(wavelength, float):
            n = float(n)

        return n

ma = Ma()

nikolov = Nikolov()

default = nikolov
