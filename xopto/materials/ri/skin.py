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

class Ding(RefractiveIndex):
    material = 'skin'

    def __init__(self, t: float = 295.0):
        '''
        Refractive index of human skin:
            "Huafeng Ding, Phys. Med. Biol.51(2006) 1479–1489,
             325 nm - 1557 nm, 22 C"

        Parameters
        ----------
        t: float
            Default temperature (K) of the medium.
        '''
        super().__init__(t=t, trange=None, wrange=(325e-9, 1557e-9))

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

        w = w*1e9

        n = 1.3549 + 1.7899e1/w - 3.5938e6/w**3.5

        if isinstance(wavelength, float):
            n = float(n)

        return n

ding = Ding()

default = ding
