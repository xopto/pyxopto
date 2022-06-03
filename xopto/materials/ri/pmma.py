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

class Beadie(RefractiveIndex):
    material = 'pmma'
    def __init__(self, t=293.1):
        '''
        Refractive index of (C5O2H8)n (Poly(methyl methacrylate), PMMA)

        Refractive index of water as given in:
            "G. Beadie, M. Brindza, R. A. Flynn, A. Rosenberg, and
            J. S. Shirk. Refractive index measurements of
            poly(methyl methacrylate) (PMMA) from 0.4-1.6 μm,
            Appl. Opt. 54, F139-F143 (2015)"

        Parameters
        ----------
        t: float
            Default temperature (K) of the medium.
        '''
        super().__init__(t=t, trange=None, wrange=(420e-9, 1620e-9))

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

        if temperature is None:
            temperature = self.t
        t = float(temperature)
        self.check_temperature(t)

        w = w*1e6

        w2 = w*w
        w4 = w2*w2
        w6 = w4*w2
        w8 = w4*w4
        n2 = 2.1778 + 6.1209e-3*w2 - 1.5004e-3*w4 + 2.3678e-2/w2 - \
            4.2137e-3/w4 + 7.3417e-4/w6 - 4.5042e-5/w8

        if t != 293.1:
            wT = np.array([0.472, 0.780, 1.0557, 1.3089])*1e-9
            kT = np.array([-1.37, -1.37, -1.30, -1.33])*1e-4
            kt = np.interp(w, wT, kT)
            n = np.sqrt(n2)*(1.0 + kt*(t - 293.1))
        else:
            n = np.sqrt(n2)

        if isinstance(wavelength, float):
            n = float(n)

        return n


beadie = Beadie()

default = beadie 
