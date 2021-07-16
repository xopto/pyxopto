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

class Naglic(RefractiveIndex):
    material = 'siliglass'

    def __init__(self, t=293):
        '''
        Refractive index of PlatSil® SiliGlass:
            "P. Naglič, Y. Zelinskyi, L. Rogelj, J. Stergar, M. Milanič,
            B. Kumperščak and M. Bũrmen. Optical properties of PlatSil
            SiliGlass tissue-mimicking phantoms, BOE, Vol. 11, Issue 7,
            3753-68 (2020)"

        Parameters
        ----------
        t: float
            Temperature of the medium.
        '''
        super().__init__(t=t, trange=(293.0, 303.0), wrange=(400e-9, 1000e-9))

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
        w = np.asarray(wavelength, dtype=np.float)
        self.check_wavelengths(w)

        if temperature is None:
            temperature = self.t
        t = np.asarray(temperature, dtype=np.float)
        self.check_temperature(t)


        inv_wn2 = 1.0/(np.asarray(w)/4.364e-07)**2
        t_c = t - 273.0

        n_20 = np.sqrt(1.0 + 0.31669505/(1.0 - 0.05941627*inv_wn2) +
                            0.31669505/(1.0 - 0.05941630*inv_wn2) +
                            0.31669506/(1.0 - 0.05941629*inv_wn2))
        n_25 = np.sqrt(1.0 + 0.31504349/(1.0 - 0.05943151*inv_wn2) +
                            0.31504349/(1.0 - 0.05943150*inv_wn2) +
                            0.31504349/(1.0 - 0.05943151*inv_wn2))
        n_30 = np.sqrt(1.0 + 0.31330351/(1.0 - 0.05907095*inv_wn2) +
                            0.31330350/(1.0 - 0.05907094*inv_wn2) +
                            0.31330351/(1.0 - 0.05907094*inv_wn2))
        
        x = (t_c - 20.0)/10.0
        
        n = (2.0*n_30 + 2.0*n_20 - 4.0*n_25)*x**2 + \
            (4.0*n_25 - 3.0*n_20 - n_30)*x + (n_20)

        if isinstance(wavelength, float) and isinstance(temperature, float):
            n = float(n)

        return n

naglic = Naglic()

default = naglic
