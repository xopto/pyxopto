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

from ..base import RefractiveIndex

class Schott(RefractiveIndex):
    material = 'nbk7'
    def __init__(self, t=293):
        '''
        Refractive index of NBK-7 glass:
            https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT

        Parameters
        ----------
        t: float
            Temperature of the medium.
        '''
        super().__init__(t=t, trange=None, wrange=(300e-9, 2500e-9))

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

        w = w*1e6
        w2 = w*w
        n2 = 1.0 + 1.03961212*w2/(w2 - 0.00600069867) + \
             0.231792344*w2/(w2 - 0.0200179144) + \
             1.01046945*w2/(w2 - 103.560653)
        
        n = np.sqrt(n2)

        if isinstance(wavelength, float):
            n = float(n)

        return n

schott = Schott()

default = schott
