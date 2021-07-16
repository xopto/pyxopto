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


class Malitson(RefractiveIndex):
    material = 'pmma'
    def __init__(self, t=293):
        '''
        Refractive index of fused silica SiO2:
            "I. H. Malitson. Interspecimen comparison of the refractive index
            of fused silica, J. Opt. Soc. Am. 55, 1205-1208 (1965)"

        Parameters
        ----------
        t: float
            Temperature of the medium.
        '''
        super().__init__(t=t, trange=None, wrange=(210e-9, 3710e-9))

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

        w = w*1e6
        w2 = w*w
        n2 = 1.0 + w2*(0.6961663/(w2 - 0.0684043**2) + \
             0.4079426/(w2 - 0.1162414**2) + 0.8974794/(w2 - 9.896161**2))
    
        n = np.sqrt(n2)

        if isinstance(wavelength, float):
            n = float(n)

        return n

malitson = Malitson()

default = malitson
