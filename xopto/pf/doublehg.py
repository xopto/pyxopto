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

from .pfbase import PfBase


class DoubleHg(PfBase):
    def __init__(self, g1: float, g2: float, b: float):
        '''
        Double Henyey-Greenstein scattering phase function constructor.

        Parameters
        ----------
        g1: float
            Anisotropy factor of the first Henyey-Greenstein scattering phase function.
        g2: float
            Anisotropy factor of the second Henyey-Greenstein scattering phase function.
        b: float
            Fractional contribution of the first Henyey-Greenstein phase function.

        Examples
        --------
        Double Henyey-Greenstein scattering phase function for anisotropy
        factors g = {0, 0.3 0.5, 0.8, 0.9, 0.95}.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for g in [0, 0.3, 0.5, 0.8, 0.9, 0.95]:
        >>>     pf = DoubleHg(g, -g, 0.5)
        >>>     pp.semilogy(cos_theta, pf(cos_theta), label='g1={}, g2={}, b=0.5'.format(g, -g))
        >>> pp.legend()
        '''

        super().__init__()
        eps = np.finfo(np.float64).eps
        self._g1 = max(min(g1, 1.0 - eps), -1 + eps)
        self._g2 = max(min(g2, 1.0 - eps), -1 + eps)
        self._b = min(max(float(b), 0.0), 1.0)

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Double Henyey-Greenstein scattering phase function.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function
            is evaluated.

        Returns
        -------
        f: float or np.ndarray
            Scattering phase function at the specified scattering angle cosines.
        '''
        return self._b*0.5*(1.0 - self._g1*self._g1)/ \
                    (1 + self._g1*self._g1 - 2*self._g1*costheta)**1.5 + \
               (1.0 - self._b)*0.5*(1.0 - self._g2*self._g2)/ \
                    (1 + self._g2*self._g2 - 2*self._g2*costheta)**1.5

    def __repr__(self):
        return 'DoubleHg(g1={}, g2={}, b={})'.format(
            self._g1, self._g2, self._b)
