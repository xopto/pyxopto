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
from .gk import Gk

class MGk(PfBase):
    def __init__(self, gg: float, a: float, b: float):
        '''
        Modified Gegenbauer kernel scattering phase function.

        Parameters
        ----------
        gg: float
            Parameter of the Gegenbauer kernel scattering phase function
            (:math:`|gg| <= 1`).
        a: float
            Parameter of the Gegenbauer kernel scattering phase function
            (:math:`a > - 1/2`).
            A value of 0.5 produces the Henyey-Greenstein scattering
            phase function.
        b: float
            Fractional contribution of the Gegenbauer kernel
            scattering phase function.

        Examples
        --------
        Modified Gegenbauer kernel scattering phase function fo anisotropy
        factors g = {0, 0.3 0.5, 0.8, 0.9, 0.95}, a = 0.5 and b=0.5.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>> a = 0.5
        >>>
        >>> pp.figure()
        >>> for g in [0.0, 0.3, 0.5, 0.8, 0.9, 0.95]:
        >>>     pp.semilogy(cos_theta, MGk(g, 0.5)(cos_theta), label='b=0.5, g={}'.format(g))
        >>> pp.legend()
        '''
        super().__init__()
        self._b = np.clip(float(b), 0.0, 1.0)
        self._gk = Gk(gg, a)

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Modified Gegenbauer kernel scattering phase function.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function
            is evaluated.

        Returns
        -------
        f: float or np.ndarray
            Phase function at the specified deflection angle cosines.
        '''
        return self._b*self._gk(costheta) + \
            (1.0 - self._b)*3.0/(2.0)*costheta**2

    def __repr__(self):
        return 'MGk({}, {}, {})'.format(self._gk._gg, self._gk._a, self._b)
