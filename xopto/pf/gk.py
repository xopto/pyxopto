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


class Gk(PfBase):
    def __init__(self, gg: float, a: float):
        '''
        Gegenbauer Kernel scattering phase function constructor.

        Parameters
        ----------
        gg: float
            Parameter of the Gegenbauer kernel phase function
            (:math:`|gg| <= 1`).
        a: float
            Parameter of the Gegenbauer kernel phase function
            (:math:`a > - 1/2`).
            A value of 0.5 produces the Henyey-Greenstein scattering
            phase function.

        Examples
        --------
        Gegenbauer kernel scattering phase function for
        gg = {0, 0.3 0.5, 0.8, 0.9, 0.95} and a=0.5.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for gg in [0.0, 0.3, 0.5, 0.8, 0.9, 0.95]:
        >>>     pp.semilogy(cos_theta, Gk(gg, 0.5)(cos_theta), label='a=0.5, gg={}'.format(gg))
        >>> pp.legend()
        '''
        super().__init__()
        self._gg = gg = min(gg, 1.0 - np.finfo(float).eps)
        self._a = a = float(a)
        if gg == 0.0:
            def pf(costheta):
                return np.tile(0.5, np.asarray(costheta).shape)
            self._pf = pf
        elif a == 0:
            def pf(costheta):
                K = gg/(np.log((1.0 + gg)/(1.0 - gg)))
                return K/(1.0 + gg*gg - 2.0*gg*costheta)
            self._pf = pf
        else:
            def pf(costheta):
                K = 2*a*gg*(1.0 - gg*gg)**(2.0*a)/((1.0 + gg)**(2.0*a) - \
                    (1.0 - gg)**(2.0*a))
                return K*(1.0 + gg*gg - 2.0*gg*costheta)**(-1.0 - a)
            self._pf = pf

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Gegenbauer kernel scattering phase function.

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
        return self._pf(costheta)

    def __repr__(self):
        return 'Gk({}, {})'.format(self._gg, self._a)

