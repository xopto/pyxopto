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
from .gk     import Gk

class Gk2(PfBase):
    def __init__(self, gg1: float, a1: float, gg2: float, a2: float, b: float):
        '''
        Two-term Gegenbauer Kernel scattering phase function constructor.

        Parameters
        ----------
        gg1: float
            Parameter of the first Gegenbauer kernel phase function
            (:math:`0 <= gg_1 <= 1`).
        a1: float
            Parameter of the first Gegenbauer kernel phase function
            (:math:`a > - 1/2`).
            A value of 0.5 produces the Henyey-Greenstein scattering
            phase function.
        gg2: float
            Parameter of the second Gegenbauer kernel phase function
            (:math:`-1.0 <= gg_2 <= 0`).
        a2: float
            Parameter of the second Gegenbauer kernel phase function
            (:math:`a > - 1/2`).
            A value of 0.5 produces the Henyey-Greenstein scattering
            phase function.
        b: float
            Contribution of the second Gk.

        Examples
        --------
        Two-term Gegenbauer kernel scattering phase function for
        gg = {0, 0.3 0.5, 0.8, 0.9, 0.95} and a=0.5.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for gg in [0.0, 0.3, 0.5, 0.8, 0.9, 0.95]:
        >>>     pp.semilogy(cos_theta, Gk2(gg, 0.5, -gg, 0.5, 0.1)(cos_theta), label='a=0.5, gg={}'.format(gg))
        >>> pp.legend()
        '''
        super().__init__()
        if gg1 < 0:
            raise ValueError('Parameter gg of the first Gk must be positive!')
        if gg2 > 0:
            raise ValueError('Parameter gg of the second Gk must be negative!')

        self._gk1 = Gk(gg1, a1)
        self._gk2 = Gk(gg2, a2)
        self._b = float(b)
        
        def pf(costheta):
            return (1.0 - self._b)*self._gk1(costheta) + self._b*self._gk2(costheta)

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

    def g(self, n, *args, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with 
        an analytical solution.
        '''
        return (1.0 - self._b)*self._gk1.g(n, *args, **kwargs) + \
               self._b*self._gk2.g(n, *args, **kwargs)

    def fastg(self, n, *args, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.fastg` method of the base class with 
        an analytical solution.
        '''
        return (1.0 - self._b)*self._gk1.fastg(n, *args, **kwargs) + \
               self._b*self._gk2.fastg(n, *args, **kwargs)

    def __repr__(self):
        return 'Gk2({}, {}, {}, {}, {})'.format(
            self._gk1._gg, self._gk1._a, self._gk2._gg, self._gk2._a, self._b)
