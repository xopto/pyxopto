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


class MPc(PfBase):
    def __init__(self, n: float, b: float):
        '''
        Modified power of cosines scattering phase function constructor.

        .. math::

            Pc(\\cos(\\theta)) &= \\frac{(n + 1)}{2}^{(n + 1)} * (1 + \\cos(\\theta))^n

            MPc(\\cos(\\theta)) &= b Pc (\\cos(\\theta)) + \\frac{3}{2}(1 - b)\\cos(\\theta)^2

        Parameters
        ----------
        n: float
            Parameter of the power of cosine scattering phase function.
        b: float
            Contribution of the Power of cosine scattering component.
            Contribution of the Rayleigh scattering component is (1 - b).

        Examples
        --------
        Modified Power of cosines scattering phase function.
        n = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0} and b=0.5

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for n in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        >>>     pp.semilogy(cos_theta, MPc(n, 0.5)(cos_theta), label='n={}, b=0.5'.format(n))
        >>> pp.legend()
        '''
        super().__init__()
        self._n = float(n)
        self._b = min(max(float(b), 0.0), 1.0)

        self._k1 = self._b*(n + 1.0)/(2**(n + 1))
        self._k2 = (1.0 - self._b)*3.0/2.0

        self._gs = np.zeros((20 + 1,))
        self._gs[1] = self._n/(self._n + 2.0)
        for moment in range(2, 20 + 1):
            self._gs[moment] = self._gs[moment - 1]*(self._n - moment + 1)/\
                               (self._n + moment + 1)
        self._gs *= self._b
        self._gs[0] = 1.0
        self._gs[2] += (1.0 - self._b)*0.4
        self._gs = np.abs(np.maximum(self._gs, 0.0))

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Gegenbauer kernel scattering phase function object.

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
        return self._k1*(1.0 + costheta)**self._n + self._k2*costheta**2

    def g(self, n: int) -> float:
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with
        an analytical solution.
        '''
        if n < self._gs.size:
            g = float(self._gs[n])
        else:
            g = np.prod(np.arange(self._n, self._n - n + 2))/ \
                np.prod(np.arange(self._n + 2, self._n + n + 2))
            g = abs(max(float(g), 0.0))
        return g

    def gs(self, last: int) -> np.ndarray:
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.
        '''
        gs = np.zeros((last + 1,))
        if last + 1 < self._gs.size:
            gs[:last + 1] = self._gs[:last + 1]
        else:
            gs[:self._gs.size] = self._gs
            for moment in range(self._gs.size, last):
                g = gs[moment - 1]*(self._n - moment + 1)/ \
                                   (self._n + moment + 1)
                gs[moment] = abs(max(g, 0.0))
        return gs

    def fastg(self, n: int, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.fastg` method of the base class with
        an analytical solution.
        '''
        return self.g(n)

    def fastgs(self, last: int, **kwargs) -> np.ndarray:
        '''
        Overloads the :py:meth:`PfBase.fastgs` method of the base class with
        an analytical solution.
        '''
        return self.gs(last)

    def __repr__(self):
        return 'MPc({}, {})'.format(self._n, self._b)
