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

        self._precalculated_gs = self._precalculate_gs()

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

    def fastg(self, n: int, *args, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.fastg` method of the base class.

        Note
        ----
            Using analytical solution for n = 0, 1, 2 or 3.
        '''
        g = None
        if n <= 3:
            g = self._precalculated_gs[n]
        if g is None:
            g = PfBase.fastg(self, n, *args, **kwargs)

        return g

    def g(self, n, *args, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with
        an analytical solution.
        '''
        g = None
        if n <= 3:
            g = self._precalculated_gs[n]
        if g is None:
            g = PfBase.g(self, n, *args, **kwargs)

        return g

    def _precalculate_gs(self):
        g, a = self._gg, self._a

        g1 = g2 = g3 = None
        if a != 1 and  a != 0 and g != 0:
            g_g = g*g

            t1, t2 = (1.0 + g)**(2.0*a), (1.0 - g)**(2.0*a)
            if t1 == t2:
                print(g, a)
                raise ValueError('Bad')

            L = (t1 + t2)/(t1 - t2)

            g1 = (2.0*g*a*L - (1.0 + g_g))/(2.0*g*(a - 1.0))

            if a != 2:
                g2 = 3.0*(1.0 + g_g)*g1/(2.0*g*(2.0 - a)) - (1.0 + a)/(2.0 - a)

                if a != 3:
                    g3 = 5.0/(2.0*(3.0 - a))*((1.0 + g_g)*g2/g - \
                        a*L + (1.0 + g_g)/(2.0*g)) - 3.0/2.0*g1

        return 1.0, g1, g2, g3

