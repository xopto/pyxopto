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
from .hg import Hg

class MHg(PfBase):
    def __init__(self, g: float, b: float):
        '''
        Modified Henyey-Greenstein scattering phase function.

        Parameters
        ----------
        g: float
            Anisotropy factor.
        b: float
            Fractional contribution of the Henyey-Greenstein scattering
            phase function.

        Examples
        --------
        Modified Henyey-Greenstein scattering phase function fo anisotropy
        factors g = {0, 0.3 0.5, 0.8, 0.9, 0.95} and b=0.5.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for g in [0.0, 0.3, 0.5, 0.8, 0.9, 0.95]:
        >>>     pp.semilogy(cos_heta, MHg(g, 0.5)(cos_theta), label='b=0.5, g={}'.format(g))
        >>> pp.legend()
        '''
        super().__init__()
        self._b = np.clip(float(b), 0.0, 1.0)
        self._hg = Hg(g)
        self._rayleigh_gs_data = np.array((1.0, 0.0, 0.4))

    def _rayleigh_g(self, n: int) -> float:
        if n <= 2:
            return self._rayleigh_gs_data[n]
        else:
            return 0.0

    def _rayleigh_gs(self, n: int) -> np.ndarray:
        gs = np.zeros((n + 1,))
        n = int(n)
        gs[:min(n, 3)] = self._rayleigh_gs_data[:min(n, 3)]
        return gs

    def g(self, n: int) -> float:
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with
        an analytical solution.
        '''
        return self._b*self._hg.g(n) + (1.0 - self._b)*self._rayleigh_g(n)

    def gs(self, last: int) -> np.ndarray:
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.
        '''
        return self._b*self._hg.gs(last) + \
               (1.0 - self._b)*self._rayleigh_gs(last)

    def fastg(self, n: int, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.fastg` method of the base class
        with an analytical solution.
        '''
        return self.g(n)

    def fastgs(self, last: int, **kwargs) -> np.ndarray:
        '''
        Overloads the :py:meth:`PfBase.fastgs` method of the base class with
        an analytical
        olution.
        '''
        return self.gs(last)

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Modified Henyey-Greenstein scattering phase function.

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
        return self._b*self._hg(costheta) + \
            (1.0 - self._b)*3.0/(2.0)*costheta**2

    def __repr__(self):
        return 'MHg({}, {})'.format(self._hg._g, self._b)
