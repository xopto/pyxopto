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


class Hg(PfBase):
    def __init__(self, g: float):
        '''
        Henyey-Greenstein scattering phase function constructor.

        Parameters
        ----------
        g: float
            Anisotropy factor.

        Examples
        --------
        Henyey-Greenstein scattering phase function for anisotropy factors
        g = {0, 0.3 0.5, 0.8, 0.9, 0.95}.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for g in [0, 0.3, 0.5, 0.8, 0.9, 0.95]:
        >>>     pp.semilogy(cos_theta, Hg(g)(cos_theta), label='g={}'.format(g))
        >>> pp.legend()
        '''

        super().__init__()
        eps = np.finfo(np.float).eps
        self._g = max(min(g, 1.0 - eps), -1.0 + eps)

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Henyey-Greenstein scattering phase function object.

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
        return 0.5*(1.0 - self._g*self._g)/ \
            (1 + self._g*self._g - 2*self._g*costheta)**1.5

    def g(self, n):
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with 
        an analytical solution.
        '''
        return self._g**n

    def gs(self, last):
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.
        '''
        return self._g**np.arange(last + 1, dtype=np.float)

    def fastg(self, n, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with
        an analytical solution.
        '''
        return self._g**n

    def fastgs(self, last, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.
        '''
        return self._g**np.arange(last + 1, dtype=np.float)

    def __repr__(self):
        return 'Hg({})'.format(self._g)
