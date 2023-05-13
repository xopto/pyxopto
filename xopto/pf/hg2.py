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

class Hg2(PfBase):
    def __init__(self, g1: float, g2: float, b: float):
        '''
        Double Henyey-Greenstein (HG) scattering phase function constructor.

        Parameters
        ----------
        g1: float
            Anisotropy factor of the first HG scattering phase function. Must
            be from interval [0, 1].
        g1: float
            Anisotropy factor of the second HG scattering phase function. Must
            be from interval [-1, 0]
        b: float
            Relative contribution of the second HG to the total scattering phase
            function. Must be a value from [0, 1].

        Examples
        --------
        Double HG scattering phase function for anisotropy factors
        g = {0, 0.3 0.5, 0.8, 0.9, 0.95}, where g1=g, g2=-g1 and b=0.5.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for g in [0, 0.3, 0.5, 0.8, 0.9, 0.95]:
        >>>     pp.semilogy(cos_theta, Hg2(g, -g, 0.5)(cos_theta), label='g={}'.format(g))
        >>> pp.legend()
        '''

        super().__init__()
        if g1 < 0.0 or g1 > 1.0:
            raise ValueError(
                'Anisotropy (g1) of the fist HG must be from [0, 1]!')
        if g2 < -1.0 or g2 > 0.0:
            raise ValueError(
                'Anisotropy (g2) of the second HG must be from [-1, 0]!')

        self._hg1 = Hg(g1)
        self._hg2 = Hg(g2)
        self._b = min(max(float(b), 0.0), 1.0)

        def pf(costheta):
            return (1.0 - self._b)*self._hg1(costheta) + \
                   self._b*self._hg2(costheta)

        self._pf = pf
    
    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the double Henyey-Greenstein scattering phase function
        object.

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

    def g(self, n):
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with 
        an analytical solution.
        '''
        return (1.0 - self._b)*self._hg1._g**n + self._b*self._hg2._g**n

    def gs(self, last):
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.
        '''
        p = np.arange(last + 1, dtype=np.float64)
        return (1.0 - self._b)*self._hg1._g**p + self._b*self._hg2._g**p

    def fastg(self, n, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with
        an analytical solution.
        '''
        return (1.0 - self._b)*self._hg1._g**n + self._b*self._hg2._g**n

    def fastgs(self, last, **kwargs):
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.
        '''
        rp = np.arange(last + 1, dtype=np.float64)
        return (1.0 - self._b)*self._hg1._g**p + self._b*self._hg2._g**p

    def __repr__(self):
        return 'Hg2({}, {}, {})'.format(self._hg1.g1, self._hg2.g2, self._b)
