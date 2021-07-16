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


class Rayleigh(PfBase):
    def __init__(self, gamma: float):
        '''
        Rayleigh scattering phase function constructor.

        Parameters
        ----------
        gamma: float
            Molecular anisotropy. Scattering is isotropic for gamma = 0.

        Examples
        --------
        Rayleigh scattering phase function for gamma [0.0, 0.1, 0.2, 0.5, 1.0].

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for gamma in [0.0, 0.1, 0.2, 0.5, 1.0]:
        >>>     pp.semilogy(cos_theta, Rayleigh(gamma)(cos_theta), label='gamma={}'.format(gamma))
        >>> pp.legend()
        '''
        super().__init__()
        self._gamma = float(gamma)
        self._gs_data = np.array((1.0, 0.0, super().g(2)))

    def g(self, n: int) -> float:
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class with
        an analytical solution.
        '''
        if n <= 2:
            return self._gs_data[n]
        else:
            return 0

    def gs(self, last: int) -> np.ndarray:
        '''
        Overloads the :py:meth:`PfBase.gs` method of the base class with
        an analytical solution.
        '''
        gs = np.zeros((last + 1,))
        last = int(last)
        gs[:min(3, last)] = self._gs_data[:min(3, last)]
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

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Rayleigh scattering phase function object.

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
        return (3.0/8.0)*((1.0 + 3.0*self._gamma) + \
                (1.0 - self._gamma)*costheta*costheta)/(1.0 + 2.0*self._gamma)

    def __repr__(self):
        return 'Rayleigh({})'.format(self._gamma)
