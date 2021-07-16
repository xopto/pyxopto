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
from scipy.integrate import quad
from scipy.interpolate import interp1d
from .pfbase import PfBase


class Discrete(PfBase):
    def __init__(self, costheta: np.ndarray, pf: np.ndarray,
                 kind: str = 'cubic', **kwargs):
        '''

        Scattering phase function defined on adiscrete grid of
        :math:`\\cos(\\theta)`. The grid should always include the -1 and +1
        scattering angle cosines!
        A cubic interpolating spline is used to estimate the scatterin
        phase function at an arbitrary :math:`\\cos(\\theta)`.

        Parameters
        ----------
        costheta: np.ndarray
            Discrete :math:`\\cos(\\theta)` points at which the scattering
            phase function values in parameter pf are defined.
        pf: np.ndarray
            The values of the scattering phase function at the costheta points.
        kind: str
            Type of interpolation used to estimate the value of the
            scattering phase function at an arbitrary :math:`\\cos(\\theta)`.
        kwargs: dict
            Optional input arguments passed to the scipy.interpolate.interp1d.

        Examples
        --------
        Discrete representation of a HG scattering phase function for
        :math:`g=0.8`.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>> cos_theta_d = np.linspace(-1.0, 1.0, 50)
        >>> pf_hg = Hg(0.8)
        >>> pf_dhg = Discrete(cos_theta_d, pf_hg(cos_theta_d))
        >>>
        >>> pp.figure()
        >>> pp.semilogy(cos_theta, pf_hg(cos_theta), label='Hg')
        >>> pp.semilogy(cos_theta, pf_dhg(cos_theta), label='Discrete')
        >>> pp.legend()
        '''
        super().__init__()
        self._pf = interp1d(costheta, np.log(pf), kind=kind, **kwargs)
        self._k = 1.0/quad(lambda x: np.exp(self._pf(x)), -1.0, 1.0)[0]
        self._costheta_values = costheta
        self._pf_values = pf
        self._kind = kind
        self._kwargs = kwargs

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the Discrete scattering phase function object.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function
            is evaluated.

        Returns
        -------
        f: float or.ndarray
            Scattering phase function at the specified scattering angle cosines.

        '''
        return np.exp(self._pf(costheta))*self._k

    def __repr__(self):
        return 'Discrete(costheta={}, pf={}, kind={}, **kwargs={})'.format(
            self._costheta_values, self._pf_values, self._kind, self._kwargs)
