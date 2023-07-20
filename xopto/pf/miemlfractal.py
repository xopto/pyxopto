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

from typing import Callable, Tuple

import numpy as np
from scipy.integrate import simps, quad

from .distribution import Fractal
from .miemlpd import MieMlPd, ComplexVector, FloatVector


class MieMlFractal(MieMlPd):
    def __init__(self, alpha: float,
                 nlayers: ComplexVector, nmedium: ComplexVector,
                 wavelength: float,
                 drange: Tuple[float, float] = (10e-9, 10e-6),
                 nd: int = 1000, limit: int = None):
        '''
        Scattering phase function of fractally distributed (number density)
        layered spherical particles:

        .. math::

            p(d) = A \\left(\\frac{1}{d}\\right)^\\alpha

        The value of parameter :math:`A` is computed so as to normalize the
        integral of the number density function on the specified interval
        drange = :math:`(d_1, d_2)` to 1:

        .. math::

            A = \\frac{d_1^{\\alpha + 1}}/{\\alpha(1 - \\left(\\frac{d_1}{d_2})^{\\alpha + 1}\\right)}

        Parameters
        ----------
        alpha: float
            Parameter alpha of the fractal distribution.
        nlayers, nmedium, wavelength, drange, nd:
            Parameters passed to the
            :py:meth:`xopto.pf.miemlpd.MieMlPd` base class
            constructor. See help of :py:class:`xopto.pf.miemlpd.MieMlPd`
            class for more details.

        Note
        ----
        If the medium or particle have a nonzero absorption coefficient, the
        refractive index becomes complex :math:`n + ik`, where :math:`k` is
        related to the absorption coefficient :math:`\\mu_{a}` as
        :math:`\\mu_{a} = 4 \\pi k / \\lambda_0`, where :math:`\\lambda_0`
        is the wavelength of light in vacuum.

        Examples
        --------
        Scattering phase function of fractally distributed hollow spherical
        particles with outer diameter from 10 nm to 10 um and parameter
        alpha=2.4.
        The wall thickness of spherical particles accounts for 5% of the
        outer particle diameter.

        >>> from matplotlib import pyplot as pp
        >>> import numpy as np
        >>>
        >>> cost_heta = np.linspace(-1.0, 1.0, 1000)
        >>> nlayers=[1.0, 1.45]
        >>> fmieml = MieMlFractal(alpha=2.4, drange=[10e-9, 10e-6], nlayers=nlayers, nmedium=1.33, wavelength=550e-9, nd=1000)
        >>>
        >>> pp.figure()
        >>> pp.semilogy(cost_heta, fmieml(cost_heta))
        >>>
        '''
        self._alpha = float(alpha)
        self._kwargs = {'nlayers':nlayers, 'nmedium':nmedium,
                        'wavelength':wavelength, 'drange':drange, 'nd':nd,
                        'limit':limit}

        drange = self._kwargs['drange']

        pd = Fractal(self._alpha, drange)

        super().__init__(pd=pd, **self._kwargs)

    def __repr__(self):
        return 'MieMlFractal(alpha={}, **{})'.format(
            self._alpha, self._kwargs)
