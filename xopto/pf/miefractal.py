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

from typing import Tuple

import numpy as np

from .pfbase import PfBase
from .miepd import MiePd
from .distribution import Fractal


class MieFractal(MiePd):
    def __init__(self, alpha: float,
                 nsphere: float or complex,
                 nmedium: float or complex,
                 wavelength: float,
                 drange: Tuple[float, float] = (10e-9, 10e-6), nd: int = 1000):
        '''
        Scattering phase function of fractally distributed (number density)
        spherical particles:

        .. math::

            p(d) = A\\frac{1}{d^{\\alpha}}

        The value of parameter :math:`A` is computed so as to normalize the
        integral of the number density function on the specified interval
        drange = :math:`(d_1, d_2)` to 1:

        .. math::

           A = \\frac{d_1^{\\alpha + 1}}{\\alpha\\left(1 - \\left(\\frac{d_1}{d_2}\\right)^{\\alpha + 1}\\right)}

        Parameters
        ----------
        alpha: float
            Parameter alpha of the fractal distribution.
        nsphere, nmedium, wavelength, drange, nd:
            Parameters passed to the :py:meth:`xopto.pf.miepd.MiePd`
            base class constructor.
            See help of :py:class:`xopto.pf.miepd.MiePd` class for more details.

        Examples
        --------
        Scattering phase function of fractally distributed microspherical
        particles with alpha=2.4 and diameter from 10 nm to 10 um.

        >>> from matplotlib import pyplot as pp
        >>> import numpy as np
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>> fmie = MieFractal(alpha=2.4, drange=[10e-9, 10e-6], nsphere=1.6, nmedium=1.33, wavelength=550e-9, nd=1000)
        >>>
        >>> pp.figure()
        >>> pp.semilogy(cos_theta, fmie(cos_theta))
        >>>
        '''
        self._alpha = float(alpha)
        self._kwargs = {'nsphere':nsphere, 'nmedium':nmedium,
                        'wavelength':wavelength, 'drange':drange, 'nd':nd}

        drange = self._kwargs['drange']

        pd = Fractal(self._alpha, drange)

        super().__init__(pd=pd, **self._kwargs)

    def __repr__(self):
        return 'MieFractal(alpha={}, **{})'.format(
            self._alpha, self._kwargs)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))
