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

from typing import Callable

import numpy as np

from .pfbase import PfBase
from .miepd import MiePd
from .distribution import Normal

class MieNormal(MiePd):
    def __init__(self, center: float, sigma: float,
                 nsphere: float or complex,
                 nmedium: float or complex,
                 wavelength: float,
                 clip: float = 5, nd: int = 100):
        '''
        Scattering phase function of Normally distributed (number density)
        spherical particles:

        .. math::

            p(d) = A e^{\\frac{-(d - center)^2}{2\\sigma^2}}

        The value of parameter :math:`A` is computed so as to normalize the
        integral of the density function on the clipped diameter interval
        :math:`[center - \\sigma clip, center + \\sigma clip]` to 1:

        Parameters
        ----------
        center: float
            Distribution/diameter mean (m).
        sigma: float
            Distribution/diameter standard deviation (m).
        nsphere, nmedium, nd:
            Parameters passed to the :py:meth:`xopto.pf.miepd.MiePd`
            base class constructor.
            See help of :py:class:`xopto.pf.miepd.MiePd` class for more details.
        clip: float
            Distribution/diameter range used to estimate the phase function
            defined as :math:`[center - clip \\sigma, center + clip \\sigma]`.

        Examples
        --------
        Scattering phase functions of Normally distributed microspherical
        particles with a mean diameter of 1 um and a standard deviation of
        100 nm, a mean diameter of 1 um and standard deviation of 50 nm, and
        a mean diameter of 1 um and a standard deviation of 25 nm compared to
        a monodisperse 1 um microspherical particles.

        >>> from matplotlib import pyplot as pp
        >>> import numpy as np
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>> nmie1 = MieNormal(1e-6, 0.1e-6, nsphere=1.6, nmedium=1.33, wavelength=550e-9, nd=1000)
        >>> nmie2 = MieNormal(1e-6, 0.05e-6, nsphere=1.6, nmedium=1.33, wavelength=550e-9, nd=1000)
        >>> nmie3 = MieNormal(1e-6, 0.025e-6, nsphere=1.6, nmedium=1.33, wavelength=550e-9, nd=1000)
        >>> mmie = Mie(nsphere=1.6, nmedium=1.33, diameter=1e-6, wavelength=550e-9)
        >>>
        >>> pp.figure()
        >>> pp.semilogy(costheta, nmie1(cos_theta), label='Normal(1 um, 100 nm)')
        >>> pp.semilogy(costheta, nmie2(cos_theta), label='Normal(1 um, 50 nm)')
        >>> pp.semilogy(costheta, nmie3(cos_theta), label='Normal(1 um, 25 nm)')
        >>> pp.semilogy(costheta, mmie(cos_theta), label='Monodisperse(1 um)')
        >>> pp.legend()
        >>>
        '''
        self._center = float(center)
        self._sigma = float(sigma)
        self._clip = float(clip)
        self._kwargs = {'nsphere':nsphere, 'nmedium':nmedium,
                        'wavelength':wavelength, 'nd':nd}

        pd = Normal(center, sigma, clip)

        super().__init__(pd=pd, drange=pd.range, **self._kwargs)

    def __repr__(self):
        return 'MieNormal(center={}, sigma={}, clip={}, **{})'.format(
            self._center, self._sigma, self._clip, self._kwargs)
