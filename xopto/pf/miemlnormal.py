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

from .distribution import Normal
from .miemlpd import MieMlPd, ComplexVector, FloatVector


class MieMlNormal(MieMlPd):
    def __init__(self, center: float, sigma: float,
                 nlayers: ComplexVector, nmedium: ComplexVector,
                 diameters: FloatVector, wavelength: float,
                 clip: float = 5, nd: int = 100, limit: int = None):
        '''
        Scattering phase function of Normally distributed (number density)
        layered spherical particles:

            pd = A e^{\\frac{-(d - center)^2}{(2\\sigma^2}}

        The value of parameter A is computed so as to normalize the integral
        of the density function on the clipped diameter interval
        :math:`[center - \\sigma clip, center + \\sigma clip]` to 1:

        Parameters
        ----------
        center: float
            Distribution/diameter mean (m).
        sigma: float
            Distribution/diameter standard deviation (m).
        nlayers, nmedium, diameters, nd:
            Parameters passed to the MieMlPd base class constructor.
            See help of MieMlPd class for more details.
        clip: float
            Distribution/diameter range used to estimate the phase function
            defined as :math:`[center - clip \\sigma, center + clip \\sigma]`.

        Note
        ----
        If the medium or particle have a nonzero absorption coefficient, the
        refractive index becomes complex :math:`n + ik`, where :math:`k` is
        related to the absorption coefficient :math:`\\mu_{a}` as
        :math:`\\mu_{a} = 4 \\pi k / \\lambda_0`, where :math:`\\lambda_0`
        is the wavelength of light in vacuum.

        Examples
        --------
        Scattering phase functions of Normally distributed hollow spherical
        particles with mean diameter 1 um and standard deviation 100 nm,
        mean diameter 1 um and standard deviation 50 nm, and
        mean diameter 1 um and standard deviation 25 nm compared to
        a monodisperse 1 um microspherical particles. The wall thickness of
        spherical particles accounts for 5% of the outer particle diameter.

        >>> from matplotlib import pyplot as pp
        >>> import numpy as np
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>> nlayers=[1.0, 1.45]
        >>> diameters=[0.9e-6, 1e-6]
        >>> nmie1 = MieMlNormal(1e-6, 0.1e-6, nlayers=nlayers, nmedium=1.33, diameters=diameters, wavelength=550e-9, nd=1000)
        >>> nmie2 = MieMlNormal(1e-6, 0.05e-6, nlayers=nlayers, nmedium=1.33, diameters=diameters, wavelength=550e-9, nd=1000)
        >>> nmie3 = MieMlNormal(1e-6, 0.025e-6, nlayers=nlayers, nmedium=1.33, diameters=diameters, wavelength=550e-9, nd=1000)
        >>> mmie = MieMl(nlayers=nlayers, nmedium=1.33, diameters=diameters, wavelength=550e-9)
        >>>
        >>> pp.figure()
        >>> pp.semilogy(cos_theta, nmie1(cos_theta), label='Normal(1 um, 100 nm)')
        >>> pp.semilogy(cos_theta, nmie2(cos_theta), label='Normal(1 um, 50 nm)')
        >>> pp.semilogy(cos_theta, nmie3(cos_theta), label='Normal(1 um, 25 nm)')
        >>> pp.semilogy(cos_theta, mmie(cos_theta), label='Monodisperse(1 um)')
        >>> pp.legend()
        >>>
        '''
        self._center = float(center)
        self._sigma = float(sigma)
        self._clip = float(clip)
        self._kwargs = {'nlayers':nlayers, 'nmedium':nmedium,
                        'diameters':diameters,
                        'wavelength':wavelength, 'nd':nd,
                        'limit':limit}

        pd = Normal(center, sigma, clip)

        super().__init__(pd=pd, drange=pd.range, **self._kwargs)

    def __repr__(self):
        return 'MieMlNormal(center={}, sigma={}, clip={}, **{})'.format(
            self._center, self._sigma, self._clip, self._kwargs)


