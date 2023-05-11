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
import warnings

import numpy as np

from .pfbase import PfBase
from .miepd import MiePd
from .distribution import Mixture


class MieMixture(PfBase):
    def __init__(self, miepds: Tuple[MiePd],
                 weights: Tuple[float] or np.ndarray):
        '''
        Scattering phase function of a mixture of MiePd scattering
        phase functions.

        Parameters
        ----------
        miepds`: Tuple[MiePd]
            A list of scattering phase functions.
        weights: Tuple[float] or numpy.ndarray
            Weights of the individual scattering phase functions.
            The sum of weights should equal 1.

        Examples
        --------
        Scattering phase functions of a mixture of Normally distributed
        microspherical particles:
        1. mean diameter 0.5 um, standard deviation 15 nm, and weight 0.98,
        2. mean diameter 1.0 um, standard deviation 25 nm, and weight 0.02.

        >>> from matplotlib import pyplot as pp
        >>> import numpy as np
        >>> from xopto import pf
        >>>
        >>> cost_theta = np.linspace(-1.0, 1.0, 1000)
        >>> normal1 = pf.MieNormal(0.5e-6, 15e-9, nsphere=1.6, nmedium=1.33, wavelength=550e-9, nd=1000)
        >>> normal2 = pf.MieNormal(1e-6, 25e-9, nsphere=1.6, nmedium=1.33, wavelength=550e-9, nd=1000)
        >>> mixture = pf.MieMixture((normal1, normal2), (0.98, 0.02))
        >>>
        >>> pp.figure()
        >>> pp.semilogy(cost_theta, mixture(cost_theta), label='Mixture 0.98*N(0.5 um, 15 nm) + 0.02*N(1 um, 25 nm)')
        >>> pp.grid()
        >>> pp.legend()
        >>> pp.show()
        >>>
        '''
        super().__init__()

        if not isinstance(miepds, (tuple, list)):
            miepds = (miepds,)

        if not isinstance(weights, (tuple, list, np.ndarray)):
            weights = (weights,)

        weights_sum = np.sum(weights)
        weights = np.asarray(weights, dtype=np.float64)/weights_sum
        if np.abs(weights_sum - 1.0) > 10*np.finfo(weights.dtype).eps:
            warnings.warn(
                'The sum of scattering phase function weights should be 1!')

        self._miepds = tuple(miepds)
        self._weights = tuple(weights.tolist())

    def distribution(self) -> Mixture:
        '''
        Creates and returns a mixture object that represents the
        number density of all constituents/components as a function
        of diameter.

        Returns
        -------
        mixture: distribution.Mixture
            A mixture object representing the number density of all
            constituents/components as a function of diameter.
        '''
        return Mixture(
            [miepd.distribution() for miepd in self._miepds],
            self._weights
        )

    def miepd(self, index: int or slice) -> MiePd or Tuple[MiePd]:
        '''
        Returns scattering phase function at the specified index.

        Parameters
        ----------
        index: int, slice
            Scattering phase function(s) at the specified index/slice.

        Returns
        -------
        miepd: MiePd
            Scattering phase function(s) at the specified index.
        '''
        return self._miepds[index]

    def weight(self, index: int or slice) -> float or Tuple[float]:
        '''
        Returns weight of scattering phase functions at the specified
        index or slice.

        Parameters
        ----------
        index: int, slice
            Index or slice of the selected scattering phase function weight.

        Returns
        -------
        weight: float or tuple
            The selected scattering phase function weight(s).
        '''
        return self._weights[index]

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the scattering phase function .

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function
            is evaluated.

        Returns
        -------
        pf: float or np.ndarray
            scattering phase function at the specified scattering angle cosines.
        '''
        f = None
        for weight, miepd in zip(self._weights, self._miepds):
            if f is None:
                f = weight*miepd(costheta)*miepd.scs()
            else:
                f += weight*miepd(costheta)*miepd.scs()

        return f/self.scs()

    def pd(self, diameter: float) -> float:
        '''
        Evaluates to the number density value at the specified particle
        diameter.

        Parameters
        ----------
        diameter: float
            Particle diameter (m).

        Returns
        -------
        pd: float
            The value of number density function at the specified
            particle diameter.
        '''
        f = None
        for weight, miepd in zip(self._weights, self._miepds):
            if f is None:
                f = weight*miepd.pd(diameter)
            else:
                f += weight*miepd.pd(diameter)
        return f

    def g(self, n: int, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class.
        Computes the n-th Legendre moment of the scattering phase function.

        Note
        ----
        If n is 1, a precalculated g1 is returned.
        '''
        g = None
        for weight, miepd in zip(self._weights, self._miepds):
            if g is None:
                g = weight*miepd.g(n, **kwargs)*miepd.scs()
            else:
                g += weight*miepd.g(n, **kwargs)*miepd.scs()

        return g/self.scs()

    def fastg(self, n: int, *args, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.fastg` method of the base class.

        Note
        ----
        If n is 1, a precalculated g1 is returned.
        '''
        g = None
        for weight, miepd in zip(self._weights, self._miepds):
            if g is None:
                g = weight*miepd.fastg(n, *args, **kwargs)*miepd.scs()
            else:
                g += weight*miepd.fastg(n, *args, **kwargs)*miepd.scs()

        return g/self.scs()

    def scs(self) -> float:
        '''
        Returns the scattering cross section.
        '''
        scs = None
        for weight, miepd in zip(self._weights, self._miepds):
            if scs is None:
                scs = weight*miepd.scs()
            else:
                scs += weight*miepd.scs()

        return scs

    def ecs(self) -> float:
        '''
        Returns the extinction cross section.
        '''
        ecs = None
        for weight, miepd in zip(self._weights, self._miepds):
            if ecs is None:
                ecs = weight*miepd.ecs()
            else:
                ecs += weight*miepd.ecs()

        return ecs

    def acs(self) -> float:
        '''
        Returns the absorption cross section.
        '''
        acs = None
        for weight, miepd in zip(self._weights, self._miepds):
            if acs is None:
                acs = weight*miepd.acs()
            else:
                acs += weight*miepd.acs()

        return acs

    def __repr__(self):
        return 'MieMixture(miepds={}, weights={})'.format(
            self._miepds, self._weights)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))
