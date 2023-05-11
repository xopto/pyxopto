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
from scipy.integrate import quad, simps

from .pfbase import PfBase
from .mie import Mie


class MiePd(PfBase):
    def __init__(self, nsphere: float or complex, nmedium: float or complex,
                 wavelength: float, pd: Callable[[float], float],
                 drange: Tuple[float, float], nd: int = 1000):
        '''
        Mie scattering phase function of an arbitrary size distribution
        of spherical particles over the specified diameter range.

        Parameters
        ----------
        nsphere: float or complex
            Refractive index (complex) of the spherical particles.
        nmedium: float or complex
            Refractive index (complex) of the surrounding medium.
        wavelength: float
            Wavelength of light (m).
        drange: Tuple[float, float]
            Diameter range of the spherical particles as [dmin, dmax], where
            dmin and dmax stand for the minimum and maximum diameters of the
            spherical particles, respectively.
        pd: Callable[[float], float]
            Particle distribution number probability density function. Integral
            of pd over drange should equal 1.0.
        nd: int
            Number of equally spaced control points between dmin and dmax that
            are used to estimate the scattering phase function.
            A fixed-step Simpson numerical integration is used to estimate
            the scattering phase function at the given scattering angle cosines.
            If nd is None, an adaptive-step numerical integration is used (note
            that the computational time might increase dramatically!!!).

        Examples
        --------
        Mie scattering phase function for a normal distribution of
        spherical particles with a mean diameter of 1 um, standard deviation
        of 100 nm, and refractive index of 1.6, for 550 nm monochromatic light
        in water surrounding (refractive index 1.33).
        A 7-sigma diameter range is used for numerical computations.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> dmean = 1e-6    # normal particle distrbution mean
        >>> dsigma = 0.1e-6 # normal particle distribution standard deviation
        >>> pd = lambda d: 1.0/(dsigma*np.sqrt(2*np.pi))*np.exp(-(d - dmean)**2/(2*dsigma**2))
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> miepd = MiePd(1.6, 1.33, 550e-9, pd, [dmean - 7*dsigma, dmean + 7*dsigma])
        >>> pf_values = miepd(cos_theta)
        >>>
        >>> pp.figure()
        >>> pp.semilogy(cos_theta, pf_values)
        '''
        super().__init__()
        self._nd = nd
        self._drange = drange
        self._pd = pd
        self._wavelength = wavelength
        self._nmedium = nmedium
        self._nsphere = nsphere

        if self._drange[0] < 0.0 or self._drange[1] <= 0.0:
            raise ValueError(
                'The range of diameters includes negative values or equals 0!')
 
        if self._nd is None:
            self._scs = quad(
                lambda d: pd(d)*Mie(nsphere, nmedium, d, wavelength).scs(),
                drange[0], drange[1])[0]
            self._ecs = quad(
                lambda d: pd(d)*Mie(nsphere, nmedium, d, wavelength).ecs(),
                drange[0], drange[1])[0]
            self._g1 = quad(
                lambda d: pd(d)*MiePd._g1_scs(nsphere, nmedium, d, wavelength),
                drange[0], drange[1])[0]/self._scs
        else:
            if self._drange is None:
                raise ValueError('Particle diameter range must not be None!')

            self._D = np.linspace(
                float(self._drange[0]), float(self._drange[1]), self._nd)
            self._dd = (self._D[-1] - self._D[0])/(self._D.size - 1)
            self._pdpts = self._pd(self._D)
            self._mie = [None]*self._D.size
            G1_mie = np.zeros((self._D.size,))
            Scs_p = np.zeros((self._D.size,))
            Ecs_p = np.zeros_like(Scs_p)
            for i in range(self._D.size):
                self._mie[i] = mie = Mie(self._nsphere, self._nmedium,
                                         self._D[i], self._wavelength)
                G1_mie[i] = mie.g(1)
                Scs_p[i] = mie.scs()*self._pdpts[i]
                Ecs_p[i] = mie.ecs()*self._pdpts[i]

            self._scs = simps(Scs_p, dx=self._dd)
            self._ecs = simps(Ecs_p, dx=self._dd)
            self._g1 = simps(G1_mie*Scs_p, dx=self._dd)/self._scs

    @staticmethod
    def _g1_scs(nsphere: float or complex, nmedium: float or complex,
                d: float, wavelength: float) -> float:
        mie = Mie(nsphere, nmedium, d, wavelength)
        return mie.scs()*mie.g(1)

    def distribution(self) -> Callable[[float], float]:
        '''
        Returns the underlying distribution object.

        Returns
        -------
        pd: callable(float) -> float
        '''
        return self._pd

    def pd(self, diameter: float) -> float:
        '''
        Evaluates the number density function at the specified particle
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
        return self._pd(diameter)

    def g(self, n: int, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class.
        Computes the n-th Legendre moment of the scattering phase function.

        Note
        ----
        If n is 1, a precalculated g1 is returned.
        '''
        if n == 1:
            return self._g1
        else:
            return super().g(n, **kwargs)

    def fastg(self, n: int, *args, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.fastg` method of the base class.
        Note
        ----
        If n is 1, a precalculated g1 is returned.
        '''
        if n == 1:
            return self._g1
        else:
            return super().fastg(n, *args, **kwargs)

    def scs(self) -> float:
        '''
        Returns the scattering cross section.
        '''
        return self._scs

    def ecs(self) -> float:
        '''
        Returns the extinction cross section.
        '''
        return self._ecs

    def acs(self) -> float:
        '''
        Returns the absorption cross section.
        '''
        return self._ecs - self._scs

    def _mie_pd(self, Pf: np.ndarray) -> np.ndarray:
        self._pdpts.shape = (self._D.size, 1)
        pf = simps(Pf*self._pdpts, dx=self._dd, axis=0)
        return pf

    def _mie_pd_quad(self, costheta: np.ndarray) -> np.ndarray:
        # Multiplication with scs is within Mie.__call__(), the second bool
        #   parameter forces multiplication.
        pf = np.asarray(
            [quad(lambda d: self._pd(d)*Mie(self._nsphere, self._nmedium, d,
                                            self._wavelength)(cost, True),
                  self._drange[0], self._drange[1])[0] for cost in costheta]
        )
        return pf

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the phase function of a spherical particle distribution.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function
            is evaluated.

        Returns
        -------
        pf: float or np.ndarray
            Scattering phase function at the specified deflection angle cosines.
        '''
        costheta = np.array(costheta, copy=False, ndmin=1)
        if self._nd is None:
            return self._mie_pd_quad(costheta)/self._scs
        else:
            costheta = np.asarray(costheta, dtype=np.float64)
            Pf = np.zeros([self._D.size, costheta.size])
            for i in range(self._D.size):
                mie_i = self._mie[i]
                Pf[i] = mie_i(costheta)*mie_i.scs()
            return self._mie_pd(Pf)/self._scs

    def __repr__(self):
        return 'MiePd(nsphere={}, nmedium={}, wavelength={}, pd={}, '\
            'drange={}, nd={})'.format(
                self._nsphere, self._nmedium, self._wavelength, self._pd,
                self._drange, self._nd)
