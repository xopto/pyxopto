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

from .pfbase import PfBase
from .mieml import MieMl, ComplexVector, FloatVector


class MieMlPd(PfBase):
    @staticmethod
    def scale(diameters: np.ndarray, diameter: float) -> np.ndarray:
        '''
        Multiplicatively scale the diameters of the layer stack to
        match the outer diameter of the stack with the specified diameter.

        Parameters
        ----------
        diameters: np.ndarray
            Diameters of the layer stack.
        diameter: float
            Target outer diameter.

        Returns
        -------
        scaled_diameters: np.ndarray
            Diameters of the layer stack scaled to the specified outer diameter.
        '''
        return np.asarray(diameters, dtype=np.float64)*(diameter/diameters[-1])
    
    def __init__(self, nlayers: ComplexVector, nmedium: ComplexVector,
                 diameters: FloatVector, wavelength: float,
                 pd: Callable[[float], float], drange: Tuple[float, float],
                 nd: int = 1000,
                 dscalefun: Callable[[np.ndarray, float], np.ndarray] = None,
                 limit: int = None):
        '''
        Mie scattering phase function for layered spherical particles that
        follow the given size probability density function within the given
        range of outer diameters.
        By default, the layer stack is multiplicatively scaled to the current
        value of the outer diameter.

        Parameters
        ----------
        nlayers: float, complex, list, tuple, numpy.ndarray vector
            Refractive indices (complex) of the concentric spherical layers,
            starting with the refractive index of the innermost layer.
        nmedium: float, complex
            Refractive index (complex) of the surrounding medium.
        diameters: float, list, tuple, numpy.ndarray vector
            Diameters of the concentric spherical layers (m), starting with
            the innermost layer. The layer stack is scaled according
            to the current diameter of the distribution. A custom scaling
            function can be passed in the scalefun parameter.
            A multiplicative scaling of the layer stack to the current
            outer diameter given by the distribution is used by default.
        wavelength: float
            Wavelength of light (m).
        pd: Callable[[float], float]
            Particle number density distribution function. Integral
            of pd over diameter range given in drange should equal 1.0.
        drange: [float, float] or (float, float)
            Layered spherical particle outer diameter range as [dmin, dmax],
            where dmin and dmax stand for the minimum and maximum diameter of
            spherical particles, respectively.
        nd: int
            Number of equally spaced control points between dmin and dmax that
            are used to estimate the phase function. A fixed-step Simpson
            numerical integration is used to estimate the phase function at
            the given deflection angle cosines. If nd is None, adaptive-step
            numerical integration is used (note that the computational time
            might increase dramatically!!!).
        scalefun: Callable[[np.ndarray, float], float]
            A function that scales the layer stack to the current outer diameter
            of the distribution.
            If None (default), multiplicative scaling is used to transform
            the layer stack so that the outermost layer diameter matches the
            current distribution diameter.
        limit: int
            An upper bound on the number of subintervals used in the quad
            adaptive algorithm. Set to None for default value.

        Note
        ----
        If the medium or particle have a nonzero absorption coefficient, the
        refractive index becomes complex :math:`n + ik`, where :math:`k` is
        related to the absorption coefficient :math:`\\mu_{a}` as
        :math:`\\mu_{a} = 4 \\pi k / \\lambda_0`, where :math:`\\lambda_0`
        is the wavelength of light in vacuum.

        Examples
        --------
        Mie scattering phase function for a normal distribution of
        hollow spherical particles with a mean outer diameter of 1 um, standard
        deviation of 100 nm, a wall thickness that equals 5% of the outer
        diameter, wall refractive index of 1.45, for 550 nm monochromatic light
        (vacuum) suspended in water (refractive index 1.33). A 7-sigma
        diameter range is used for numerical computation.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>> from xopto.pf.distribution import Normal
        >>>
        >>> dmean = 1e-6    # normal particle distrbution mean
        >>> dsigma = 0.1e-6 # normal particle distribution standard deviation
        >>> pd = Normal(dmean, dsigma)
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> miemlpd = MieMlPd([1.0, 1.45], 1.33, [0.9e-6, 1.0e-6], 550e-9, pd, pd.range)
        >>> pf_values = miemlpd(cos_theta)
        >>>
        >>> pp.figure()
        >>> pp.semilogy(cos_theta, pf_values)
        '''
        super().__init__()

        if isinstance(nlayers, (int, float, complex)):
            nlayers = (nlayers,)
        if isinstance(diameters, (int, float, complex)):
            diameters = (diameters,)
        if len(nlayers) != len(diameters):
            raise ValueError('The number of values/layers in the nlayers and '\
                             'diameters arguments must be the same!')
        if limit is None:
            limit = 50

        if dscalefun is None:
            dscalefun = MieMlPd.scale

        nmedium = np.asarray(nmedium, dtype=np.complex128)
        nlayers = np.asarray(nlayers, dtype=np.complex128)
        diameters = np.asarray(diameters, dtype=np.float64)

        self._nd = nd
        self._drange = drange
        self._pd = pd
        self._wavelength = wavelength
        self._nmedium = nmedium
        self._nlayers = nlayers
        self._diameters = diameters
        self._dscalefun = dscalefun

        if self._nd is None:
            self._scs = quad(
                lambda d: pd(d)*MieMl(
                    nlayers, nmedium, dscalefun(diameters, d),
                    wavelength).scs(), drange[0], drange[1], limit=limit)[0]
            self._ecs = quad(
                lambda d: pd(d)*MieMl(
                    nlayers, nmedium, dscalefun(diameters, d),
                    wavelength).ecs(), drange[0], drange[1], limit=limit)[0]
            self._g1 = quad(
                lambda d: pd(d)*MieMlPd._g1_scs(
                    nlayers, nmedium, dscalefun(diameters, d), wavelength),
                    drange[0], drange[1], limit=limit)[0]/self._scs
        else:
            if self._drange is None:
                raise ValueError('Particle diameter range must not be None!')

            self._D = np.linspace(
                float(self._drange[0]), float(self._drange[1]), self._nd)
            self._dd = (self._D[-1] - self._D[0])/(self._D.size - 1)
            self._pdpts = self._pd(self._D)
            self._mieml = [None]*self._D.size
            G1_mieml = np.zeros((self._D.size,))
            Scs_p = np.zeros((self._D.size,))
            Ecs_p = np.zeros_like(Scs_p)
            for i in range(self._D.size):
                self._mieml[i] = mieml = MieMl(nlayers, nmedium,
                                               dscalefun(diameters, self._D[i]),
                                               self._wavelength)
                G1_mieml[i] = mieml.g(1)
                Scs_p[i] = mieml.scs()*self._pdpts[i]
                Ecs_p[i] = mieml.ecs()*self._pdpts[i]

            self._scs = simps(Scs_p, dx=self._dd)
            self._ecs = simps(Ecs_p, dx=self._dd)
            self._g1 = simps(G1_mieml*Scs_p, dx=self._dd)/self._scs

    @staticmethod
    def _g1_scs(nlayers, nmedium: ComplexVector, diameters: FloatVector,
                wavelength: float) -> float:
        mieml = MieMl(nlayers, nmedium, diameters, wavelength)
        return mieml.scs()*mieml.g(1)

    def distribution(self) -> Callable[[float], float]:
        '''
        Returns the underlaying size distribution object.

        Returns
        -------
        pd: Callable[[float], float]
            A callable size distribution that takes the size
            parameter/diameter and returns the corresponding probability
            density.
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

    def g(self, n: float, **kwargs) -> float:
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

    def fastg(self, n: float, *args, **kwargs) -> float:
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

    def _mieml_pd(self, Pf):
        self._pdpts.shape = (self._D.size, 1)
        pf = simps(Pf*self._pdpts, dx=self._dd, axis=0)
        return pf

    def _mieml_quad_pd(self, costheta):
        # Multiplication with scs is within MieMl.__call__(), the second bool
        #   parameter forces multiplication.
        pf = np.asarray(
            [quad(lambda d: self._pd(d)*MieMl(
                self._nlayers, self._nmedium, self._dscalefun(self._diameters, d),
                self._wavelength)(cost, True),
                self._drange[0], self._drange[1])[0] for cost in costheta]
        )
        return pf

    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the scattering phase function.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function is evaluated.

        Returns
        -------
        pf: float or np.ndarray
            The Scattering phase function at the specified scattering
            angle cosines.
        '''
        costheta = np.array(costheta, copy=False, ndmin=1)
        if self._nd is None:
            return self._mieml_quad_pd(costheta)/self._scs
        else:
            costheta = np.asarray(costheta, dtype=np.float64)
            Pf = np.zeros([self._D.size, costheta.size])
            for i in range(self._D.size):
                mie_i = self._mieml[i]
                Pf[i] = mie_i(costheta)*mie_i.scs()
            return self._mieml_pd(Pf)/self._scs

    def __repr__(self):
        return 'PdMie(nlayers={}, nmedium={}, diameters={}, wavelength={}, '\
               'pd={}, drange={}, nd={}, dscalefun={})'.format(
                   self._nlayers, self._nmedium, self._diameters,
                   self._wavelength, self._pd, self._drange,
                   self._nd, self._dscalefun)
