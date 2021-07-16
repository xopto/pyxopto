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

from typing import List, Tuple

import numpy as np
from scipy.integrate import simps, quad
from scattnlay import scattnlay

from .pfbase import PfBase


ComplexVector = float or complex or List[float or complex] or Tuple[float or complex] or np.ndarray
FloatVector = float or List[float] or Tuple[float] or np.ndarray


class MieMl(PfBase):
    def __init__(self, nlayers: ComplexVector, nmedium: ComplexVector, 
                 diameters: List[float], wavelength: float):
        '''
        Constructor of a Mie scattering phase function for multilayered
        spherical particles.

        Parameters
        ----------
        nlayers:  float, complex, list, tuple, numpy.ndarray vector
            Refractive indices (complex) of the concentric spherical layers,
            starting with the refractive index of the innermost layer.
        nmedium: float, complex
            Refractive index (complex) of the surrounding medium.
        diameters: float, list, tuple, numpy.ndarray vector
            Diameters of the concentric spherical layers (m), starting with
            the innermost layer.
        wavelength: float
            Wavelength of light (m).

        Examples
        --------
        Scattering phase function at 550 nm for a hollow SiO2 spherical
        particle with a wall thickness that matches 5% of the outer diameter.
        Refractive index of SiO2 is set to 1.45 and refractive index of water
        to 1.33. The hollow space is given a refractive index of 1.0.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for d in [0.1, 0.5, 1, 2, 5]:
        >>>   pf = MieMl([1.0, 1.45], 1.33, [0.9*d*1e-6, d*1e-6], 550e-9)
        >>>   pp.semilogy(cos_theta, pf(cos_theta), label='Hollow d={:1f} um'.format(d))
        >>> pp.legend()
        '''
        super().__init__()

        if isinstance(nlayers, (int, float, complex)):
            nlayers = (nlayers,)
        if isinstance(diameters, (int, float, complex)):
            diameters = (diameters,)
        if len(nlayers) != len(diameters):
            raise ValueError('The number of values/layers in the nlayers and '\
                             'diameters arguments must be the same!')

        nmedium = np.asarray(nmedium, dtype=np.complex128)
        nlayers = np.asarray(nlayers, dtype=np.complex128)
        diameters = np.asarray(diameters, dtype=np.float)

        wavelength_medium = wavelength/nmedium
        m = nlayers/nmedium
        m.shape = (1, m.size)

        x = np.pi*diameters/wavelength_medium
        x.shape = (1, x.size)

        terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
            x.real, m)

        self._qext = float(Qext)
        self._qsca = float(Qsca)
        self._qabs = float(Qabs)
        self._qbk = float(Qbk)
        self._qpr = float(Qpr)
        self._g1 = float(g)
        self._albedo = float(Albedo)
        self._s1 = S1
        self._S2 = S2

        s = np.pi*(0.5*diameters[-1])**2
        self._sigma_scs = self._qsca*s
        self._sigma_ecs = self._qext*s

        self._x = x
        self._m = m

        self._nmedium = nmedium
        self._nlayers = nlayers
        self._wavelength = wavelength
        self._wavelength_medium = wavelength_medium
        self._diameters = diameters

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
        Overloads the :py:meth:`fastg` method of the base class.
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
        return self._sigma_scs

    def ecs(self) -> float:
        '''
        Returns the extinction cross section.
        '''
        return self._sigma_ecs

    def acs(self) -> float:
        '''
        Returns the absorption cross section.
        '''
        return self._sigma_ecs - self._sigma_scs

    def __call__(self, costheta: float or np.ndarray, scsmul: bool = False) \
            -> float or np.ndarray:
        '''
        Call method of the MieMl scattering phase function.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function is
            evaluated.
        scsmul: bool
            If nonzero, the scattering phase function is multiplied by the
            scattering cross section.

        Returns
        -------
        pf: float or np.ndarray
            Scattering phase function at the specified scattering angle cosines.
        '''
        if isinstance(costheta, float):
            costheta = np.asarray((costheta,), dtype=np.float)

        terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
            self._x.real, self._m, theta=np.arccos(costheta))

        x = np.pi*self._diameters[-1]/self._wavelength_medium.real
        spf = 2*np.pi*(np.abs(S1)**2+np.abs(S2)**2)/ \
            (2*np.pi*x**2.0*Qsca)

        if scsmul:
            return spf*self._sigma_scs
        else:
            return spf

    def __repr__(self):
        return 'MieMl(nlayers={}, nmedium={}, diameters={}, '\
                    'wavelength={})'.format(self._nlayers, self._nmedium,
                                            self._diameters, self._wavelength)
