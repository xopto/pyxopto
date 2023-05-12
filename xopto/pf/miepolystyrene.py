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

from .mie import Mie
from .miefractal import MieFractal
from .mienormal import MieNormal

from xopto.materials import ri

class MiePolystyrene(Mie):
    def __init__(self, diameter: float, wavelength: float,
                ripolystyrene: float or complex = None,
                riliquid: float or complex = None):
        '''
        Scattering phase function of monodisperse polystyrene microspherical
        particles.

        Parameters
        ----------
        diameter: float
            Diameter of the Microspherical particle (m).

        wavelength: float or complex
            Wavelength of light (m).

        ripolystyrene: float or complex
            Refractive index of polystyrene at the specified wavelength.
            If None, the builtin polystyrene refractive index is used.

        riliquid: float
            Refrective index of the liquid phase (water by default).
            If None, the builtin refractive index of water is used.

        Note
        ----
        If the medium or particle have a nonzero absorption coefficient, the
        refractive index becomes complex :math:`n + ik`, where :math:`k` is
        related to the absorption coefficient :math:`\\mu_{a}` as
        :math:`\\mu_{a} = 4 \\pi k / \\lambda_0`, where :math:`\\lambda_0`
        is the wavelength of light in vacuum.
        '''
        if ripolystyrene is None:
            ripolystyrene = ri.polystyrene.default(wavelength)

        if riliquid is None:
            riliquid = ri.water.default(wavelength)

        super().__init__(ripolystyrene, riliquid, diameter, wavelength)

    def __repr__(self):
        return 'MiePolystyrene(diameter={}, wavelength={})'\
            .format(self._diameter, self._wavelength)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))


class MieNormalPolystyrene(MieNormal):
    def __init__(self, center: float, sigma:float, wavelength: float,
                 ripolystyrene: float or complex = None,
                 riliquid: float or complex = None,
                 clip: float = 5,nd: int = 100):
        '''
        Fractal scattering phase function for spherical polystyrene
        microparticles.

        Parameters
        ----------
        center: float
            Distribution/diameter mean (m).
        sigma: float
            Distribution/diameter standard deviation (m).
        wavelength: float
            Wavelength of light (m).
        ripolystyrene: float
            Refrective index of polystyrene at the specified wavelength.
            If None, the builtin polystyrene refractive index is used.
        riliquid: float
            Refrective index of the liquid phase (water by default).
            If None, the builtin refractive index of water is used.
        clip: float
            Distribution/diameter range used to estimate the phase function
            defined as [center - clip*sigma, center + clip*sigma].
        nd: int
            Number of equally spaced control points between dmin and dmax that
            are used to estimate the phase function. A fixed-step Simpson
            numerical integration is used to estimate the phase function at
            the given deflection angle cosines. If nd is None, adaptive-step
            numerical integration is used (note that the computational time
            might increase dramatically!!!).
        '''
        if ripolystyrene is None:
            ripolystyrene = ri.polystyrene.default(wavelength)

        if riliquid is None:
            riliquid = ri.water.default(wavelength)

        super().__init__(cnter=center, sigma=sigma,
                         nsphere=ripolystyrene,
                         nmedium=riliquid,
                         wavelength=wavelength,
                         clip=clip, nd=nd)

    def __repr__(self):
        return 'MieNormalPolystyrene(center={}, sigma={}, wavelength={})'\
            .format(self._center, self._sigma, self._wavelength)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))


class MieFractalPolystyrene(MieFractal):
    def __init__(self, alpha: float, wavelength: float,
                 drange: Tuple[float, float] = None,
                 ripolystyrene: float or complex = None,
                 riliquid: float or complex = None, nd: int = 1000):
        '''
        Fractal scattering phase function for spherical polystyrene
        microparticles.

        Parameters
        ----------
        alpha: float
            Parameter alpha of the fractal distribution.

        wavelength: float
            Wavelength (m).

        drange: list, tuple of two float
            Finite range (m) of the particle diameter used to compute the
            scattering phase function given as (dmin, dmax). If None, default
            a range [5e-9, 30e-6] m is used.

        ripolystyrene: float
            Refrective index of polystyrene at the specified wavelength.
            If None, the builtin polystyrene refractive index is used.

        riliquid: float
            Refrective index of the liquid phase (water by default).
            If None, the builtin refractive index of water is used.

        nd: int
            Number of equally spaced control points between dmin and dmax that
            are used to estimate the phase function. A fixed-step Simpson
            numerical integration is used to estimate the phase function at
            the given deflection angle cosines. If nd is None, adaptive-step
            numerical integration is used (note that the computational time
            might increase dramatically!!!).
        '''
        # drange = np.array([10e-9, 20e-6])  # Calabro, JBO, 2014
        if drange is None:
            drange = [5e-9, 30e-6]  # Naglic, OL, 2017

        if ripolystyrene is None:
            ripolystyrene = ri.polystyrene.default(wavelength)

        if riliquid is None:
            riliquid = ri.water.default(wavelength)

        super().__init__(alpha=alpha, drange=drange,
                         nsphere=ripolystyrene,
                         nmedium=riliquid,
                         wavelength=wavelength, nd=nd)

    def __repr__(self):
        return 'MieFractalPolystyrene(alpha={}, wavelength={})'\
            .format(self._alpha, self._wavelength)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))
