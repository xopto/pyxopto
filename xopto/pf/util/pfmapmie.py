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
import os.path

import numpy as np

from .pfmapbase import PfMap2DBase
from xopto.pf.miepolystyrene import MieFractalPolystyrene, MiePolystyrene
from xopto import DATA_PATH


class MieFractalPolystyreneMap(PfMap2DBase):
    DEFAULT_MAP_FILE = 'fractal_mie_polystyrene_map.npz'
    XLABEL = '$\\alpha$'
    YLABEL = '$\\lambda (nm)$'
    PLOTSCALEFACTORX = 1.0
    PLOTSCALEFACTORY = 1e9

    @classmethod
    def precalculate(cls, n: int = 100, filename: str = None,
                     verbose: bool = True):
        '''
        Precalculate fractal Mie polystyrene scattering phase function lookup table.

        Parameters
        ----------
        n: int
            Number of steps along the scattering phase function parameters.
        filename: str
            Output file or None to save as the default lookup table.
        verbose: bool
            Turn on verbose progress report.
        '''
        if verbose:
            print('\nCreating MieFractalPolystyreneMap map:')

        fracmiepolymap = \
            MieFractalPolystyreneMap(alpha=np.linspace(2, 6, 30),
                                    wavelength=np.linspace(0.4e-6, 1e-6, 60),
                                    ng=15, ncostheta=5000, nd=1000)

        if filename is None:
            filename = cls.default_data_file()

        fracmiepolymap.save(filename=filename)

    def __init__(self, alpha: np.ndarray = None, wavelength: np.ndarray = None,
                 ng: int = 15, drange: Tuple[float, float] = (5e-9, 30e-6),
                 nd: int = 1000, filename: str = None, ncostheta: int = None):
        '''
        Prepares maps of the first ng Legendre moments of the
        MieFractalPolystyrene (FMIE) phase function, over the specified
        range of the FMIE parameters alpha and wavelength. If ng >= 2, a map of
        gamma is preapred and if ng >= 3 maps of delta and sigma are
        prepared as well. The maps are used to obtain an initial estimate
        when calculateng the FMIE phase function parameters from given
        Legendre moments, gamma and/or delta.

        Parameters
        ----------
        alpha: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the alpha parameter of the Fractal Mie phase function.
        wavelength: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the a wavelength parameter of Fractal Mie phase function [m].
        drange: list, tuple of two float
            Finite range (m) of the particle diameter used to compute the
            scattering phase function given as (dmin, dmax). If None, default
            a range [5e-9, 30e-6] m is used.
        nd: int
            Number of equally spaced control points between dmin and dmax that
            are used to estimate the phase function. A fixed-step Simpson
            numerical integration is used to estimate the phase function at
            the given deflection angle cosines. If nd is None, adaptive-step
            numerical integration is used (note that the computational time
            might increase dramatically!!!).
        ng: int
            Maps are created for the first ng Legendre moments. If ng >= 2,
            a map of gamma is prepared and if g >= 3 maps of delta and sigma
            are prepared as well.
        filename: str
            File with saved data. The values of all the other parameters are
            ignored and restored from the file.
        ncostheta: int
            Number of nodes used to compute the Legendre moments. Use a
            large number (> 1000) for accurate results. If None, adaptive
            step integration is used, which is accurate but can become slow.

        Note
        ----
        The value of parameter ng should be >> 3 to accurately
        estimate the value of parameter sigma.

        Examples
        --------
        Prepares maps of gamma and delta, and estimates the FMIE parameters
        given a) g and gamma are known b) gamma and delta are known.

        >>> import numpy as np
        >>>
        >>> m = MieFractalPolystyreneMap(np.linspace(2, 6, 30), np.linspace(0.4e-6, 1e-6, 60), ng=3)
        >>> alpha, wavelength = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'alpha:', alpha, 'wavelength:', wavelength)
        >>> alpha, wavelength = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'alpha:', alpha, 'wavelength:', wavelength)
        >>>

        Load maps from the default file included in the data/pf folder.

        >>> m = MieFractalPolystyreneMap.fromfile()
        >>> alpha, wavelength = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'alpha:', alpha, 'wavelength:', wavelength)
        >>> alpha, wavelength = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'alpha:', alpha, 'wavelength:', wavelength)
        >>>
        '''
        if alpha is None:
            alpha = np.linspace(2.0, 6.0, 30)

        if wavelength is None:
            np.linspace(0.4e-6, 1e-6, 60)

        super().__init__(param1=alpha, param2=wavelength,
                         ng=ng, pf=MieFractalPolystyrene,
                         filename=filename, ncostheta=ncostheta,
                         drange=drange, nd=nd)

    def alpha(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of alpha parameter.
        '''
        return self.param1()

    def wavelength(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of wavelengths.
        '''
        return self.param2()

    def alpha_grid(self) -> np.ndarray:
        '''
        Returns a 2D map (meshgrid) of the first scattering phase function
        parameter alpha.
        '''
        return self.grid1()

    def wavelength_grid(self) -> np.ndarray:
        '''
        Returns a 2D map (meshgrid) of the second parameter of the
        scatterin phase function wavelength.
        '''
        return self.grid2()

class MiePolystyreneMap(PfMap2DBase):
    DEFAULT_MAP_FILE = 'mie_polystyrene_map.npz'
    XLABEL = '$d (nm)$'
    YLABEL = '$\\lambda (nm)$'
    PLOTSCALEFACTORX = 1e9
    PLOTSCALEFACTORY = 1e9

    @classmethod
    def precalculate(cls, n: int = 100, filename: str = None,
                     verbose: bool = True):
        '''
        Precalculate monodisperse Mie polystyrene scattering phase function
        lookup table.

        Parameters
        ----------
        n: int
            Number of steps along the scattering phase function parameters.
        filename: str
            Output file or None to save as the default lookup table.
        verbose: bool
            Turn on verbose progress report.
        '''
        if verbose:
            print('\nCreating MiePolystyrene map:')

        miepolymap = MiePolystyreneMap(diameter=np.linspace(0.2e-6, 5e-6, n),
                                    wavelength=np.linspace(0.4e-6, 1e-6, n),
                                    ng=15)
        if filename is None:
            filename = cls.default_data_file()

        miepolymap.save(filename=filename)

    def __init__(self, diameter: np.ndarray = None,
                 wavelength: np.ndarray = None, ng: int = 15,
                 filename: str = None):
        '''
        Prepares maps of the first ng Legendre moments of the
        MiePolystyrene scattering phase function, over the specified
        range of spherical particle diameters and wavelengths.
        If ng >= 2, a map of gamma is prepared and if ng >= 3 maps of delta
        and sigma are prepared as well. The maps are used to obtain an
        initial estimate when calculating the MIE scattering phase function
        parameters from the given Legendre moments, gamma and/or delta.

        Parameters
        ----------
        diameter: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the micro-sphere diameter parameter of Mie phase function[m].
        wavelength: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the a wavelength parameter of Mie phase function[m].
        ng: int
            Maps are created for the first ng Legendre moments. If ng >= 2,
            a map of gamma is prepared and if g >= 3 maps of delta and sigma
            are prepared as well.
        filename: str
            File with saved data. The values of all the other parameters are
            ignored and restored from the file.
        ncostheta: int
            Number of nodes used to compute the Legendre moments. Use a
            large number (> 1000) for accurate results. If None, adaptive
            step integration is used, which is accurate but can become slow.

        Note
        ----
        The value of parameter ng should be >> 3 to accurately
        estimate the value of parameter sigma.

        Examples
        --------
        Prepares maps of gamma and delta, and estimates the MIE parameters
        given a) g and gamma are known b) gamma and delta are known.

        >>> import numpy as np
        >>>
        >>> m = MiePolystyreneMap(np.linspace(0.2e-6, 5e-6, 100), np.linspace(0.4e-6, 1e-6, 100), ng=3)
        >>> diameter, wavelength = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'diameter:', diameter, 'wavelength:', wavelength)
        >>> diameter, wavelength = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'diameter:', diameter, 'wavelength:', wavelength)
        >>>

        Load maps from the default file included in the data/pf folder.

        >>> m = MiePolystyreneMap.fromfile()
        >>> diameter, wavelength = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'diameter:', diameter, 'wavelength:', wavelength)
        >>> diameter, wavelength = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'diameter:', diameter, 'wavelength:', wavelength)
        >>>
        '''
        if diameter is None:
            diameter = np.linspace(0.2e-6, 5e-6, 100)

        if wavelength is None:
            wavelength = np.linspace(0.4e-6, 1e-6, 100)

        super().__init__(param1=diameter, param2=wavelength, ng=ng,
                         pf=MiePolystyrene, filename=filename)

    def diameter(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of diameters.
        '''
        return self.param1()

    def wavelength(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of wavelengths.
        '''
        return self.param2()

    def diameter_grid(self):
        '''
        Returns a 2D map (meshgrid) of the first parameter of the scattering
        phase function (diameter).
        '''
        return self.grid1()

    def wavelength_grid(self):
        '''
        Returns a 2D map (meshgrid) of the second parameter of the scattering
        phase function (wavelength).
        '''
        return self.grid2()
