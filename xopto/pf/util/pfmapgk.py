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

import os.path

import numpy as np

from .pfmapbase import PfMap2DBase
from ..gk import Gk
from xopto import DATA_PATH


class GkMap(PfMap2DBase):
    DEFAULT_MAP_FILE = 'gk_map.npz'
    XLABEL = '$g_{gk}$'
    YLABEL = '$\\alpha_{gk}$'
    PLOTSCALEFACTORX = 1.0
    PLOTSCALEFACTORY = 1.0

    @classmethod
    def precalculate(cls, n: int = 100, filename: str = None,
                     verbose: bool = True):
        '''
        Precalculate the Gk scattering phase function lookup table.

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
            print('Creating Gk map:')

        gkmap = GkMap(ggk=np.linspace(0.0, 0.98, n),
                    a=np.linspace(-0.5, 10.0, n),
                    ng=15)

        if filename is None:
            filename = cls.default_data_file()

        gkmap.save(filename=filename)

    def __init__(self, ggk: np.ndarray = None, a: np.ndarray = None, ng:
                 int = 15, filename: str = None, ncostheta: int = None):
        '''
        Prepares maps of the first ng Legendre moments of the
        Gegenbauer Kernel (GK) scattering phase function, over the specified
        range of the GK parameters gg and a. If ng >= 2, a map of
        gamma is prepared and if ng >= 3 maps of delta and sigma are
        prepared as well. The maps are used to obtain an initial estimate
        when calculating the GK scattering phase function parameters
        from given Legendre moments, gamma and/or delta.

        Parameters
        ----------
        ggk: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the gg parameter of the GK scattering phase function.

        a: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the a parameter of the GK scattering phase function.

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
        Prepares maps of gamma and delta, and estimates the GK parameters
        given a) g and gamma are known b) gamma and delta are known.

        >>> import numpy as np
        >>>
        >>> m = GkMap(np.linspace(0.0, 0.98, 100), np.linspace(-0.5, 5.0, 100), ng=3)
        >>> ggk, a = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'ggk:', ggk, 'a:', a)
        >>> ggk, a = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'ggk:', ggk, 'a:', a)
        >>>

        Load maps from the default file included in the data/pf folder.

        >>> m = GkMap.fromfile()
        >>> ggk, a = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'ggk:', ggk, 'a:', a)
        >>> ggk, a = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'ggk:', ggk, 'a:', a)
        >>>
        '''
        if ggk is None:
            ggk = np.linspace(0.0, 0.98, 100)

        if a is None:
            a = np.linspace(-0.5, 10.0, 100)

        super().__init__(param1=ggk, param2=a, ng=ng,
                         pf=Gk, filename=filename, ncostheta=ncostheta)

    def ggk(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of gg GK parameter
        '''
        return self.param1()

    def a(self)  -> np.ndarray:
        '''
        Returns a vector of points defining the grid of parameter a of the GK
        scattering phase function.
        '''
        return self.param2()

    def ggk_grid(self) -> np.ndarray:
        '''
        Returns a 2D map (meshgrid) of the first parameter ggk of the GK
        scattering phase function input parameter.
        '''
        return self.grid1()

    def a_grid(self) -> np.ndarray:
        '''
        Returns a 2D map (meshgrid) of the second parameter a of the GK
        scattering phase function input parameter.
        '''
        return self.grid2()
