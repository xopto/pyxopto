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
from ..mpc import MPc
from xopto import DATA_PATH


class MPcMap(PfMap2DBase):
    DEFAULT_MAP_FILE = 'mpc_map.npz'
    XLABEL = 'n'
    YLABEL = 'b'
    PLOTSCALEFACTORX = 1.0
    PLOTSCALEFACTORY = 1.0

    @classmethod
    def precalculate(cls, n: int = 500, filename: str = None,
                     verbose: bool = True):
        '''
        Precalculate MPc scattering phase function lookup table.

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
            print('\nCreating MPc map:')

        mpcmap = MPcMap(n=np.linspace(0.0, 50, n),
                        b=np.linspace(0.0, 1.0, n),
                        ng=15)

        if filename is None:
            filename = cls.default_data_file()

        mpcmap.save(filename=filename)

    def __init__(self, n=None, b=None, ng=15, filename=None,
                 ncostheta=None):
        '''
        Prepares maps of the first ng Legendre moments of the
        Modified Power of Cosine (MPC) scattering phase function, over the
        specified range of the MPC parameters n and b. If ng >= 2, a map of
        gamma is prepared and if ng >= 3 maps of delta and sigma are
        prepared as well. The maps are used to obtain an initial estimate
        when calculating the MPC scattering phase function parameters from given
        Legendre moments, gamma and/or delta.

        Parameters
        ----------
        n: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the n MPC parameter.
        b: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            the b MPC parameter (Power of cosine scattering phase function
            contribution).
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
        Prepares maps of gamma and delta, and estimates the MHG parameters
        given a) g and gamma are known b) gamma and delta are known.

        >>> import numpy as np
        >>>
        >>> m = MPcMap(np.linspace(0.0, 10, 100), np.linspace(0.0, 0.99, 100), ng=3)
        >>> n, b = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'n:', n, 'b:', b)
        >>> n, b = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'n:', n, 'b:', b)
        >>>

        Load maps from the default file included in the data/pf folder.

        >>> m = MPcMap.fromfile()
        >>> n, b = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'n:', n, 'b:', b)
        >>> n, b = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'n:', n, 'b:', b)
        >>>
        '''
        if n is None:
            n = np.linspace(0.0, 10.0, 100)

        if b is None:
            b = np.linspace(0.0, 1.0, 100)

        super().__init__(param1=n, param2=b, ng=ng,
                         pf=MPc, filename=filename, ncostheta=ncostheta)

    def n(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of the first MPC
        parameter n.
        '''
        return self.param1()

    def b(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of the first MPC
        parameter b.
        '''
        return self.param2()

    def n_grid(self) -> np.ndarray:
        '''
        Returns a 2D map (meshgrid) of the first parameter n of the MPC
        scattering phase function.
        '''
        return self.grid1()

    def b_grid(self) -> np.ndarray:
        '''
        Returns a 2D map (meshgrid) of the second parameter b of the MPC
        scattering phase function.
        '''
        return self.grid2()
