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
from ..mhg import MHg
from .helpers import mhg_inverse_g1gamma, mhg_inverse_gammadelta
from xopto import DATA_PATH, USER_DATA_PATH


class MHgMap(PfMap2DBase):
    DEFAULT_MAP_FILE = 'mhg_map.npz'
    XLABEL = '$g_{mhg}$'
    YLABEL = '$b_{mgk}$'
    PLOTSCALEFACTORX = 1.0
    PLOTSCALEFACTORY = 1.0

    @classmethod
    def precalculate(cls, n: int = 100, filename: str = None,
                     verbose: bool = True):
        '''
        Precalculate MHg scattering phase function lookup table.

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
            print('\nCreating MHg map:')

        mhgmap = MHgMap(gmhg=np.linspace(0.0, 0.99, n),
                        b=np.linspace(0.0, 0.99, n),
                        ng=15)

        if filename is None:
            filename = cls.default_data_file()

        mhgmap.save(filename=filename)

    def __init__(self, gmhg: np.ndarray = None, b: np.ndarray = None,
                 ng: int = 15, filename: str = None, ncostheta: int = None):
        '''
        Prepares maps of the first ng Legendre moments of the
        Modified Henyey-Greenstein (MHG) scattering phase function, over the
        specified range of the MHG parameters g and b. If ng >= 2, a map of
        gamma is preapred and if ng >= 3 maps of delta and sigma are
        prepared as well. The maps are used to obtain an initial estimate
        when calculating the MHG scattering phase function parameters from
        Legendre moments, gamma and/or delta.

        Parameters
        ----------
        gmhg: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            parameter gg of the MHG scattering phase function.
        b: np.ndarray vector
            A Vector of equally spaced values defining the grid density of
            parameter b of the MHG scattering phase function (HG contribution).

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
        >>> m = MHgMap(np.linspace(0.0, 0.99, 100), np.linspace(0.0, 0.99, 100), ng=3)
        >>> gmhg, b = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'gmhg:', gmhg, 'b:', b)
        >>> gmhg, b = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'gmhg:', gmhg, 'b:', b)
        >>>

        Load maps from the default file included in the data/pf folder.

        >>> m = MHgMap.fromfile()
        >>> gmhg, b = m.invgammadelta(gamma=2.2, delta=3.5)
        >>> print('gamma=2.2, delta=3.5 ==>', 'gmhg:', gmhg, 'b:', b)
        >>> gmhg, b = m.invgamma(g=0.8, gamma=2.2)
        >>> print('g=0.8, gamma=2.2 ==>', 'gmhg:', gmhg, 'b:', b)
        >>>
        '''
        if gmhg is None:
            gmhg = np.linspace(0.0, 0.99, 100)

        if b is None:
            b = np.linspace(0.0, 0.99, 100)

        super().__init__(param1=gmhg, param2=b, ng=ng,
                         pf=MHg, filename=filename, ncostheta=ncostheta)

    def gmhg(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of g parameter of the
        MHg scattering phase function
        '''
        return self.param1()

    def b(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of b parameter of the
        MHg scattering phase function.
        '''
        return self.param2()

    def gmhg_grid(self):
        '''
        Returns a 2D map (meshgrid) of the first scattering phase function
        parameter (gmhg).
        '''
        return self.grid1()

    def b_grid(self):
        '''
        Returns a 2D map (meshgrid) of the second scattering phase function
        parameter (b).
        '''
        return self.grid2()

    def invgamma(self, g: float, gamma: float) -> float:
        '''
        Overloading the base class method with analytical solution to the
        inverse problem.

        Parameters
        ----------
        g: float
            Target value of the first Legendre moment.
        gamma: float
            Target value of parameter gamma.

        Returns
        -------
        param1, param2: float, float
            Phase function input parameters.
        '''
        mhg_g, mhg_b, valid = mhg_inverse_g1gamma(g, gamma)

        return float(mhg_g), float(mhg_b)

    def invgammadelta(self, gamma, delta, **kwargs):
        '''
        Overloading the base class method with analytical solution to the
        inverse problem.

        Parameters
        ----------
        gamma: float
            Target value of parameter gamma.
        delta: float
            Target value of parameter gamma.
        kwargs: dict
            Keyword arguments passed to the fmin_l_bfgs_b optimization
            function.

        Returns
        -------
        param1, param2: float, float
            Phase function input parameters.
        '''
        mhg_g, mhg_b, valid = mhg_inverse_gammadelta(gamma, delta)

        return float(mhg_g), float(mhg_b)
