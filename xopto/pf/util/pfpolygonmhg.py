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
import pickle
import sys

import numpy as np

from .helpers import g2gamma, g2delta
from .pfpolygonbase import GammaDeltaPolygonBase
from xopto.pf import MHg
from xopto import DATA_PATH


class MHgPolygon(GammaDeltaPolygonBase):
    DEFAULT_POLYGON_FILE = 'mhg_polygon.npz'
    NUM_POINTS = 100000

    def __init__(self, g_lim: Tuple[float, float] = (0, 0.99),
                 beta_lim: Tuple[float, float] = (0, 1), n: int = 5000,
                 verbose: bool = False, filename: str = None):
        '''
        Creates an object for fast polygon-based validation of the MHG
        scattering phase function domain.

        Parameters
        ----------
        g_lim: (float, float)
            Range of the g parameter of the MHG scattering phase function.
        beta_lim: (float, float)
            Range of the beta parameter of the MHG scattering phase function.
        n: int
            Number of points along each border of the
            scattering phase function domain.
        verbose: bool
            Print progress information to stdout.
        filename: str
            Load polygon data from a file. See the static method
            :py:meth:`~MHgPolygon.fromfile`.
        '''
        data = {'gamma': None, 'delta': None, 'filename': filename}

        if filename is None:
            data = self._calculate_polygon_pts(g_lim, beta_lim, n, verbose)

        super().__init__(**data)

    def _calculate_polygon_pts(
            self, g_lim: Tuple[float, float], beta_lim: Tuple[float, float], n: int,
            verbose: bool) -> dict:
        def _run_pf(g, beta):
            g = np.asarray(g).flatten()
            beta = np.asarray(beta).flatten()

            gamma = np.zeros(g.shape)
            delta = np.zeros(g.shape)
            for index in range(beta.size):
                gs = MHg(g[index], beta[index]).gs(3)
                gamma[index] = g2gamma(gs[1], gs[2])
                delta[index] = g2delta(gs[1], gs[3])
                if verbose:
                    sys.stdout.write(
                        '\rPreparing MHg gamma-delta polygon segment ... ' \
                        '{:.1f}% done'.format(100.0*(index + 1)/g.size))
                    sys.stdout.flush()

            return gamma, delta

        # bottom left polygon segment moving from left to right (increasing gamma)
        beta1 = np.linspace(beta_lim[0], beta_lim[1], n)
        g1 = np.tile(g_lim[0], (beta1.size,))
        gamma1mhg, delta1mhg = _run_pf(g1, beta1)

        # top polygon segment moving from left to right (increasing gamma)
        beta2 = np.linspace(beta_lim[0], beta_lim[1], n)
        g2 = np.tile(g_lim[1], (beta2.size,))
        gamma2mhg, delta2mhg = _run_pf(g2, beta2)

        # bottom right polygon segment moving from left to right (decreasing gamma)
        g3 = np.linspace(g_lim[0], g_lim[1], n)
        beta3 = np.tile(beta_lim[1], (g3.size,))
        gamma3mhg, delta3mhg = _run_pf(g3, beta3)


        # merge the three polygon segments in counter clockwise direction
        # starting with the bottom left segment followed by the bottom right and
        # top segments, kick out the duplicate junction points, except
        # the first and last points
        gamma = np.hstack(
            (
                gamma1mhg,
                gamma3mhg[1:],
                np.flip(gamma2mhg, axis=0)[1:],
            )
        )
        delta = np.hstack(
            (
                delta1mhg,
                delta3mhg[1:],
                np.flip(delta2mhg, axis=0)[1:],
            )
        )

        return {'gamma':gamma, 'delta':delta}

    def show(self, marker: str = '-', step: int = 1):
        '''
        Show the polygon using matplotlib.

        Parameters
        ----------
        marker: str
            Standard matplotlib line style, color and and point marker.
        step: int
            Sample the boundary points with the given step.
        '''
        MHgPolygon.show(self, marker, step, 'MHg')
