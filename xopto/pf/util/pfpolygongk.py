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

from xopto.pf import Gk
from .helpers import g2gamma, g2delta
from .pfpolygonbase import GammaDeltaPolygonBase
from xopto import DATA_PATH


class GkPolygon(GammaDeltaPolygonBase):
    DEFAULT_POLYGON_FILE = 'gk_polygon.npz'
    NUM_POINTS = 10000

    def __init__(self, ggk_lim: Tuple[float, float] = (0, 0.97),
                 a_lim: Tuple[float, float] = (-0.5, 10), n: int = 5000,
                 verbose: bool = False, filename: str = None):
        '''
        Creates an object for fast polygon-based validation of the Gk
        scattering phase function domain.

        Parameters
        ----------
        ggk_lim: tuple
            Range of the g_gk parameter of the scattering phase function.
        a_lim: tuple
            Range of the a_gk parameter of the scattering phase function.
        n: int
            Number of points along each border of the domain.
        verbose: bool
            Print progress information to stdout.
        filename: str
            Load polygon data from a file. See the static method
            :py:meth:`~GkPolygon.fromfile`.
        '''
        data = {'gamma': None, 'delta': None, 'filename': filename}

        if filename is None:
            data = self._calculate_polygon_pts(ggk_lim, a_lim, n, verbose)

        super().__init__(**data)

    def _calculate_polygon_pts(
            self, ggk_lim: Tuple[float, float], a_lim: Tuple[float, float], n: int,
            verbose: bool) -> dict:
        def _run_pf(g_gk, a_gk):
            g_gk = np.asarray(g_gk).flatten()
            a_gk = np.asarray(a_gk).flatten()

            gamma = np.zeros(a_gk.shape)
            delta = np.zeros(a_gk.shape)
            for index in range(a_gk.size):
                gs = Gk(g_gk[index], a_gk[index]).gs(3)
                gamma[index] = g2gamma(gs[1], gs[2])
                delta[index] = g2delta(gs[1], gs[3])
                if verbose:
                    sys.stdout.write(
                        '\rPreparing GK gamma-delta polygon segment ... ' \
                        '{:.1f}% done'.format(100.0*(index + 1)/a_gk.size))
                    sys.stdout.flush()

            return gamma, delta

        # top left polygon segment moving from left to right (increasing gamma)
        g1 = np.linspace(ggk_lim[0], ggk_lim[1], n)
        a1 = np.tile(a_lim[0], (g1.size,))
        gamma1gk, delta1gk = _run_pf(g1, a1)

        # top right polygon segment moving from left to right (increasing gamma)
        a2 = np.linspace(a_lim[0], a_lim[1], n)
        g2 = np.tile(ggk_lim[1], (a2.size,))
        gamma2gk, delta2gk = _run_pf(g2, a2)

        # bottom polygon segment moving from left to right (decreasing gamma)
        g3 = np.linspace(ggk_lim[0], ggk_lim[1], n)
        a3 = np.tile(a_lim[1], (g3.size,))
        gamma3gk, delta3gk = _run_pf(g3, a3)

        # merge the three polygon segments in counter clockwise direction
        # starting with the bottom segment followed by the top right and
        # top left segments, kick out the duplicate junction points except
        # the first and last points
        gamma = np.hstack(
            (
                gamma3gk,
                np.flip(gamma2gk, axis=0)[1:],
                np.flip(gamma1gk, axis=0)[1:],
            )
        )
        delta = np.hstack(
            (
                delta3gk,
                np.flip(delta2gk, axis=0)[1:],
                np.flip(delta1gk, axis=0)[1:],
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
        GammaDeltaPolygonBase.show(self, marker, step, 'Gk')
