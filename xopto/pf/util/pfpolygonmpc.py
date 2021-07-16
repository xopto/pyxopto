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
from xopto.pf import MPc
from xopto import DATA_PATH


class MPcPolygon(GammaDeltaPolygonBase):
    DEFAULT_POLYGON_FILE = 'mpc_polygon.npz'
    NUM_POINTS = 100000

    def __init__(self, n_lim=(0.0, 100.0), beta_lim=(0.0, 1.0), n=5000,
                 verbose=False, filename=None):
        '''
        Creates an object for fast polygon-based validation of the MPc
        scattering phase function domain.

        Parameters
        ----------
        n_lim: tuple
            Range of parameter n of the MPC scattering phase function.
        beta_lim: tuple
            Range of parameter beta of the MPC scattering phase function.
        n: int
            Number of points along each border of the domain.
        verbose: bool
            Print progress information to stdout.
        filename: str
            Load polygon data from a file. See the static method
            :py:meth:`~MPcPolygon.fromfile`.
        '''
        data = {'gamma': None, 'delta': None, 'filename': filename}

        if filename is None:
            data = self._calculate_polygon_pts(n_lim, beta_lim, n, verbose)

        super().__init__(**data)

    def _calculate_polygon_pts(
            self, n_lim: Tuple[float, float], beta_lim: Tuple[float, float], n: int,
            verbose: bool) -> dict:
        def _run_pf(n, beta):
            n = np.asarray(n).flatten()
            beta = np.asarray(beta).flatten()

            gamma = np.zeros(n.shape)
            delta = np.zeros(n.shape)
            for index in range(beta.size):
                gs = MPc(n[index], beta[index]).gs(3)
                gamma[index] = g2gamma(gs[1], gs[2])
                delta[index] = g2delta(gs[1], gs[3])
                if verbose:
                    sys.stdout.write(
                        '\rPreparing MPc gamma-delta polygon segment ... ' \
                        '{:.1f}% done'.format(100.0*(index + 1)/n.size))
                    sys.stdout.flush()

            return gamma, delta

        # bottom left polygon segment moving from left to right (increasing gamma)
        beta1 = np.linspace(beta_lim[0], beta_lim[1], n)
        g1 = np.tile(n_lim[0], (beta1.size,))
        gamma1mpc, delta1mpc = _run_pf(g1, beta1)

        # top segment polygon moving from left to right (increasing gamma)
        beta2 = np.linspace(beta_lim[0], beta_lim[1], n)
        g2 = np.tile(n_lim[1], (beta2.size,))
        gamma2mpc, delta2mpc = _run_pf(g2, beta2)

        # bottom right polygon segment moving from left to right (increasing gamma)
        g3 = np.linspace(n_lim[0], n_lim[1], n)
        beta3 = np.tile(beta_lim[1], (g3.size,))
        gamma3mpc, delta3mpc = _run_pf(g3, beta3)

        # merge the three polygon segments in counter clockwise direction
        # starting with the bottom left segment followed by the bottom right and
        # top segments, kick out the duplicate junction points, except
        # the first and last points
        gamma = np.hstack(
            (
                gamma1mpc,
                gamma3mpc[1:],
                np.flip(gamma2mpc, axis=0)[1:]
            )
        )
        delta = np.hstack(
            (
                delta1mpc,
                delta3mpc[1:],
                np.flip(delta2mpc, axis=0)[1:],
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
            Smaple the boundary points with the given step.
        '''
        GammaDeltaPolygonBase.show(self, marker, step, 'MPc')
