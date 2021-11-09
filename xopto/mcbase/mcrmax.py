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

import numpy as np
from scipy.integrate import simps

from xopto.diffusion.srr import SRDA


class RmaxConst:
    def __init__(self, rmax: float = np.inf):
        '''
        Simple constant simulation radius constructor.
        Parameters
        ----------
        rmax: float
            Maximum simulation radius
        '''
        self._rmax = float(rmax)

    def __call__(self, *args, **kwargs) -> float:
        '''
        Returns the simulation radius passed to the constructor
        '''
        return self._rmax

    def __str__(self):
        return 'RmaxConst(rmax={})'.format(self._rmax)

    def __repr__(self):
        return '{:s} # {}'.format(self.__str__(), id(self))


class RmaxDiffusion:
    STEPS = 100000
    R_MAX = 500e-3

    def __init__(self, nmedium: float, noutside: float,
                 acceptance: float,
                 rmin: float = 10e-3, dr: float = 5e-3,
                 rtol: float = 1e-3):
        '''
        Initiates the diffusion approximation solution for spatially resolved
        reflectance, which is then used to calculate minimal acceptable radius
        rmin beyond which the acquired residual reflectance represents less than
        rtol of the total acquired reflectance.

        Parameters
        ----------
        nmedium: float
            Refractive index of the medium.
        noutside: float
            Refractive index above the medium.
        acceptance: float
            Acceptance angle of the detector (radians).
        rmin: float
            Lowest value of the minimal acceptable radius (m) (includes dr)
            to be used in the Monte Carlo simulations for different optical
            properties (i.e., for all of the cases the simulation radius will
            be equal or greater than rmin).
        dr: float
            Iteration step for finding minimal acceptable radius (m).
        rtol: float
            Maximum acceptable amount of residual to total reflectance
            calculated at the minimal acceptable radius.
        '''

        self._rmin = float(rmin)
        self._dr = float(dr)
        self._rtol = float(rtol)
        self._nmedium = float(nmedium)
        self._noutside = float(noutside)
        self._acceptance = float(acceptance)

        self._da = SRDA(self._nmedium, self._noutside, self._acceptance)

    def __call__(self, mua: float, musr: float) -> float:
        '''
        Calculates the minimal acceptable radius beyond which the acquired
        residual reflectance represents less than rtol of the total acquired
        reflectance.

        Parameters
        ----------
        mua: float
            Absorption coefficient of the medium in 1/m.
        musr: float
            Reduced scattering coefficient of the medium in 1/m.

        Returns
        -------
        rmax: float
            Minimum simulation radius to be used in Monte Carlo simulations
            (already increased by dr).
        '''

        self._mua = mua
        self._musr = musr

        return self._find_rmax()

    def _r_total(self) -> float:
        r = np.linspace(0, self.R_MAX, self.STEPS)
        return simps(
            2.0*np.pi*self._da(r, self._mua, self._musr)*r, dx=r[1] - r[0])

    def _residual(self, i) -> float:
        r = np.linspace(self._rmin + (i - 1)*self._dr, self.R_MAX, self.STEPS)
        return simps(
            2.0*np.pi*self._da(r, self._mua, self._musr)*r, dx=r[1] - r[0])

    def _find_rmax(self) -> float:
        i = 0
        r_total = self._r_total()
        residual = self._residual(i)/r_total

        while residual > self._rtol/2.0:
            # print(i)
            i += 1
            residual = self._residual(i)/r_total

        return self._rmin + i * self._dr

    def __str__(self):
        return 'RmaxDiffusion(nmedium={:f}, noutside={:f}, acceptance={:f}, ' \
               'rmin={:f}, dr={:f}, rtol={:f})'.format(
                   self._nmedium, self._noutside, self._acceptance,
                   self._rmin, self._dr, self._rtol)

    def __repr__(self):
        return '{:s} # {}'.format(self.__str__(), id(self))


if __name__ == "__main__":
    import time
    from matplotlib import pyplot as pp

    mua_vector = np.linspace(0.01e2, 5.0e2, 30)
    musr_vector = np.linspace(5.0e2, 35.0e2, 50)
    mua_array, musr_array = np.meshgrid(mua_vector, musr_vector, indexing='ij')
    mua_array_flat = mua_array.flatten()
    musr_array_flat =musr_array.flatten()
    rmax_array_flat = np.empty_like(mua_array_flat)

    rmin = 10.0e-3
    dr = 5.0e-3
    rtol = 1.0e-3

    rmax_estimator = RmaxDiffusion(
        1.33, 1.0, np.deg2rad(10.0), rmin=10e-3, dr=5e-3, rtol=1e-3)

    t1 = time.perf_counter()
    for index, (mua, musr) in enumerate(zip(mua_array_flat, musr_array_flat)):
        rmax_array_flat[index] = rmax_estimator(mua, musr)
    t2 = time.perf_counter()

    dt_sample = (t2 - t1)*1e3/mua_array.size
    print('Estimation time per sample: {:.3f} ms'.format(dt_sample))

    rmax_array = np.reshape(rmax_array_flat, (mua_vector.size, musr_vector.size))

    pp.imshow(np.flipud(1e3 * rmax_array),
              extent=[musr_vector[0]*1e-2, musr_vector[-1]*1e-2,
                      mua_vector[0]*1e-2, mua_vector[-1]*1e-2],
              aspect='auto', cmap='hot_r')
    pp.xlabel('Musr (1/cm)')
    pp.ylabel('Mua (1/cm)')
    cbar = pp.colorbar()
    cbar.ax.set_ylabel('Radius (mm)')
    pp.title('Minimum simulation radius for MC')
    pp.show()
