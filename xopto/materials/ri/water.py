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
from scipy.interpolate import interp1d

from .base import RefractiveIndex


class Daimon(RefractiveIndex):
    material = 'water'
    def __init__(self, t: float = 293.15):
        '''
        Refractive index of water as given in:
            Daimon and Masumura, Appl. Opt. 46, 3811-3820 (2007)

        Parameters
        ----------
        t: float
            Default temperature (K) of the medium.
        '''
        super().__init__(t=t, trange=(292, 297), wrange=(180e-9, 1129e-9))

    def __call__(self, wavelength: float or np.ndarray,
                 temperature: float or None = None) -> float or np.ndarray:
        '''
        Calculate the refractive index of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)
        temperature: float or None
            Use this temperature instead of the default temperature.

        Returns
        -------
        ri: float or np.ndarray
            Refractive index at the given wavelength(s).
        '''
        w = np.asarray(wavelength, dtype=np.float64)
        self.check_wavelengths(w)

        if temperature is None:
            temperature = self.t
        t = float(temperature)
        self.check_temperature(t)

        # wavelength must be in units of um
        inv_w2 = 1.0/(w*1e6)**2

        n_19 = np.sqrt(1.0 + 5.672526103e-1/(1.0 - 5.085550461e-3*inv_w2) +
                             1.736581125e-1/(1.0 - 1.814938654e-2*inv_w2) +
                             2.121531502e-2/(1.0 - 2.617260739e-2*inv_w2) +
                             1.138493213e-1/(1.0 - 1.073888649e1*inv_w2))

        n_20 = np.sqrt(1.0 + 5.684027565e-1/(1.0 - 5.101829712e-3*inv_w2) +
                             1.726177391e-1/(1.0 - 1.821153936e-2*inv_w2) +
                             2.086189578e-2/(1.0 - 2.620722293e-2*inv_w2) +
                             1.130748688e-1/(1.0 - 1.069792721e1*inv_w2))

        n_21p5 = np.sqrt(1.0 + 5.689093832e-1/(1.0 - 5.110301794e-3*inv_w2) +
                            1.719708856e-1/(1.0 - 1.825180155e-2*inv_w2) +
                            2.062501582e-2/(1.0 - 2.624158904e-2*inv_w2) +
                            1.123965424e-1/(1.0 - 1.067505178e1*inv_w2))

        n_24 = np.sqrt(1.0 + 5.666959820e-1/(1.0 - 5.084151894e-3*inv_w2) +
                             1.731900098e-1/(1.0 - 1.818488474e-2*inv_w2) +
                             2.095951857e-2/(1.0 - 2.625439472e-2*inv_w2) +
                             1.125228406e-1/(1.0 - 1.073842352e1*inv_w2))

        result = interp1d(
            [19.0, 20.0, 21.5, 24.0],
            np.vstack([n_19, n_20, n_21p5, n_24]).T,
            kind='linear',
            assume_sorted=True,
            axis=-1,
            fill_value='extrapolate'
        )(t - 273.15)

        if isinstance(wavelength, float):
            result = float(result)

        return result


class Schiebener(RefractiveIndex):
    material = 'water'
    def __init__(self, t: float = 293.15):
        '''
        Refractive index of water as given in:
            Schiebener and Straub, J.Phys. Chem. Ref. Data, (1990).

        Parameters
        ----------
        t: float
            Default temperature (K) of the medium.
        '''
        super().__init__(
            t=t, trange=(273, 498), wrange=(200e-9, 2500e-9))

    def __call__(self, wavelength: float or np.ndarray,
                 temperature: float or None = None) -> float or np.ndarray:
        '''
        Calculate the refractive index of the medium.

        Parameters
        ----------
        wavelength: float or np.ndarray
            Wavelength of light (m)
        temperature: float or None
            Use this temperature instead of the default temperature.

        Returns
        -------
        ri: float or np.ndarray
            Refractive index at the given wavelength(s).
        '''
        w = np.asarray(wavelength, dtype=np.float64)
        self.check_wavelengths(w)

        if temperature is None:
            temperature = self.t
        self.check_temperature(temperature)
        t = float(temperature)

        # wavelength must be in (nm)
        w = w*1e9

        rho_w = 1000 # density of water (kg/m^3)

        lambda_UV = 0.2292020
        lambda_IR = 5.432937

        a_0 = 0.243905091
        a_1 = 9.53518094e-3
        a_2 = -3.64358110e-3
        a_3 = 2.65666426e-4
        a_4 = 1.59189325e-3
        a_5 = 2.45733798e-3
        a_6 = 0.897478251
        a_7 = -1.63066183e-2

        lorentz = a_0 + a_1*rho_w/1e3 + a_2*t/273.15 + \
            a_3*(w/589.)**2*t/273.15 + \
            a_4*(w/589.)**-2 + a_5*((w/589.)**2-lambda_UV**2)**(-1) + \
            a_6*((w/589.)**2 - lambda_IR**2)**(-1) + a_7*(rho_w/1e3)**2

        return np.sqrt((2*lorentz+1)/(1-lorentz))

# Global instances for immediate use.
daimon = Daimon()
schiebener = Schiebener()

default = daimon
