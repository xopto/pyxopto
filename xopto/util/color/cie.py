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
from typing import Tuple

import numpy as np
from scipy.integrate import simps
from scipy.interpolate import interp1d

from xopto import DATA_PATH


STANDARD_WAVELENGTH_STEP = 5e-9
'''
Standard wavelength step/resolution (m) used in the datasets.
'''

STANDARD_WAVELENGTH_RANGE = (380e-9, 740e-9)
'''
Standard wavelength range (m) used for computations.
'''

STANDARD_WAVELENGTHS = np.arange(
    STANDARD_WAVELENGTH_RANGE[0],
    STANDARD_WAVELENGTH_RANGE[1] + STANDARD_WAVELENGTH_STEP*0.5,
    STANDARD_WAVELENGTH_STEP)
'''
Vector of standard wavelengths (m) that can be used by external libraries
to better match the internal resolution/grid of datasets.
'''

STANDARD_WAVELENGTHS_5_NM = STANDARD_WAVELENGTHS
'''
Vector of standard wavelengths (m) with resolution of 5 nm that can be used
by external libraries to better match the internal resolution/grid of datasets.
'''

STANDARD_WAVELENGTHS_10_NM = STANDARD_WAVELENGTHS[::2]
'''
Vector of standard wavelengths (m) with resolution of 10 nm that can be used
by external libraries to better match the internal resolution/grid of datasets.
'''

STANDARD_WAVELENGTHS_20_NM = STANDARD_WAVELENGTHS[::4]
'''
Vector of standard wavelengths (m) with resolution of 20 nm that can be used
by external libraries to better match the internal resolution/grid of datasets.
'''


RGB_COLOR_SPACES = {
    'adobe':         {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd65', 'r': (0.6400, 0.3300, 0.297361), 'g': (0.2100, 0.7100, 0.627355), 'b': (0.1500, 0.0600, 0.075285)},
    'apple':         {'companding': 'gamma', 'gamma': 1.8, 'illuminant': 'd65', 'r': (0.6250, 0.3400, 0.244634), 'g': (0.2800, 0.5950, 0.672034), 'b': (0.1550, 0.0700, 0.083332)},
    'best':          {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd50', 'r': (0.7347 ,0.2653, 0.228457), 'g': (0.2150, 0.7750, 0.737352), 'b': (0.1300, 0.0350, 0.034191)},
    'beta':          {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd50', 'r': (0.6888, 0.3112, 0.303273), 'g': (0.1986, 0.7551, 0.663786), 'b': (0.1265, 0.0352, 0.032941)},
    'bruce':         {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd65', 'r': (0.6400, 0.3300, 0.240995), 'g': (0.2800, 0.6500, 0.683554), 'b': (0.1500, 0.0600, 0.075452)},
    'cie':           {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'e',   'r': (0.7350, 0.2650, 0.176204), 'g': (0.2740, 0.7170, 0.812985), 'b': (0.1670, 0.0090, 0.010811)},
    'color_match':   {'companding': 'gamma', 'gamma': 1.8, 'illuminant': 'd50', 'r': (0.6300, 0.3400, 0.274884), 'g': (0.2950, 0.6050, 0.658132), 'b': (0.1500, 0.0750, 0.066985)},
    'don':           {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd50', 'r': (0.6960, 0.3000, 0.278350), 'g': (0.2150, 0.7650, 0.687970), 'b': (0.1300, 0.0350, 0.033680)},
    'eci_v2':        {'companding': 'l*',    'gamma': 'l*','illuminant': 'd50',	'r': (0.6700, 0.3300, 0.320250), 'g': (0.2100, 0.7100, 0.602071), 'b': (0.1400, 0.0800, 0.077679)},
    'ekta_space_ps5':{'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd50', 'r': (0.6950, 0.3050, 0.260629), 'g': (0.2600, 0.7000, 0.734946), 'b': (0.1100, 0.0050, 0.004425)},
    'ntsc':          {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'c',   'r': (0.6700, 0.3300, 0.298839), 'g': (0.2100, 0.7100, 0.586811), 'b': (0.1400, 0.0800, 0.114350)},
    'pal_secam':     {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd65', 'r': (0.6400, 0.3300, 0.222021), 'g': (0.2900, 0.6000, 0.706645), 'b': (0.1500, 0.0600, 0.071334)},
    'prophoto':      {'companding': 'gamma', 'gamma': 1.8, 'illuminant': 'd50', 'r': (0.7347, 0.2653, 0.288040), 'g': (0.1596, 0.8404, 0.711874), 'b': (0.0366, 0.0001, 0.000086)},
    'smpte-c':       {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd65', 'r': (0.6300, 0.3400, 0.212395), 'g': (0.3100, 0.5950, 0.701049), 'b': (0.1550, 0.0700, 0.086556)},
    'srgb':          {'companding': 'srgb',  'gamma': 2.2, 'illuminant': 'd65', 'r': (0.6400, 0.3300, 0.212656), 'g': (0.3000, 0.6000, 0.715158), 'b': (0.1500, 0.0600, 0.072186)},
    'wide_gamut':    {'companding': 'gamma', 'gamma': 2.2, 'illuminant': 'd50', 'r': (0.7350, 0.2650, 0.258187), 'g': (0.1150, 0.8260, 0.724938), 'b': (0.1570, 0.0180, 0.016875)},
}
'''
Specifications of some standard RGB color spaces.
'''

class _Data:
    ''' Internal class used as estimator/interpolator for all the datasets. '''
    class Estimator:
        def __init__(self, x, y, *args, **kwargs):
            self._x = x
            self._y = y
            self._interp_fun = interp1d(x, y, *args, **kwargs)
            self._range = (np.min(x), np.max(x))

        def __call__(self, x=None):
            if x is None:
                return self._x, self._y

            return self._interp_fun(x)

        def range(self):
            return self._range

        def checkRange(self, x):
            return np.all(x >= self._range[0]) and np.all(x <= self._range[1])

    def __init__(self):
        ROOT_DIR = os.path.join(DATA_PATH, 'color')

        self.d50_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_d50.npz'))['data']
        self.d55_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_d55.npz'))['data']
        self.d65_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_d65.npz'))['data']
        self.d75_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_d75.npz'))['data']
        self.a_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_a.npz'))['data']
        self.c_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_c.npz'))['data']
        self.e_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_e.npz'))['data']

        self.fl_1_12_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_fl_1_12.npz'))['data']
        self.fl3_1_15_source = \
            np.load(os.path.join(ROOT_DIR, 'cie_source_fl3_1_15.npz'))['data']

        self.cmf_1931_xyz = \
            np.load(os.path.join(ROOT_DIR, 'cie_cmf_1931_xyz.npz'))['data']

        self.cmf_1931_rgb = \
            np.load(os.path.join(ROOT_DIR, 'cie_cmf_1931_rgb.npz'))['data']

        self.cmf_1964_xyz = \
            np.load(os.path.join(ROOT_DIR, 'cie_cmf_1964_xyz.npz'))['data']

        self.cmf_1964_rgb = \
            np.load(os.path.join(ROOT_DIR, 'cie_cmf_1964_rgb.npz'))['data']

_DATA = _Data()

_d50_source_estimator = _Data.Estimator(
    _DATA.d50_source[:, 0], _DATA.d50_source[:, 1],
    bounds_error=False, fill_value=(_DATA.d50_source[0, 1], _DATA.d50_source[-1, 1])
)
_d55_source_estimator = _Data.Estimator(
    _DATA.d55_source[:, 0], _DATA.d55_source[:, 1],
    bounds_error=False, fill_value=(_DATA.d55_source[0, 1], _DATA.d55_source[-1, 1])
)
_d65_source_estimator = _Data.Estimator(
    _DATA.d65_source[:, 0], _DATA.d65_source[:, 1],
    bounds_error=False, fill_value=(_DATA.d65_source[0, 1], _DATA.d65_source[-1, 1])
)
_d75_source_estimator = _Data.Estimator(
    _DATA.d75_source[:, 0], _DATA.d75_source[:, 1],
    bounds_error=False, fill_value=(_DATA.d75_source[0, 1], _DATA.d75_source[-1, 1])
)
_a_source_estimator = _Data.Estimator(
    _DATA.a_source[:, 0], _DATA.a_source[:, 1],
    bounds_error=False, fill_value=(_DATA.a_source[0, 1], _DATA.a_source[-1, 1])
)
_c_source_estimator = _Data.Estimator(
    _DATA.c_source[:, 0], _DATA.c_source[:, 1],
    bounds_error=False, fill_value=(_DATA.c_source[0, 1], _DATA.c_source[-1, 1])
)
_e_source_estimator = _Data.Estimator(
    _DATA.e_source[:, 0], _DATA.e_source[:, 1],
    bounds_error=False, fill_value=(_DATA.e_source[0, 1], _DATA.e_source[-1, 1])
)
_cmf_1931_xyz_estimator = _Data.Estimator(
    _DATA.cmf_1931_xyz[:, 0], _DATA.cmf_1931_xyz[:, 1:].T,
    bounds_error=False,
    fill_value=(_DATA.cmf_1931_xyz[0, 1:], _DATA.cmf_1931_xyz[-1, 1:])
)
_cmf_1964_xyz_estimator = _Data.Estimator(
    _DATA.cmf_1964_xyz[:, 0], _DATA.cmf_1964_xyz[:, 1:].T,
    bounds_error=False,
    fill_value=(_DATA.cmf_1964_xyz[0, 1:], _DATA.cmf_1964_xyz[-1, 1:])
)
_cmf_1931_rgb_estimator = _Data.Estimator(
    _DATA.cmf_1931_rgb[:, 0], _DATA.cmf_1931_rgb[:, 1:].T,
    bounds_error=False,
    fill_value=(_DATA.cmf_1931_rgb[0, 1:], _DATA.cmf_1931_rgb[-1, 1:])
)
_cmf_1964_rgb_estimator = _Data.Estimator(
    _DATA.cmf_1964_rgb[:, 0], _DATA.cmf_1964_rgb[:, 1:].T,
    bounds_error=False,
    fill_value=(_DATA.cmf_1964_rgb[0, 1:], _DATA.cmf_1964_rgb[-1, 1:])
)

for i in range(12):
    estimator_name = '_fl_{}_source_estimator'.format(i + 1) 
    globals()[estimator_name] = _Data.Estimator(
        _DATA.fl_1_12_source[:, 0], _DATA.fl_1_12_source[:, i + 1],
        bounds_error=False, fill_value=(_DATA.fl_1_12_source[0, i + 1], _DATA.fl_1_12_source[-1, i + 1])
    )

for i in range(15):
    estimator_name = '_fl3_{}_source_estimator'.format(i + 1) 
    globals()[estimator_name] = _Data.Estimator(
        _DATA.fl3_1_15_source[:, 0], _DATA.fl3_1_15_source[:, i + 1],
        bounds_error=False, fill_value=(_DATA.fl3_1_15_source[0, i + 1], _DATA.fl3_1_15_source[-1, i + 1])
    )

CMF_RGB_2_XYZ = np.array(
    [
        [2.7688, 1.7517, 1.1301],
        [1.0000, 4.5906, 0.0601],
        [0.0000, 0.0565, 5.5942],
    ],
    dtype=np.float64
)

'''
Linear transformation from RGB-color matching functions to XYZ color-matching
functions.
'''

CMF_RGB_2_XYZ = np.linalg.inv(CMF_RGB_2_XYZ)
'''
Linear transformation from RGB color-matching functions to XYZ color-matching
functions.
'''

D50_1931_XYZ = (96.4448, 100.0, 82.6079)
'''
Standard CIE A source X, Y, Z color coordinates for CIE 1931 color-matching
functions.
'''

D50_1931_xy = (0.34562, 0.35836)
'''
Standard CIE A source X, Y, Z color coordinates for CIE 1931 color-matching
functions.
'''

D50_1964_XYZ = (96.7191, 100.0, 81.4233)
'''
Standard CIE A source X, Y, Z color coordinates for CIE 1931 color-matching
functions.
'''

D50_1964_xy = (0.34773, 0.35953)
'''
Standard CIE A source X, Y, Z color coordinates for CIE 1931 color-matching
functions.
'''

D65_1931_XYZ = (95.0489, 100.0, 108.884)
'''
Standard CIE D65 source X, Y, Z color coordinates for CIE 1931 color-matching
functions.
'''

D65_1931_xy = (0.31271, 0.32902)
'''
Standard CIE D65 source x, y chromaticities for CIE 1931 color-matching
functions.
'''

D65_1964_XYZ = (94.8110, 100.0, 107.304)
'''
Standard CIE A source X, Y, Z color coordinates for CIE 1931 color-matching
functions.
'''

D65_1964_xy = (0.31382, 0.33100)
'''
Standard CIE A source x, y chromaticities for CIE 1931 color-matching functions.
'''

class Observer:
    def __init__(self, estimator: _Data.Estimator,
                 rgb_estimator: _Data.Estimator,
                 fov_deg: float, name: str = 'Observer'):
        '''
        Parameters
        ----------
        estimator: _Data.Estimator
            Estimator of color matching functions for the X, Y and Z color
            components as a function of wavelength.
        rgb_estimator: _Data.Estimator
            Estimator of color matching functions for the linear CIE
            R, G and B color components as a function of wavelength.
        fov_deg: float
            Field of fow (degree).
        name: str
            Observer name.
        '''
        self._xyz_estimator = estimator
        self._rgb_estimator = rgb_estimator
        self._name = name
        self._fov_deg = float(fov_deg)
        self._fov_rad = float(np.deg2rad(fov_deg))

    def _get_fov_deg(self) -> float:
        return self._fov_deg
    fov_deg = property(_get_fov_deg, None, None,
                       'Observer field of view (deg)')

    def _get_fov_rad(self) -> float:
        return self._fov_rad
    fov = property(_get_fov_rad, None, None,
                  'Observer field of view (rad)')

    def _get_name(self) -> str:
        return self._name
    name = property(_get_name, None, None, 'Observer name.')

    def xyz_cmf(self, wavelengths: None or float or np.ndarray = None) \
            -> np.ndarray:
        '''
        Returns the values of the XYZ color matching functions at the specified
        wavelengths. If no wavelengths are specified, the values at
        standard wavelengths :py:attr:`STANDARD_WAVELENGTHS` are returned.

        Parameters
        ----------
        wavelengths: None or float or np.ndarray
            Wavelengths of light (m) at which to compute the color matching
            functions. If no wavelengths are specified,
            the values at standard wavelengths :py:attr:`STANDARD_WAVELENGTHS`
            are returned. 

        Returns
        -------
        XYZ: np.ndarray
            Color matching functions X, Y and Z as at the specified wavelengths
            as a numpy array of shape (wavelengths.size, 3). The values of
            color matching functions are in rows of the returned array,
            i.e. X in the first, Y in the second and Z in the third row.
        '''
        if wavelengths is None:
            wavelengths = STANDARD_WAVELENGTHS

        return self._xyz_estimator(wavelengths)

    def rgb_cmf(self, wavelengths: None or float or np.ndarray = None) \
            -> np.ndarray:
        '''
        Returns the values of the RGB color matching functions at the specified
        wavelengths. If no wavelengths are specified, the values at
        standard wavelengths :py:attr:`STANDARD_WAVELENGTHS` are returned.

        Parameters
        ----------
        wavelengths: None or float or np.ndarray
            Wavelengths of light (m) at which to compute the color matching
            functions. If no wavelengths are specified,
            the values at standard wavelengths :py:attr:`STANDARD_WAVELENGTHS`
            are returned. 

        Returns
        -------
        RGB: np.ndarray
            Color matching functions R, G and B as at the specified wavelengths
            as a numpy array of shape (wavelengths.size, 3). The values of
            color matching functions are in rows of the returned array,
            i.e. R in the first, G in the second and B in the third row.
        '''
        if wavelengths is None:
            wavelengths = STANDARD_WAVELENGTHS

        return self._rgb_estimator(wavelengths)

    def xyz_cmf_raw(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns raw XYZ color matching functions and corresponding wavelengths.
        The values of the color matching function for component X are in the
        first row, for Y in the second row and for Z in the third row.

        Returns
        -------
        wavelengths: np.ndarray
            Wavelengths at which the raw spectral power densities are defined.
        cmf: np.ndarray
            Raw color matching functions defined at the returned wavelengths.
            The shape of the array is (wavelengths.size, 3).
            The values of the color matching function for component X are
            in the first row, for Y in the second row and for Z
            in the third row.
        '''
        return self._xyz_estimator()

    def rgb_cmf_raw(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns raw RGB color matching functions and corresponding wavelengths.
        The values of the color matching function for component X are in the
        first row, for Y in the second row and for Z in the third row.

        Returns
        -------
        wavelengths: np.ndarray
            Wavelengths at which the raw spectral power densities are defined.
        cmf: np.ndarray
            Raw color matching functions defined at the returned wavelengths.
            The shape of the array is (wavelengths.size, 3).
            The values of the color matching function for component R are
            in the first row, for G in the second row and for B
            in the third row.
        '''
        return self._rgb_estimator()

    def XYZ(self, wavelengths: np.ndarray, spectrum: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Compute CIE X, Y, Z color components from the given spectrum.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths of light (m) in the spectrum.

        spectrum: np.ndarray
            Intensity at the given wavelengths

        Returns
        -------
        X, Y, Z: (np.ndarray, np.ndarray, np.ndarray)
            Values of the color components as a tuple of numpy arrays.
        '''
        cmf = self.xyz_cmf(wavelengths)

        return \
            simps(spectrum*cmf[0], wavelengths), \
            simps(spectrum*cmf[1], wavelengths), \
            simps(spectrum*cmf[2], wavelengths)

    def rgb(self, wavelengths: np.ndarray, spectrum: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Compute linear CIE R, G, B color components from the given spectrum.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths of light (m) in the spectrum.

        spectrum: np.ndarray
            Intensity at the given wavelengths

        Returns
        -------
        R, G, B: (np.ndarray, np.ndarray, np.ndarray)
            Values of the color components as a tuple of numpy arrays.
        '''
        rgb_cmf = self.rgb_cmf(wavelengths)

        return \
            simps(spectrum*rgb_cmf[0], wavelengths), \
            simps(spectrum*rgb_cmf[1], wavelengths), \
            simps(spectrum*rgb_cmf[2], wavelengths)

    def xyY(self, wavelengths: np.ndarray, spectrum: np.ndarray) -> np.ndarray:
        '''
        Compute CIE x, y chromaticities and Y luminosity from the given
        spectrum.

        Parameters
        ----------
        wavelength: np.ndarray
            Wavelengths of light (m) in the spectrum.

        spectrum: np.ndarray
            Intensity at the given wavelengths

        Returns
        -------
        x, y, Y: (np.ndarray, np.ndarray, np.ndarray)
            Values of the chromaticities x and y, and luminosity Y.
        '''
        X, Y, Z = self.XYZ(wavelengths, spectrum)
        sum_XYZ = X + Y + Z

        return X/sum_XYZ, Y/sum_XYZ, Y

    @staticmethod
    def _f_lab(t: np.ndarray) -> np.ndarray:
        delta = 6.0/29.0

        if hasattr(t, '__len__'):
            mask = t > delta**3
            inv_mask = np.logical_not(mask)

            result = np.zeros_like(t)

            result[mask] = np.cbrt(t[mask])
            result[inv_mask] = t[inv_mask]/(3.0*delta**2) + 4.0/29.0
        else:
            if t > delta**3:
                result = np.cbrt(t)
            else:
                result = t/(3.0*delta**2) + 4.0/29.0

        return result

    def XYZ2Lab(self, xyz: np.ndarray, illuminant: 'Illuminant') -> np.ndarray:
        '''
        Compute CIE L*, a*, b* color components from CIE X, Y, Z color
        components for the given illuminant.

        Parameters
        ----------
        xyz: np.ndarray
            CIE X, Y, Z color components that were computed with this observer
            or transformed to this observer. Multiple entries must be stored
            in the row-wise.  Do NOT use coordinates that were
            computed with a different observer!

        Returns
        -------
        lab: np.ndarray
            Values of the L*, a* and b* color components in rows.

        Note
        ----
        Illuminant is used to compute the required CIE XYZ tristimulus values
        of the reference white point.
        '''
        Xn, Yn, Zn = illuminant.XYZ(self)

        f_y_yn = self._f_lab(Y/Yn)
        L = 116.0*f_y_yn - 16.0
        a = 500.0*(self._f_lab(X/Xn) - f_y_yn)
        b = 200.0*(f_y_yn - self._f_lab(Z/Zn))

        return L, a, b

    def __str__(self):
        return '{} Observer'.format(self._name)

    def __repr__(self):
        return self.__str__()


CIE1931 = Observer(_cmf_1931_xyz_estimator, _cmf_1931_rgb_estimator,
                   fov_deg=2.0, name='CIE 1931')
''' Instance of CIE 1931 observer (2 degree field of view). '''

CIE1964 = Observer(_cmf_1964_xyz_estimator, _cmf_1964_rgb_estimator,
                   fov_deg=10.0, name='CIE 1964')
''' Instance of CIE 1964 observer (10 degree field of view). '''


def standard_observer(observer: str or Observer) -> Observer:
    '''
    Returns the requested observer.

    Parameters
    ----------
    observer: str or Observer
        Selected standard observer "1931" ("cie1931") or "1964" ("cie1964").

    Returns
    -------
    Observer: observer
        The requested observer.
    '''
    obj = {
        '1931': CIE1931, 'cie1931': CIE1931, CIE1931: CIE1931,
        '1964': CIE1964, 'cie1964': CIE1964, CIE1964: CIE1964,
    }.get(str(observer).lower())

    if observer is None:
        raise ValueError('Unknown observer "{}"!'.format(observer))

    return obj

class RGB:
    @staticmethod
    def gamma_companding(v: np.ndarray, gamma: float) -> np.ndarray:
        '''
        Nonlinear intensity transformation that should be applied to
        all the linear RGB coordinates computed from XYZ by a
        linear transformation. See the individual RGB color spaces for details
        on the value of gamma. Note that some color spaces such as the
        sRGB and ECI RGB require different companding.
        '''
        return v**(1.0/gamma)

    @staticmethod
    def srgb_companding(v: np.ndarray) -> np.ndarray:
        '''
        Nonlinear intensity transformation that should be applied to
        all the linear sRGB coordinates computed from XYZ by a
        linear transformation. Note that some implementations use gamma
        companding with :math:`\\gamma=2.2`.
        '''
        if hasattr(v, '__len__'):
            v = np.asarray(v, dtype=np.float64)
            out = np.zeros_like(v)
            mask = v <= 0.0031308
            out[mask] = 12.92*v[mask]
            inv_mask = np.logical_not(mask)
            out[inv_mask] = 1.055*v[inv_mask]**(1.0/2.4) - 0.055
        else:
            if v <= 0.0031308:
                out = 12.92*v
            else:
                out = 1.055*v**(1.0/2.4) - 0.055

        return out

    @staticmethod
    def l_companding(v: np.ndarray) -> np.ndarray:
        '''
        Nonlinear intensity transformation that should be applied to
        all the linear ECI RGB coordinates.
        '''
        e = 0.008856 # 216/24389
        k = 903.3    # 24389/27

        if hasattr(v, '__len__'):
            v = np.asarray(v, dtype=np.float64)
            out = np.zeros_like(v)
            mask = v <= e
            out[mask] = v[mask]*k*0.01
            inv_mask = np.logical_not(mask)
            out[inv_mask] = 1.16*np.cbrt(v[inv_mask]) - 0.16
        else:
            if v <= e:
                out = v*k*0.01
            else:
                out = 1.16*np.cbrt(v) - 0.16

        return out

    def __init__(self, name: str, data: dict, pretty_name: str ='RGB'):
        '''
        Initializes a RGB color space from data.

        Parameters
        ----------
        name: str
            Name of the RGB color space as given in :py:attr:`RGB_COLOR_SPACES`.
        data: dict
            Color space data under the specified name as a dict. See
            :py:attr:`RGB_COLOR_SPACES`.
        pretty_name: str
            Display name.
        '''
        self._data = data
        self._pretty_name = pretty_name

    def _get_illuminant(self) -> 'Illuminant':
        '''
        Returns:
        Illuminant:
            Illuminant of the RGB color space.
        '''
        return standard_illuminant(self._data['illuminant'])

    illuminant = property(_get_illuminant, None ,None,
                                   'Color space illuminant.')

    def to_xyz_transformation(self, observer: Observer,
                              normalize: bool or float = True):
        '''
        Compute transformation from normalized linear RGB coordinates
        to XYZ coordinates of the specified observer.

        Parameters
        ----------
        observer: Observer
            Target observer.
        normalize: bool of float
            If True, the transformation matrix is computed from the XYZ
            coordinates of the illuminant the are normalized
            by the illuminant luminosity Y as
            :math:`X_n, Y_n, Z_n = X/Y, 1.0, Z/Y`.
            If a floating-point value, the transformation matrix is computed
            from the XYZ coordinates of the illuminant the are weighted
            by the given value.

        Returns
        -------
        T: np.ndarray
            A 3x3 transformation matrix as a numpy array.
        '''
        xyzw = self.illuminant.XYZ(observer, normalize=normalize)

        xr, yr = self._data['r'][:2]
        xg, yg = self._data['g'][:2]
        xb, yb = self._data['b'][:2]

        Xr, Yr, Zr = xr/yr, 1.0, (1.0 - xr - yr)/yr
        Xg, Yg, Zg = xg/yg, 1.0, (1.0 - xg - yg)/yg
        Xb, Yb, Zb = xb/yb, 1.0, (1.0 - xb - yb)/yb

        X = np.array([[Xr, Xg, Xb],
                      [Yr, Yg, Yb],
                      [Zr, Zg, Zb]], dtype=np.float64)

        S = np.linalg.solve(
            X,
            np.array([[xyzw[0]],
                      [xyzw[1]],
                      [xyzw[2]]], dtype=np.float64)
        )

        S.shape = (1, S.size)
        M = X*S

        return M

    def to_xyz(self, rgb: np.ndarray, observer: Observer,
               normalize: bool or float = True) -> np.ndarray:
        '''
        Transforms normalized linear R, G, B color coordinates to X, Y, Z
        coordinates.

        Parameters
        ----------
        rgb: np.ndarray
            Normalized linear coordinates R, G and B to transform. The
            (R, G, B) values must be stored in rows of the input array if
            multiple colors are transformed.
        normalize: bool
            If True, the transformation matrix is computed from the XYZ
            coordinates of the illuminant that are normalized
            by the illuminant luminosity Y as
            :math:`X_n, Y_n, Z_n = X/Y, 1.0, Z/Y`.
            If a floating-point value, the transformation matrix is computed
            from the XYZ coordinates of the illuminant the are weighted
            by the given value.

        Returns
        -------
        xyz: np.ndarray
            X, Y, Z coordinates stored in rows.
        '''
        T = self.to_xyz_transformation(observer, normalize=normalize)

        return np.dot(T, rgb)

    def from_xyz_transformation(self, observer: Observer,
                                normalize: bool or float = True) -> np.ndarray:
        '''
        Compute transformation from XYZ coordinates of the specified
        observer to normalized linear RGB coordinates.

        Parameters
        ----------
        observer: Observer
            Target observer.
        normalize: bool or float
            If True, the transformation matrix is computed from the XYZ
            coordinates of the illuminant that are normalized
            by the illuminant luminosity Y as
            :math:`X_n, Y_n, Z_n = X/Y, 1.0, Z/Y`.
            If a floating-point value, the transformation matrix is computed
            from the XYZ coordinates of the illuminant the are weighted
            by the given value.

        Returns
        -------
        T: np.ndarray
            A 3x3 transformation matrix as a numpy array.
        '''
        return np.linalg.inv(
            self.to_xyz_transformation(observer, normalize=normalize))

    def from_xyz(self, xyz: np.ndarray, observer: Observer,
                 normalize: bool or float = True) -> np.ndarray:
        '''
        Transforms X, Y, Z color coordinates computed with the specified
        observer to normalized R, G, B coordinates.

        Parameters
        ----------
        xyz: np.ndarray
            Coordinates X, Y and Z to transform.  The (X, Y, Z) values must be
            stored in rows of the input array if multiple colors are
            transformed.
        normalize: bool or float
            If True, the transformation matrix is computed from the XYZ
            coordinates of the illuminant that are normalized
            by the illuminant luminosity Y as
            :math:`X_n, Y_n, Z_n = X/Y, 1.0, Z/Y`.
            If a floating-point value, the transformation matrix is computed
            from the XYZ coordinates of the illuminant the are weighted
            by the given value.

        Returns
        -------
        rgb: np.ndarray
            Normalized RGB coordinates stored in rows.
        '''
        T = self.from_xyz_transformation(observer, normalize=normalize)

        return np.dot(T, xyz)

    def from_spectrum(self, wavelengths: np.ndarray, spectrum: np.ndarray,
                      observer: Observer, companding: bool = True,
                      normalize: bool or float = True):
        '''
        A convenience method that first computes XYZ coordinates and from there
        linear RGB coordinates. If companding is True, nonlinear transformation
        is applied to the linear RGB components, which yields the "standard"
        nonlinear RGB components normalized to range :math:`[0, 1]`.

        Parameters
        ----------
        wavelengths: np.ndarray
            Wavelengths od light (m) at which the spectrum intensities are
            defined.
        spectrum: np.ndarray
            Spectrum defined at specified intensities.
        companding: bool
            If True, nonlinear transformation is applied to the linear
            RGB components, which yields the "standard" nonlinear RGB
            components normalized to range :math:`[0, 1]`.
        normalize: bool or float
            If True, the transformation matrix from XYZ to RGB is computed
            from the XYZ coordinates of the illuminant that are normalized
            by the illuminant luminosity Y as
            :math:`X_n, Y_n, Z_n = X/Y, 1.0, Z/Y`. If True, the XYZ
            components computed from the spectrum are also normalized by
            the luminosity Y of the illuminant.
            If a floating-point value, the transformation matrix is computed
            from the XYZ coordinates of the illuminant the are weighted
            by the given value.

        Returns
        -------
        rgb: np.ndarray
            Normalized RGB components as a numpy vector of length 3. 
        '''
        xyz = observer.XYZ(wavelengths, spectrum)
        rgb = self.from_xyz(xyz, observer, normalize=normalize)
        if companding:
            rgb = self.apply_companding(rgb)
        return rgb

    def apply_companding(self, linear_rgb: np.ndarray):
        '''
        Apply nonlinear companding function to normalized linear RGB
        coordinates.

        Parameters
        ----------
        linear_rgb: np.ndarray
            Normalized linear RGB coordinates in rows.

        Returns
        -------
        rgb: np.ndarray
            Normalized transformed RGB coordinates.  
        '''
        companding = self._data['companding']
        if companding == 'srgb':
            result = self.srgb_companding(linear_rgb)
        elif companding == 'l*':
            result = self.l_companding(linear_rgb)
        else:
            gamma = self._data['gamma']
            result = self.gamma_companding(linear_rgb, gamma=gamma)

        return result

    def _get_name(self) -> str:
        return self._pretty_name
    name = property(_get_name, None, None, 'Name of the RGB color space.')

    def __str__(self):
        return '{}'.format(self._pretty_name)

    def __repr__(self):
        return self.__str__()

# creating instances of implemented RGB color spaces
AppleRGB = RGB('apple', RGB_COLOR_SPACES['apple'],
               pretty_name='Apple RGB')
BestRGB = RGB('best', RGB_COLOR_SPACES['best'],
              pretty_name='Best RGB')
BetaRGB = RGB('beta', RGB_COLOR_SPACES['beta'],
              pretty_name='Beta RGB')
BruceRGB = RGB('bruce', RGB_COLOR_SPACES['bruce'],
               pretty_name='Bruce RGB')
CieRGB = RGB('cie', RGB_COLOR_SPACES['cie'],
             pretty_name='CIE RGB')
ColorMatchRGB = RGB('color_match', RGB_COLOR_SPACES['color_match'],
                    pretty_name='ColorMatch RGB')
DonRGB = RGB('don', RGB_COLOR_SPACES['don'],
             pretty_name='DON RGB')
EciRGB = RGB('eci_v2', RGB_COLOR_SPACES['eci_v2'], pretty_name='ECI RGB')
EktaSpacePs5RGB = RGB('ekta_space_ps5', RGB_COLOR_SPACES['ekta_space_ps5'],
                      pretty_name='Ekta Space PS5 RGB')
NtscRGB = RGB('ntsc', RGB_COLOR_SPACES['ntsc'],
              pretty_name='NTSC RGB')
PalSecamRGB = RGB('pal_secam', RGB_COLOR_SPACES['pal_secam'],
                  pretty_name='PAL/SECAM RGB')
ProPhotoRGB = RGB('prophoto', RGB_COLOR_SPACES['prophoto'],
                  pretty_name='ProPhoto RGB')
SmpteCRGB = RGB('smpte-c', RGB_COLOR_SPACES['smpte-c'],
                pretty_name='SMPTE-C RGB')
SRGB = RGB('srgb', RGB_COLOR_SPACES['srgb'], pretty_name='SRGB')
WideGamutRGB = RGB('wide-gamut', RGB_COLOR_SPACES['wide_gamut'],
                   pretty_name='WideGamut RGB')


class Illuminant:
    def __init__(self, estimator, name: str = 'Illuminant'):
        self._estimator = estimator
        self._name = name

    def spd(self, wavelengths: None or float or np.ndarray = None) -> np.ndarray:
        '''
        Returns the spectral power distribution of the illuminant at the
        specified wavelengths. If no wavelengths are specified,
        the values at standard wavelengths :py:attr:`STANDARD_WAVELENGTHS`
        are returned.

        Parameters
        ----------
        wavelengths: None or float or np.ndarray
            Wavelengths of light (m) at which to compute the spectral power
            distribution of the illuminant. If no wavelengths are specified,
            the values at standard wavelengths :py:attr:`STANDARD_WAVELENGTHS`
            are returned. 
        '''
        if wavelengths is None:
            wavelengths = STANDARD_WAVELENGTHS

        return self._estimator(wavelengths)

    def raw_spd(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns raw spectral power distribution and corresponding wavelengths.

        Returns
        -------
        wavelengths: np.ndarray
            Wavelengths at which the raw spectral power densities are defined.
        spd: np.ndarray
            Raw spectral power densities defined at the returned wavelengths. 
        '''
        return self._estimator()

    def XYZ(self, observer: str or Observer, normalize: bool or float = False):
        '''
        Returns the X, Y, Z color components of the illuminant for the
        specified observer. Luminosity Y is not normalized by default.

        Parameters
        ----------
        observer: str or Observer
            Selected standard observer :py:class`CIE1931` ("cie1931") or
            :py:class;`CIE1964` ("cie1964").

        normalize: bool or float
            If True, the computed X, Y, and Z components are divided by the
            luminance Y (the returned components will be X/Y, 1.0, Z/Y).
            If a floating-point value, the computed X, Y and Z components are
            weighted with the provided value.

        Returns
        -------
        XYZ: np.ndarray vector of length 3
            X, Y and Z component of the standard source computed for the
            given standard observer.
        '''
        if isinstance(observer, str):
            observer_ = {'cie1931': CIE1931, 'cie1964': CIE1964}.get(observer)
            if observer_ is None:
                raise ValueError('Unknown observer "{}"!'.format(observer))
            observer = observer_

        XYZ = simps(
            observer.xyz_cmf(STANDARD_WAVELENGTHS)*
                self.spd(STANDARD_WAVELENGTHS),
            dx=STANDARD_WAVELENGTH_STEP
        )

        if isinstance(normalize, (bool)) and normalize:
            XYZ = XYZ/XYZ[1]
        elif isinstance(normalize, float):
            XYZ *= normalize
    
        return XYZ

    def xyY(self, observer: Observer) -> Tuple[float, float, float]:
        '''
            Compute CIE x, y chromaticities and Y luminosity of the illuminant
            for the specified observer.

            Parameters
            ----------
            observer: Observer or str
                Color-matching functions that are used to compute the x and y
                chromaticities and luminosity Y.
                Should be one of "1931" or "1964" for CIE 1931
                (2 deg observer) or CIE 1964 (10 deg observer), respectively.
                Default value is set to "1931".

            Returns
            -------
            x,: float
                Chromaticity y of the illuminant: math:`y=\\frac{Y}{X + Y + Z}`.
            y: float
                Chromaticity x of the illuminant: math:`x=\\frac{X}{X + Y + Z}`.
            Y: float
                Luminosity of the illuminant.
        '''
        X, Y, Z = self.XYZ(observer)
        sum_XYZ = X + Y + Z

        return X/sum_XYZ, Y/sum_XYZ, Y

    def Y(self, observer: Observer) -> float:
        '''
            Compute CIE luminosity (Y) of the illuminant for the specified
            observer.

            Parameters
            ----------
            observer: Observer or str
                Color-matching functions that are used to compute the x and y
                chromaticities and luminosity Y.
                Should be one of "1931" or "1964" for CIE 1931
                (2 deg observer) or CIE 1964 (10 deg observer), respectively.
                Default value is set to "1931".

            Returns
            -------
            Y: float
                Luminosity of the illuminant.

            Note
            ----
            Luminosity of the illuminant depends on its spectral power
            density and the color/luminosity matching/sensitivity of
            the observer. 
        '''
        return self.XYZ(observer, normalize=False)[1]

    def transformation_to_illuminant(
            self, illuminant: 'Illuminant', observer: Observer) -> np.ndarray:
        '''
        Compute a matrix that will map X, Y, Z from this illuminant to the
        specified target illuminant.

        Parameters
        ----------
        illuminant: Illuminant
            Target illuminant of the transformation/map.
        observer: Observer
            Selected standard observer :py:class:`CIE1931` or
            :py:class:`CIE1964`.

        Returns
        -------
        M: numpy.ndarray
            A linear transformation that maps X, Y, Z coordinates to
            the target illuminant, i.e. `XYZ2 = np.dot(M, XYZ_this)`
        '''
        XYZw1 = self.XYZ(observer, normalize=True)
        XYZw2 = illuminant.XYZ(observer, normalize=True)

        M = np.diag(XYZw2/XYZw1)

        return M

    def _get_name(self) -> str:
        return self._name
    name = property(_get_name, None, None, 'Illuminant name.')

    def __str__(self):
        return '{} illuminant'.format(self._name)

    def __repr__(self):
        return self.__str__()


D50 = Illuminant(_d50_source_estimator, name='D50')
''' Instance of standard illuminant D50. '''

D55 = Illuminant(_d55_source_estimator, name='D55')
''' Instance of standard illuminant D55. '''

D65 = Illuminant(_d65_source_estimator, name='D65')
''' Instance of standard illuminant D65. '''

D75 = Illuminant(_d75_source_estimator, name='D75')
''' Instance of standard illuminant D75. '''

A = Illuminant(_a_source_estimator, name='A')
''' Instance of standard illuminant A. '''

C = Illuminant(_c_source_estimator, name='C')
''' Instance of standard illuminant C. '''

E = Illuminant(_e_source_estimator, name='E')
''' Instance of standard illuminant E. '''

FL1 = Illuminant(_fl_1_source_estimator, name='FL1')
''' Fluorescence lamp source FL1.  '''
FL2 = Illuminant(_fl_2_source_estimator, name='FL2')
''' Fluorescence lamp source FL2.  '''
FL3 = Illuminant(_fl_3_source_estimator, name='FL3')
''' Fluorescence lamp source FL3.  '''
FL4 = Illuminant(_fl_4_source_estimator, name='FL4')
''' Fluorescence lamp source FL4.  '''
FL5 = Illuminant(_fl_5_source_estimator, name='FL5')
''' Fluorescence lamp source FL5.  '''
FL6 = Illuminant(_fl_6_source_estimator, name='FL6')
''' Fluorescence lamp source FL6.  '''
FL7 = Illuminant(_fl_7_source_estimator, name='FL7')
''' Fluorescence lamp source FL7.  '''
FL8 = Illuminant(_fl_8_source_estimator, name='FL8')
''' Fluorescence lamp source FL8.  '''
FL9 = Illuminant(_fl_9_source_estimator, name='FL9')
''' Fluorescence lamp source FL9.  '''
FL10 = Illuminant(_fl_10_source_estimator, name='FL10')
''' Fluorescence lamp source FL10.  '''
FL11 = Illuminant(_fl_11_source_estimator, name='FL11')
''' Fluorescence lamp source FL11.  '''
FL12 = Illuminant(_fl_12_source_estimator, name='FL12')
''' Fluorescence lamp source FL12.  '''

FL3_1 = Illuminant(_fl3_1_source_estimator, name='FL3.1')
''' Fluorescence lamp source FL3.1.  '''
FL3_2 = Illuminant(_fl3_2_source_estimator, name='FL3.2')
''' Fluorescence lamp source FL3.2.  '''
FL3_3 = Illuminant(_fl3_3_source_estimator, name='FL3.3')
''' Fluorescence lamp source FL3.3.  '''
FL3_4 = Illuminant(_fl3_4_source_estimator, name='FL3.4')
''' Fluorescence lamp source FL3.4.  '''
FL3_5 = Illuminant(_fl3_5_source_estimator, name='FL3.5')
''' Fluorescence lamp source FL3.5.  '''
FL3_6 = Illuminant(_fl3_6_source_estimator, name='FL3.6')
''' Fluorescence lamp source FL3.6.  '''
FL3_7 = Illuminant(_fl3_7_source_estimator, name='FL3.7')
''' Fluorescence lamp source FL3.7.  '''
FL3_8 = Illuminant(_fl3_8_source_estimator, name='FL3.8')
''' Fluorescence lamp source FL3.8.  '''
FL3_9 = Illuminant(_fl3_9_source_estimator, name='FL3.9')
''' Fluorescence lamp source FL3.9.  '''
FL3_10 = Illuminant(_fl3_10_source_estimator, name='FL3.10')
''' Fluorescence lamp source FL3.10.  '''
FL3_11 = Illuminant(_fl3_11_source_estimator, name='FL3.11')
''' Fluorescence lamp source FL3.11.  '''
FL3_12 = Illuminant(_fl3_12_source_estimator, name='FL3.12')
''' Fluorescence lamp source FL3.12.  '''
FL3_13 = Illuminant(_fl3_13_source_estimator, name='FL3.13')
''' Fluorescence lamp source FL3.13.  '''
FL3_14 = Illuminant(_fl3_14_source_estimator, name='FL3.14')
''' Fluorescence lamp source FL3.14.  '''
FL3_15 = Illuminant(_fl3_15_source_estimator, name='FL3.15')
''' Fluorescence lamp source FL3.15.  '''

def standard_illuminant(illuminant: str) -> Illuminant:
    '''
    Returns standard illuminant instance.

    Parameters
    ----------
    illuminant: str
        Standard illuminant "d50", "d55", "d65", "d75", "a", "c", or "e".

    Returns
    -------
        The requested standard illuminant instance.
    '''
    obj = {
        'd50': D50,
        'd55': D55,
        'd65': D65,
        'd75': D75,
        'a': A,
        'c': C,
        'e': E,
    }.get(str(illuminant).lower())

    if obj is None:
        raise ValueError('Unknown illuminant "{}"!'.format(illuminant))

    return obj

if __name__ == '__main__':
    from matplotlib import pyplot as pp

    pp.figure()
    pp.title('Standard illuminants')
    pp.plot(D50.raw_spd()[0]*1e9, D50.raw_spd()[1], label='D50')
    pp.plot(D55.raw_spd()[0]*1e9, D55.raw_spd()[1], label='D55')
    pp.plot(D65.raw_spd()[0]*1e9, D65.raw_spd()[1], label='D65')
    pp.plot(D75.raw_spd()[0]*1e9, D75.raw_spd()[1], label='D75')
    pp.plot(A.raw_spd()[0]*1e9, A.raw_spd()[1], label='A')
    pp.plot(C.raw_spd()[0]*1e9, C.raw_spd()[1], label='C')
    pp.plot(E.raw_spd()[0]*1e9, E.raw_spd()[1], label='E')
    pp.grid()
    pp.xlabel('Wavelength (nm)')
    pp.ylabel('Spectral power density (a.u.)')
    pp.legend()

    pp.figure()
    pp.title('Fluorescence lamp sources F1-F12')
    for i in range(12):
        fl = globals().get('FL{}'.format(i + 1))
        pp.plot(fl.raw_spd()[0]*1e9, fl.raw_spd()[1], label='FL{}'.format(i + 1))
    pp.xlabel('Wavelength (nm)')
    pp.ylabel('Spectral power density (a.u.)')
    pp.legend()

    pp.figure()
    pp.title('Fluorescence lamp sources F3.1-F3.15')
    for i in range(15):
        fl = globals().get('FL3_{}'.format(i + 1))
        pp.plot(fl.raw_spd()[0]*1e9, fl.raw_spd()[1], label='FL3.{}'.format(i + 1))
    pp.xlabel('Wavelength (nm)')
    pp.ylabel('Spectral power density (a.u.)')
    pp.legend()

    pp.figure()
    pp.title('CIE 1931 XYZ color-matching functions')
    pp.plot(CIE1931.xyz_cmf_raw()[0]*1e9, CIE1931.xyz_cmf_raw()[1][0],
            '-r', label='X')
    pp.plot(CIE1931.xyz_cmf_raw()[0]*1e9, CIE1931.xyz_cmf_raw()[1][1],
            '-g', label='Y')
    pp.plot(CIE1931.xyz_cmf_raw()[0]*1e9, CIE1931.xyz_cmf_raw()[1][2],
            '-b', label='Z')
    pp.grid()
    pp.xlabel('Wavelength (nm)')
    pp.legend()

    pp.figure()
    pp.title('CIE 1964 XYZ color-matching functions')
    pp.plot(CIE1964.xyz_cmf_raw()[0]*1e9, CIE1964.xyz_cmf_raw()[1][0],
            '-r', label='X')
    pp.plot(CIE1964.xyz_cmf_raw()[0]*1e9, CIE1964.xyz_cmf_raw()[1][1],
            '-g', label='Y')
    pp.plot(CIE1964.xyz_cmf_raw()[0]*1e9, CIE1964.xyz_cmf_raw()[1][2],
            '-b', label='Z')
    pp.grid()
    pp.xlabel('Wavelength (nm)')
    pp.legend()

    pp.figure()
    pp.title('CIE 1931 RGB color-matching functions')
    pp.plot(CIE1931.rgb_cmf_raw()[0]*1e9, CIE1931.rgb_cmf_raw()[1][0],
            '-r', label='R')
    pp.plot(CIE1931.rgb_cmf_raw()[0]*1e9, CIE1931.rgb_cmf_raw()[1][1],
            '-g', label='G')
    pp.plot(CIE1931.rgb_cmf_raw()[0]*1e9, CIE1931.rgb_cmf_raw()[1][2],
            '-b', label='B')
    pp.grid()
    pp.xlabel('Wavelength (nm)')
    pp.legend()

    pp.figure()
    pp.title('CIE 1964 RGB color-matching functions')
    pp.plot(CIE1964.rgb_cmf_raw()[0]*1e9, CIE1964.rgb_cmf_raw()[1][0],
            '-r', label='R')
    pp.plot(CIE1964.rgb_cmf_raw()[0]*1e9, CIE1964.rgb_cmf_raw()[1][1],
            '-g', label='G')
    pp.plot(CIE1964.rgb_cmf_raw()[0]*1e9, CIE1964.rgb_cmf_raw()[1][2],
            '-b', label='B')
    pp.grid()
    pp.xlabel('Wavelength (nm)')
    pp.legend()

    pp.show()
