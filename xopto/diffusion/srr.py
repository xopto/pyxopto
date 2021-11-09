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

def fresnel_s_polarized(theta:float, n1: float, n2: float) -> float:
    '''
    Fresnel's reflectance for S-polarized light at a refractive index boundary
    n1 => n2.

    Parameters
    ----------
    n1: float
        Refractive index of the incidence medium.
    n2: float
        Refractive index of the medium across the boundary.
    theta: float
        Incidence angle (rad).

    Returns
    -------
    reflectance: float
        Fresnel's reflectance for S-polarized light.
    '''
    theta = np.array(theta)
    frs = lambda theta: (
        (n1*np.cos(theta) - n2*np.sqrt(1.0 - (n1/n2*np.sin(theta))**2))/
        (n1*np.cos(theta) + n2*np.sqrt(1.0 - (n1/n2*np.sin(theta))**2)))**2

    if n1/n2 > 1.0:
        theta_total = np.arcsin(n2/n1)
        mask = theta < theta_total
        return np.concatenate(
            (frs(theta[mask]), np.ones_like(theta[np.logical_not(mask)])))
    else:
        return frs(theta)


def fresnel_p_polarized(theta: float, n1: float, n2: float) -> float:
    '''
    Fresnel's reflectance for P-polarized light at a refractive index boundary
    n1 => n2.

    Parameters
    ----------
    n1: float
        Refractive index of the incidence medium.
    n2: float
        Refractive index of the medium across the boundary.
    theta: float
        Incidence angle (rad).

    Returns
    -------
    reflectance: float
        Fresnel's reflectance for P-polarized light.
    '''
    theta = np.array(theta)
    frp = lambda theta: (
        (n1*np.sqrt(1.0 - (n1/n2*np.sin(theta))**2) - n2*np.cos(theta))/
        (n1*np.sqrt(1.0 - (n1/n2*np.sin(theta))**2) + n2*np.cos(theta)))**2

    if n1/n2 > 1.0:
        theta_total = np.arcsin(n2/n1)
        mask = theta < theta_total
        return np.concatenate(
            (frp(theta[mask]), np.ones_like(theta[np.logical_not(mask)])))
    else:
        return frp(theta)


def fresnel_unpolarized(theta: float, n1: float, n2: float) -> float:
    '''
    Fresnel's reflectance for unpolarized light at a refractive index boundary
    n1 => n2.

    Parameters
    ----------
    n1: float
        Refractive index of the incidence medium.
    n2: float
        Refractive index of the medium across the boundary.
    theta: float
        Incidence angle (rad).

    Returns
    -------
    reflectance: float
        Fresnel's reflectance for unpolarized light.
    '''
    return 0.5*(fresnel_s_polarized(theta, n1, n2) +
                fresnel_p_polarized(theta, n1, n2))


def _fresnel_phi(theta: float, nmedium: float, noutside: float) -> float:
    return 2.0*np.sin(theta)*\
           np.cos(theta)*\
           fresnel_unpolarized(theta, nmedium, noutside)


def _fresnel_j(theta: float, nmedium: float, noutside: float) -> float:
    return 3.0*np.sin(theta)*\
           np.cos(theta)**2*\
           fresnel_unpolarized(theta, nmedium, noutside)


def r_eff(nmedium: float, noutside: float, steps: int = 10000) -> float:
    theta = np.linspace(0, np.pi/2.0, steps)

    r_phi = simps(_fresnel_phi(theta, nmedium, noutside), dx=theta[1] - theta[0])
    r_j = simps(_fresnel_j(theta, nmedium, noutside), dx=theta[1] - theta[0])

    return (r_phi + r_j)/(2.0 - r_phi + r_j)


class SRDA:
    def __init__(self, nmedium: float, noutside: float,
                 acceptance: float, steps: int = 10000):
        '''
        Prepares the diffusion approximation solution for spatially resolved
        reflectance.

        Parameters
        ----------
        nmedium: float
            Refractive index of the medium.
        noutside: float
            Refractive index above the medium.
        acceptance: float
            Acceptance (half) angle of the detector (radians).
        steps: int
            The number of integration steps to calculate coefficients in the
            diffusion approximation.
        '''

        theta = np.linspace(0, acceptance, steps)

        back_hemi_coeff1 = 0.5*(
            1.0 - fresnel_unpolarized(theta, nmedium, noutside))*\
            np.cos(theta)*np.sin(theta)

        self._back_hemi_coeff1 = simps(back_hemi_coeff1, dx=theta[1]-theta[0])

        back_hemi_coeff2 = 1.5*(
            1.0 - fresnel_unpolarized(theta, nmedium, noutside))*\
            np.cos(theta)**2*np.sin(theta)

        self._back_hemi_coeff2 = simps(back_hemi_coeff2, dx=theta[1] - theta[0])

        self._R_eff = r_eff(nmedium, noutside, steps)

    def __call__(self, r: np.ndarray, mua: float, musr: float) -> np.ndarray:
        '''
        Calculates the reflectance profile for selected optical properties.

        Parameters
        ----------
        r: np.ndarray
            Radial points at which the reflectance is evaluated.
        mua: float
            Absorption coefficient of the medium in 1/m.
        musr: float
            Reduced scattering coefficient of the medium in 1/m.

        Returns
        -------
        reflectance: np.ndarray
            Reflectance profile evaluated at r.
        '''

        mueff = np.sqrt(3.0*mua*(mua + musr))
        z0 = 1.0/(mua + musr)
        D = 1.0/(3.0*(mua + musr))
        zB = 2.0*D*(1.0 + self._R_eff)/(1.0 - self._R_eff)

        return self._back_hemi_coeff1 * self._phi(r, mueff, z0, zB, D) + \
               self._back_hemi_coeff2 * self._Rf(r, mueff, z0, zB)

    def _Rf(self, r, mueff, z0, zB):
        r1 = np.sqrt(z0**2 + r**2)
        r2 = np.sqrt((z0 + 2.0*zB)**2 + r**2)

        return 1.0/(4.0*np.pi)*(
            z0*(mueff + 1.0/r1)*np.exp(-mueff*r1)/r1**2 +
            (z0 + 2.0*zB)*(mueff + 1.0/r2)*np.exp(-mueff*r2)/r2**2)

    def _phi(self, r, mueff, z0, zB, D):
        r1 = np.sqrt(z0**2 + r**2)
        r2 = np.sqrt((z0 + 2.0*zB)**2 + r**2)

        return 1.0/(4.0*np.pi*D)*(np.exp(-mueff*r1)/r1 - np.exp(-mueff*r2)/r2)
