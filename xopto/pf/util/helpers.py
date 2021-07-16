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
import numpy as np


def g2gamma(g1:float, g2:float) -> float:
    '''
    Computes and return gamma from the first and second Legendre moment.

        gamma = (1 - g2)/(1 - g1)

    Parameters
    ----------
    g1: float
        First Legendre moment of the scattering phase function.
    g2: float
        Second Legendre moment of the scattering phase function.

    Returns
    -------
    gamma: float
        Similarity parameter gamma.
    '''
    return (1.0 - g2)/(1.0 - g1)


def g2delta(g1: float, g3: float) -> float:
    '''
    Computes and delta from the first and third Legendre moment.

        delta = (1 - g3)/(1 - g1)

    Parameters
    ----------
    g1: float
        First Legendre moment of the scattering phase function.
    g3: float
        Third Legendre moment of the scattering phase function.

    Returns
    -------
    delta: float
        Similarity parameter delta.
    '''
    return (1.0 - g3)/(1.0 - g1)


def gamma2g2(gamma: float, g1:float) -> float:
    '''
    Computes the second Legendre moment from gamma and the first Legendre
    moment of the scattering phase function.

        gamma = (1 - g2)/(1 - g1)

        g2 = 1 - gamma*(1 - g1)

    Parameters
    ----------
    g1: float
        First Legendre moment of the scattering phase function.
    gamma: float
        Parameter gamma.

    Returns
    -------
    g2: float
        Second Legendre moment of the scattering phase function.
    '''
    return 1.0 - gamma*(1.0 - g1)


def delta2g3(delta: float, g1: float) -> float:
    '''
    Computes the third Legendre moment from delta and the first Legendre
    moment of the scattering phase function.

        delta = (1 - g3)/(1 - g1)

        g3 = 1 - delta*(1 - g1)

    Parameters
    ----------
    g1: float
        First Legendre moment of the scattering phase function.
    gamma: float
        Parameter gamma.

    Returns
    -------
    g3: float
        Third Legendre moment of the scattering phase function.
    '''
    return 1.0 - delta*(1.0 - g1)


def gammadelta2g(g1: float, gamma: float, delta: float) \
            -> Tuple[float, float, float]:
    '''
    Computes the second and third Legendre moment from gamma and delta
    parameters.
    Returns (g1, g2, g3).

        gamma = (1 - g2)/(1 - g1)

        delta = (1 - g3)/(1 - g1)

        g2 = 1 - gamma*(1 - g1)

        g3 = 1 - delta*(1 - g1)

    Parameters
    ----------
    g1: float
        First Legendre moment of the scattering phase function.

    gamma: float
        Parameter gamma.

    delta: float
        Parameter delta.

    Returns
    -------
    g: tuple of 3 elements
        The first three Legendre moments (g1, g2, g3).
    '''
    return g1, 1.0 - gamma*(1.0 - g1), 1.0 - delta*(1.0 - g1)


def g2gammadelta(g1: float, g2: float, g3: float) -> Tuple[float, float, float]:
    '''
    Computes the gamma and delta from the first three Legendre moments of the
    scattering phase function. Returns (g1, gamma, delta).

        gamma = (1 - g2)/(1 - g1)

        delta = (1 - g3)/(1 - g1)

    Parameters
    ----------
    g1: float
        First Legendre moment of the scattering phase function.

    gamma: float
        Parameter gamma.

    delta: float
        Parameter delta.

    Returns
    -------
    g: tuple of 3 elements
        The first three Legendre moments (g1, gamma, delta).
    '''
    return g1, (1.0 - g2)/(1.0 - g1), (1.0 - g3)/(1.0 - g1)


def sigma(gs: list or tuple) -> float:
    '''
    Computes parameter sigma from the given Legendre moments.

        sigma = sum( (-0.5)^(i - 2) (1 - gi)/(1 - g1), for i = 2, 3, ...

        where gi is the i-th Legendre moment.

    References:

    Bodenschatz et al, JBO, 21(3) 2016.

    Parameters
    ----------
    gs: tuple, list, ndarray
        First n Legendre moments as (g0, g1, g2, g3, ..., gn)

    Returns
    -------
    sigma: float
        Parameter sigma
    '''
    gs = np.asarray(gs, dtype=np.float64)
    i = np.arange(1, gs.size - 1)
    return ((-0.5)**(i + 1 - 2)*(1 - gs[2:]) / (1 - gs[1])).sum()


def mhg_inverse_g1gamma(
        g1: float or np.ndarray, gamma: float or np.ndarray) \
            -> Tuple[float or np.ndarray, float or np.ndarray, bool or np.ndarray]:
    '''
    Estimates the parameters g and b of the modified Henyey-Greenstein (MHG)
    scattering phase function that best fit the given values of g1 and gamma,
    where g1 is the first Legendre moment of the MHG phase function and gamma:

        gamma = (1 - g2) / (1 - g1).

    Parameters g and b of the scattering phase function MHG can be computed
    analytically.

    Forward analytical model (Calabro et al. JBO 19(7) 2014):

    g1 = b * g

    g2 = b * g^2 + 2/5 (1 - b)

    g3 = b * g^3

    g_k = b * g^k    where   k > 2

    ...

    Parameters g and b can be derived from a quadratic equation. The
    obtained values are valid only if:

    b < 1  (if b == 1 ... MHG simplifies to HG) or gamma < 1 + g1.

    Parameters
    ---------
    g1: float
        Target value of the first Legendre moment.

    gamma: float
        Target value of parameter gamma = (1 - g2)/(1 - g1).

    Returns
    -------
    g: float
        Estimated values of parameter g.

    b: float
        Estimated values of parameter b.

    valid: bool
        A binary mask indicating valid values of g and b.
    '''
    g1 = np.asarray(g1, dtype=np.float).flatten()
    gamma = np.asarray(gamma, dtype=np.float).flatten()

    n = max(g1.size, gamma.size)
    g_hg = np.zeros((n,))
    b = np.ones((n,))

    # special case for g1 = 0
    g1_zero_mask = g1 == 0.0
    b[g1_zero_mask] = 2.5*gamma - 1.5
    g_hg[g1_zero_mask] = 0.0

    g1_nonzero_mask = np.logical_not(g1_zero_mask)

    #inverse is root of quadratic equation: k_a*g^2 + k_b*g + k_c = o
    g1_nonzero = g1[g1_nonzero_mask]
    k_a = g1_nonzero
    k_b = gamma*(1.0 - g1_nonzero) - 0.6
    k_c = -0.4*g1_nonzero
    D = k_b*k_b - 4.0*k_a*k_c

    # only b1 is a valid solution - that is verified numerically
    # (b2 is always negative or zero)
    g_hg[g1_nonzero_mask] = (-k_b + np.sqrt(D))/(2.0*k_a)
    b[g1_nonzero_mask] = g1/g_hg[g1_nonzero_mask]

    valid_values = b < 1

    return g_hg, b, valid_values

def mhg_inverse_gammadelta(
        gamma: float or np.ndarray, delta: float or np.ndarray) \
            -> Tuple[float or np.ndarray, float or np.ndarray, bool or np.ndarray]:
    '''
    Estimates the parameters g and b of the modified Henyey-Greenstein (MHG)
    scattering phase function that best fit the given values of gamma and delta,
    where gamma is

        gamma = (1 - g2) / (1 - g1),

    and delta is:

        delta = (1 - g3) / (1 - g1),

    Parameters g and b of the scattering phase function MHG can be computed
    analytically.

    Forward analytical model (Calabro et al. JBO 19(7) 2014):

    g1 = b * g

    g2 = b * g^2 + 2/5 (1 - b)

    g3 = b * g^3

    g_k = b * g^k    where   k > 2

    ...

    Parameters g and b can be derived from a third-order equation. The
    obtained values are valid only if 0 <= b <= 1.

    Parameters
    ---------
    gamma: float or np.ndarray
        Target value of parameter gamma = (1 - g2)/(1 - g1).

    delta: float or np.ndarray
        Target value of parameter delta = (1 - g3)/(1 - g1).

    Returns
    -------
    g: float or np.ndarray
        Estimated values of parameter g.

    b: float or np.ndarray
        Estimated values of parameter b.

    valid: bool or np.ndarray
        A binary mask indicating valid values of g and b.
    '''
    gamma = np.asarray(gamma, dtype=np.float).flatten()
    delta = np.asarray(delta, dtype=np.float).flatten()

    fp_tol = 1e-6

    n = max(gamma.size, delta.size)
    g_hg = np.zeros((n,))
    b = np.ones((n,))

    for index, pair in enumerate(zip(gamma, delta)):
        gamma_value, delta_value = pair

        # special handling of delta = 1.0
        if np.abs(delta_value - 1.0) < fp_tol:
            if np.abs(gamma_value - 0.6) < fp_tol:
                # b is not unique ... can be any value from [0, 1]
                b[index] = 0.0
                g_hg[index] = 1.0
            else:
                g_hg[index] = 0.0
                b[index] = 2.5*gamma_value - 1.5
        else:

            #inverse is root of quadratic equation: k_a*g^3 + k_b*g^2 + k_c*g + k_d = o
            k_a = 0.6 - gamma_value
            k_b = delta - 1.0
            k_c = gamma_value - 0.6*delta_value
            k_d = 0.4 - 0.4*delta_value

            # print('Coefficients:', k_a, k_b, k_c, k_d)
            roots = np.roots((k_a, k_b, k_c, k_d))
            # print('Roots:', roots)
            valid_roots = roots[np.isreal(roots)]
            g_hg_value = valid_roots[
                np.logical_and(valid_roots >= 0.0, valid_roots < 1.0 - fp_tol)
            ]
            if g_hg_value.size == 1:
                if abs(g_hg_value) <= fp_tol:
                    b[index] = 2.5*gamma_value - 1.0
                else:
                    b[index] = (delta_value - 1.0)/ \
                            (g_hg_value*(delta_value - g_hg_value**2))
                g_hg[index] = g_hg_value

    valid_values = b < 1

    return g_hg, b, valid_values
