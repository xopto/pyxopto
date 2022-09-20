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
from scipy.integrate import quad, simps
from scipy.optimize import minimize
from scipy.interpolate import interp1d

#%% Support functions
def fastg(n: int, pf: np.ndarray, costheta: np.ndarray = None) -> float:
    '''
    Computes the n-th Legendre moment of the scattering phase function
    defined on equally spaced grid of deflection angle cosines.
    Note that if the scattering phase function is defined only at
    a few points, the computed Legendre moments might be inaccurate.

    Parameters
    ----------
    n: int
        The Legendre moment to calculate. (:math:`n=1` yields the anisotropy
        factor :math:`g_1`).
    pf: np.ndarray
        A vector containing the values of the scattering phase function at
        equally spaced deflection angle cosines from [-1, 1].
    costheta: np.ndarray
        A vector of equally spaced numbers (deflection angle cosines) on
        the interval [-1, 1] at which the phase function is estimated. If
        None, costheta is calculated as np.linspace(-1.0, 1.0, pf.size)

    Returns
    -------
    p: float
        Value of the n-th Legendre moment.
    '''
    if costheta is None:
        costheta = np.linspace(-1.0, 1.0, pf.size)

    return simps(pf*np.polynomial.legendre.Legendre.basis(n)(costheta),
                 dx=costheta[1] - costheta[0])

def fastgs(last: int, pf: np.ndarray, costheta: np.ndarray = None,
           out: np.ndarray = None):
    '''
    Computes the first n Legendre moments of the scattering phase function
    defined on equally spaced grid of deflection angle cosines.
    Note that if the scattering phase function is defined only at
    a few points, the computed Legendre moments might be inaccurate.

    Parameters
    ----------
    last: int
        The last Legendre moment to be calculate.
    pf: np.ndarray
        A vector containing the values of the scattering phase function at
        equally spaced deflection angle cosines from [-1, 1].
    cotheta: np.ndarray
        A vector of equally spaced numbers (deflection angle cosines) on
        the interval [-1, 1] at which the phase function is estimated. If
        None, costheta is calculated as :code:`np.linspace(-1.0, 1.0, pf.size)`.
    out: np.ndarray
        A vector of size last + 1 for storing the Legendre moments.

    Returns
    -------
    p: np.ndarray
        Values of Legendre moments from 0 to 'last' as a vector
        :math:`([g_0, g_1, ..., g_{last}])`.
    '''
    if costheta is None:
        costheta = np.linspace(-1.0, 1.0, pf.size)

    if out is None:
        out = np.zeros([last + 1])

    for n in range(last + 1):
        out[n] = simps(pf*np.polynomial.legendre.Legendre.basis(n)(costheta),
                       dx=costheta[1] - costheta[0])
    return out


def lut_function(costheta: float or np.ndarray, params: list or tuple) \
        -> float or np.ndarray:
    '''
    Parametric approximation of the cumulative probability density function
    of the deflection angle cosine :math:`F(cos(\\theta))`:

    .. math::

        Fap(\\cos(\\theta)) = \\frac{a}{\\cos(\\theta) + b} + c,
        where Fap(-1) = 0 and Fap(1) = 1.

    Given the th
    wo constraints :math:`Fap(-1) = 0` and :math:`Fap(1) = 1`,
    the approximation function is fully defined by parameter :math:`c`:

    .. math::

        b = 2 c - 1
        a = 2 c (1 - c)

    Depending on the convexity of :math:`Fap`, :math:`c` can be from 
    :math:`[-inf, 0)` or :math:`(1, \\infty]`.

    Parameters
    ----------
    params: list or tuple
        Approximation function parameter list. This implementation requires
        only parameter c. Parameters a and b are calculated according
        to the above equations.
    costheta: float or np.ndarray
        A vector of deflection angle cosines at which the
        function is evaluated.

    Returns
    -------
    approx: np.ndarray
        Lut function approximation at the specified deflection angle cosines.

    Examples
    --------
    Approximation function for several different values of parameter a.

    >>> from matplotlib import pyplot as pp
    >>> import numpy as np

    >>> costheta = np.linspace(-1.0, 1.0, 1000)
    >>> C = [-0.01, -0.1, -1, -10]
    >>> for c in C:
            pp.plot(costheta, lut_function(costheta, [c]),
            label='a={}'.format(a))
    >>> pp.xlabel('cos(theta)')
    >>> pp.ylabel('F')
    >>> pp.legend()
    '''
    params = np.asarray(params)
    c = params[0]
    b = 2*c - 1
    a = 2*c*(1 - c)

    return a/(b + costheta) + c

def ilut_function(randnum: float or np.ndarray, params: list or tuple) \
        -> float or np.ndarray:
    '''
    Inverse of the parametric approximation of the cumulative probability
    density function Fap (python implementation lut_function):

    .. math::

        Fap^{-1}(F) &= a/(F - c) - b,

        where Fap^{-1}(-1) &= 0 and Fap^{-1}(1) = 1.

    Parameters
    ----------
    costheta: float or np.ndarray
        A vector of values (cumulative distribution) from [0.0, 1.0] for
        which the deflection angle cosine is computed.
    params: list or tuple
        Approximation function parameter list. This implementation requires
        only one parameter c. Parameters a and b are calculated according
        to the above equations.

    Returns
    -------
    approx: float or np.ndarray
        Deflection angle cosines approximations computed at the given
        values.

    Examples
    --------
    Inverse of the approximation function for several different values
    of parameter a.

    >>> from matplotlib import pyplot as pp
    >>> import numpy as np

    >>> F = np.linspace(0.0, 1.0, 1000)
    >>> C = [-0.01, -0.1, - 1, -10]
    >>> for c in C:
            pp.plot(F, ilut_function(F, [c]), label='a={}'.format(a))
    >>> pp.ylabel('cos(theta)')
    >>> pp.xlabel('F')
    >>> pp.legend()
    '''
    params = np.asarray(params)
    c = params[0]
    b = 2*c - 1
    a = 2*c*(1 - c)
    return a/(randnum - c) - b


class PfBase:
    '''
    The base class of all scattering phase functions.
    '''
    def __call__(self, costheta: float or np.ndarray) -> float or np.ndarray:
        '''
        Call method of the scattering phase function.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function
            is evaluated.

        Returns
        -------
        p: float or np.ndarray
            Scattering phase function at the specified scattering angle cosines.

        '''
        raise RuntimeError('The __call__ is not implemented!')

    def g(self, n: int, **kwargs) -> float:
        '''
        Computes the n-th Legendre moment of the scattering phase function.

        Parameters
        ----------
        n: int
            The Legendre moment to calculate. (n=1 is g_1 - anisotropy factor)
        kwargs: dict
            Optional keyword arguments passed to scipy.integrate.quad function.

        Returns
        -------
        p: float
            Value of the n-th Legendre moment.
        '''
        lp = np.polynomial.legendre.Legendre.basis(n)
        return quad(lambda x: self(x)*lp(x), -1.0, 1.0, **kwargs)[0]

    def gs(self, last: int, **kwargs) -> np.ndarray:
        '''
        Computes the first n Legendre moments of the scattering phase function.

        Parameters
        ----------
        last: int
            The last Legendre moment to be calculate.
        kwargs: dict
            Optional keyword arguments passed to scipy.integrate.quad function.

        Returns
        -------
        p: np.ndarray
            Legendre moments of the scatteing phase function
            from 0 to 'last' as a vector a numpy vector
            ([g_0, g_1, ..., g_last]).
        '''
        G = [self.g(g, **kwargs) for g in range(last + 1)]
        return np.asarray(G)

    def fastg(self, n: int, npts: int = 1000, pf: np.ndarray = None,
              costheta: np.ndarray = None) -> float:
        '''
        Computes the n-th Legendre moment of the scattering phase function.
        A fast fixed-step simpson quadrature is used to compute the Legendre
        moment. Note that the default number of points "npts"
        at which the value of scattering phase function is computed might be
        insufficient for accurate estimation of higher order Legendre moments.

        Parameters
        ----------
        n: int
            The Legendre moment to calculate. (:math:`n=1` yields the
            anisotropy factor  :math:`g_1`)
        npts: int
            Number of points (scattering angle cosines) at which the scattering
            phase function is estimated. The value is ignored if the value of
            costheta parameter is not None.
        pf: np.ndarray
            A vector containing the values of the scattering phase function at
            equally spaced scattering angle cosines from [-1, 1]. If None,
            the values are calculated at scattering angle cosines defined by
            argument costheta.
        costheta: np.ndarray
            A vector of equally spaced numbers (scattering angle cosines) on
            the interval [-1, 1] at which the scattering phase function
            is defined. If None, the values are calculated as:

            - :code:`np.linspace(-1.0, 1.0, npts)` if pf is None
            - :code:`np.linspace(-1.0, 1.0, pf.size)` if pf is not None

        Returns
        -------
        p: float
            Value of the n-th Legendre moment.
        '''
        if costheta is None:
            if pf is None:
                costheta = np.linspace(-1.0, 1.0, npts)
            else:
                costheta = np.linspace(-1.0, 1.0, pf.size)

        if pf is None:
            pf = self(costheta)

        return simps(pf*np.polynomial.legendre.Legendre.basis(n)(costheta),
                     dx=costheta[1] - costheta[0])

    def fastgs(self, last: int, npts: int = 1000, pf: np.ndarray = None,
               costheta: np.ndarray =None) -> np.ndarray:
        '''
        Computes the first n Legendre moments of the scattering phase function.
        A fast fixed-step simpson quadrature is used to compute the Legendre
        moments. Note that the default number of points "npts" at which the
        value of scattering phase function is computed might be insufficient
        for accurate estimation of higher order Legendre moments.

        Parameters
        ----------
        last: int
            The last Legendre moment to be calculate.
        npts: int
            Number of points (scattering angle cosines) at which the scattering
            phase function is estimated. The value is ignored if the value of
            costheta parameter is not None.
        pf: np.ndarray
            A vector containing the values of the scattering phase function at
            equally spaced scattering angle cosines from [-1, 1]. If None,
            the values are calculated at scattering angle cosines defined
            by argument costheta.
        costheta: np.ndarray
            A vector of equally spaced numbers (scattering angle cosines) on
            the interval [-1, 1] at which the scattering phase function is
            defined. If None, the values are calculated as:

            - :code:`np.linspace(-1.0, 1.0, npts)` if pf is None
            - :code:`np.linspace(-1.0, 1.0, pf.size)` if pf is not None

        Returns
        -------
        p: np.ndarray
            Values of Legendre moments from 0 to 'last' as a vector ([g_0, g_1, ..., g_last]).
        '''
        if costheta is None:
            costheta = np.linspace(-1.0, 1.0, npts)
        if pf is None:
            pf = self(costheta)

        G = np.zeros([last + 1])
        for n in range(last + 1):
            G[n] = simps(pf*np.polynomial.legendre.Legendre.basis(n)(costheta),
                         dx=costheta[1] - costheta[0])
        return G

    def cdf(self, costheta: np.ndarray, meth: str = 'simpson',
            npts: int = 10000, **kwargs) -> np.ndarray:
        '''
        Cumulative probability density function calculated at the specified
        deflection angle cosines.

        Parameters
        ---------
        costheta: np.ndarray
            A vector of scattering angle cosines at which the cumulative
            distribution function is to be computed.
        meth: str
            Numerical integration method:
                - 'trapez'  - trapezoidal
                - 'simpson' - default
                - 'quad'    - adaptive step; accurate but slower

        npts: int
            Number of control points used by the trapezoidal and simpson
            numerical integration methods.
        kwargs: dict
            Optional keyword arguments passed to the scipy.integrate.quad
            function.

        Returns
        -------
        pfcum: np.ndarray
            Returns the cumulative probability density values at the given
            scattering angle cosines.

        Note
        ----
        Computation method 'quad' might be slower on some platforms, however
        produces accurate results. Use 'quad' whenever possible. Otherwise
        check the accuracy of results obtained by 'simpson' or 'trapez' and
        adjust parameter npts if required (increase to improve accuracy,
        decrease for faster computation and increased computational error).

        Examples
        --------

        >>> from matplotlib import pyplot as pp
        >>> import numpy as np
        >>>
        >>> costheta = np.linspace(-1.0, 1.0, 1000)
        >>> # define a phase function
        >>> miepf = Mie(1.6, 1.33, 4.5e-6, 550e-9)
        >>> # accurate but sometimes slow adaptive step-size integration by quad
        >>> cdpf = miepf.cdf(costheta, 'quad')
        >>> # fast computation using simpson integration with 10000 control points
        >>> cpffast = miepf.cdf(costheta, 'simpson', 10000)
        >>>
        >>> pp.semilogy(costheta, cdpf)
        >>> pp.semilogy(costheta, cpffast)
        '''
        costheta = np.asarray(costheta, dtype=np.float64)
        pfcum = np.zeros([costheta.size])
        n = len(costheta)

        if meth == 'trapez':
            # trapez
            npts = max(npts, costheta.size)
            costhetai = np.linspace(-1, 1, npts)
            p = self(costhetai)
            pfcumi = np.zeros(int(p.size))
            dx = (2.0/(npts - 1.0))
            pfcumi[1:] = 0.5*(p[:-1] + p[1:])*dx
            pfcumi = pfcumi.cumsum()
            pfcum = np.interp(costheta, costhetai, pfcumi)
            pfcum1 = pfcumi[-1]
        elif meth == 'simpson':
            # simpson
            npts = int(max(npts, costheta.size)/2) * 2 + 1
            costhetai = np.linspace(-1, 1, npts)
            p = self(costhetai)
            pfcumi = np.zeros(p.size//2 + 1)
            dx = 2.0/(npts - 1.0)
            pfcumi[1:] = dx/3.0*(p[:-2:2] + 4.0*p[1:-1:2] + p[2::2])
            pfcumi = pfcumi.cumsum()
            pfcum = np.interp(costheta, costhetai[::2], pfcumi)
            pfcum1 = pfcumi[-1]
        else:
            # directly call the probability density function
            for i in range(n - 1):
                pfcum[i + 1] = pfcum[i] + \
                    quad(self, costheta[i], costheta[i + 1])[0]
            pfcum1 = 1.0
        # print(pfcum[0], pfcum1)
        pfcum *= 1.0/pfcum1
        return pfcum

    def mclut(self, n: int = 1000, ncd: int = 10000, **kwargs) \
            -> Tuple[Tuple[float, float, float], np.ndarray]:
        '''
        Prepares a nonlinear Monte Carlo lookup table for sampling the
        scattering angles of the scattering phase function.

        Parameters
        ----------
        n: int
            The LUT size.
        ncd: int
            The number of equally spaced points (scattering angles) at which
            the cumulative probability density of the scattering phase function
            is calculated.
        kwargs: dict
            Parameters passed to the :py:meth:`PfBase.cdf` method of the
            scattering phase function.

        Returns
        -------
        params: Tuple[float, float, float]
            List or tuple of three parameters [a, b, c] that can be used
            to estimate the scattering angle cosine from a uniform random number
            F from [0, 1]:

            .. code-block:: python

                index = 0.5*(params[0]/(F - params[2]) - params[1] + 1.0)*(lut.size - 1)
                first = np.clip(int(np.floor(index)), 0, lut.size - 1)
                second = np.clip(first + 1, 0, lut.size - 1)
                dx = index - first
                costheta = lut[first]*(1.0 - dx) + lut[second]*dx

        lut: np.ndarray
            LUT data as a vector of n elements.

        Note
        ----
        The returned values can be used to create a Monte Carlo lookup table
        by calling constructor :code:`mcbase.pf.Lut(*pf.mclut())`.

        Examples
        --------
        Prepares a Monte Carlo lookup table-based scattering phase function
        of a water dispersion containing 1 um microspherical particles with a
        refractive index 1.6 (medium/water refractive index is 1.33) at 550 nm.

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>> from xopto.mcbase import mcpf
        >>> from xopto import pf
        >>>
        >>> mie_pf = pf.Mie(1.6, 1.33, 1e-6, 550e-9)
        >>> mcpf = mcpf.Lut(*mie_pf.mclut())
        >>>
        '''
        g = self.g(1)
        if g > 0.99:
            print('Warning: Lookup table-based approximation of the scattering '\
                  'phase functions with anisotropy factor above 0.99 {:.3f} '\
                  'can become inaccurate!'.format(g))
        # Use at least 10x the size of LUT for CDF approximation.
        if ncd is not None:
            ncd = max(10*int(n), ncd)
        # Increase the number of points for scattering phase functions
        # with a high g.
        if g > 0.95:
            ncd = max(50000, ncd)

        # uniform grid on cos(theta)
        costheta = np.linspace(-1, 1, n)
        # Cumulative probability density of the scattering phase function
        # pdf(cos(theta))
        costhetai = np.linspace(-1, 1, ncd)
        cumpfi = self.cdf(costhetai, **kwargs)

        # Find optimal fit of the LUT function/model to the scattering phase
        # function CDF.
        res1 = minimize(
            lambda x: np.linalg.norm(lut_function(costhetai, x) - cumpfi),
            [2], bounds=((1 + np.finfo(np.float64).eps, np.inf),))

        res2 = minimize(
            lambda x: np.linalg.norm(lut_function(costhetai, x) - cumpfi),
            [-1], bounds=((-np.inf, -np.finfo(np.float64).eps),))

        if res1.fun < res2.fun:
            res = res1
        else:
            res = res2

        # Transform cos(theta) using the optimal fit of the LUT function.
        randlut = lut_function(costheta, res.x)
        randlut[randlut < 0.0] = 0.0
        randlut[randlut > 1.0] = 1.0
        # Compute cos(theta) of the transformed values - inverse of the
        # scattering phase function CDF.
        lut = interp1d(cumpfi, costhetai, bounds_error=False, fill_value=0.0,
                       kind='linear')(randlut)

        c = res.x[0]
        # Parameters passed to the pf.Lut scattering phase function must be
        # organized as [a, b, c], where a = 2*c*(1 - c) and b = 2*c - 1.
        # See lut_function function for more details.
        mcparams = [2*c*(1 - c), 2*c - 1, c]

        return mcparams, lut

    def mcluterr(self, params: Tuple[float, float, float], lut: np.ndarray,
                npts: int = 1000) -> np.ndarray:
        '''
        Evaluate the performance of the lookup table created by the
        :py:meth:`PfBase.mclut` method.

        The parameters of this function are as returned by the
        :py:meth:`PfBase.mclut` method)

        Parameters
        ----------
        lut: np.ndarray
            Lookup table data as a numpy vector of length n.
        params: Tuple[float, float, float]
            A tuple of three parameters [a, b, c] that can be used to estimate
            the scattering angle cosine from a unifor random number
            F from [0, 1]:

            .. code-block: python

                index = 0.5*(params[0]/(F - params[2]) - params[1] + 1.0)*(lut.size - 1)
                first = np.clip(int(np.floor(index)), 0, lut.size - 1)
                second = np.clip(first + 1, 0, lut.size - 1)
                dx = index - first
                costheta = lut[first]*(1.0 - dx) + lut[second]*dx

        npts: int
            Number of equally spaced points from [0.0, 1.0] at which the
            scattering phase function is evaluated.

        Returns
        -------
        costheta: np.ndarray
            The estimated scattering angle cosines.
        flut: np.ndarray
            Cumulative probability density of the scattering phase function at
            the given scattering angle cosines computed from the lookup table.
        ftrue: np.ndarray vector
            True values of the cumulative probability density of the scattering
            phase function at the given scattering angle cosines.

        Examples
        --------

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> pf = Mie(1.6, 1.33, 5e-6, 550e-9)
        >>> costheta, flut, ftrue = pf.mcluterr(*pf.mclut(512))
        >>>
        >>> pp.figure()
        >>> pp.subplot(121)
        >>> pp.semilogy(costheta[1:], flut[1:], label='CDF Estimated by LUT')
        >>> pp.semilogy(costheta[1:], ftrue[1:], label='True CDF')
        >>> pp.grid()
        >>> pp.legend()
        >>>
        >>> pp.subplot(122)
        >>> pp.plot(costheta[1:], (flut[1:] - ftrue[1:])/ftrue[1:])
        >>> pp.grid()
        >>> pp.title('Relative CDF error')
        '''
        randnum = np.linspace(0.0, 1.0, npts)

        index = 0.5*(params[0]/(randnum - params[2]) - params[1] + 1.0)*(lut.size - 1)
        first = np.clip((np.floor(index)).astype(np.int), 0, lut.size - 1)
        second = np.clip(first + 1, 0, lut.size - 1)
        dx = index - first

        costhetalut = lut[first]*(1.0 - dx) + lut[second]*dx
        true = self.cdf(costhetalut)

        return costhetalut, randnum, true

    def lut(self, n: int = 1000, ncd: int = 10000, **kwargs):
        # uniform grid on cos(theta)
        costheta = np.linspace(-1, 1, n)
        # cum. phase function pdf(cos(theta))
        costhetai = np.linspace(-1, 1, ncd)
        cumpfi = self.cdf(costhetai, **kwargs)

        # find optimal fit of the cumulative PDF by the lUT function
        res = minimize(
            lambda x: np.linalg.norm(lut_function(costhetai, x) - cumpfi), [2],
            bounds=((1, None),))
        # transform cos(theta) by the optimal LUT function
        randlut = lut_function(costheta, res.x)
        randlut[randlut < 0.0] = 0
        randlut[randlut > 1.0] = 1.0
        # compute cos(theta) of the transformed values - inverse phase function cum. PDF
        lut = interp1d(cumpfi, costhetai, bounds_error=False, fill_value=0.0,
                       kind='linear')(randlut)

        return lut, res.x, lambda values: \
            np.interp(ilut_function(values, res.x), np.linspace(-1, 1, lut.size), lut)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))
