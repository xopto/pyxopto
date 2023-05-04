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

from typing import Callable, Tuple

import numpy as np
from scipy.stats import norm
from scipy.integrate import quad


class Uniform(object):
    def __init__(self, lower: float, upper: float):
        '''
        Create a uniform distribution instance. The object is callable as
        obj(x) and returns the value of distribution at x.

        Parameters
        ----------
        lower: float
            Lower bound of the uniform distribution.
        upper: float
            Upper bound of the uniform distribution.

        Examples
        --------
        Creates an object representing a uniform distribution on [0, 2].

        >>> ud = Uniform(0.0, 2.0)
        >>> from matplotlib import pyplot as pp
        >>> x = np.linspace(-1.0, 3.0, 1000)
        >>> pp.plot(x, ud(x))
        >>>
        '''
        self._range = (0.0, 1.0)
        self._inv_range = 1.0
        self._set_range([lower, upper])

    def raw_moment(self, n: int) -> float:
        '''
        Computes the n-th raw moment of the distribution p(x) as:
            int(p(x)*x**n) on the interval [range[0], range[1]].

        Parameters
        ----------
        n: int
            The requested raw moment.

        Returns
        -------
        n-th raw moment of the distribution defined on the interval
        on the interval [range[0], range[1]].
        '''
        p = np.arange(0, n + 1)
        k = 1.0/(n + 1.0)
        #return quad(lambda x: self(x)*(x**n), *self._range)[0]
        return np.sum((self._range[0]**p)*(self._range[1]**p[::-1]))*k

    def _update(self):
        self._inv_delta_range = 1.0/(self._range[1] - self._range[0])

    def __call__(self, x: float):
        if np.isscalar(x):
            if x < self._range[0] or x > self._range[1]:
                f = 0
            else:
                f = 1.0*self._inv_delta_range
        else:
            x = np.asarray(x)
            f = np.tile(self._inv_delta_range, x.shape)
            f[x < self._range[0]] = 0
            f[x > self._range[1]] = 0

        return f

    def cdf(self, x: float) -> float:
        return np.clip((x - self._range[0])*self._inv_delta_range, 0.0, 1.0)

    def _get_lower(self) -> float:
        return self._range[0]
    def _set_lower(self, lower: float):
        lower = float(lower)
        self._set_range((lower, self._range[1]))
    lower = property(_get_lower, _set_lower, None,
                     'Lower bound of the uniform distribution range.')

    def _get_upper(self) -> float:
        return self._range[1]
    def _set_upper(self, upper: float):
        upper = float(upper)
        self._set_range((self._range[0], upper))
    upper = property(_get_upper, _set_upper, None,
                     'Upper bound of the uniform distribution range.')

    def _get_range(self) -> Tuple[float, float]:
        return self._range
    def _set_range(self, bounds: Tuple[float, float]):
        lower = float(bounds[0])
        upper = float(bounds[1])
        if lower >= upper:
            raise ValueError('Lower bound >= upper bound!')
        self._range = (lower, upper)
        self._update()
    range = property(_get_range, _set_range, None,
                     'The range of uniform distribution as (lower, upper).')

    def _get_parameters(self) -> dict:
        return {'lower':self._range[0], 'upper':self._range[1]}
    parameters = property(
        _get_parameters, None, None,
        'All Uniform distribution parameters as a dict object.')

    def __repr__(self):
        return 'Uniform(lower={}, upper={})'.format(*self._range)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))


class Normal(object):
    def __init__(self, mean: float, sigma: float, clip: float = 5):
        '''
        Create a Normal distribution instance. The object is callable as
        obj(x) and returns the value of distribution at x.

        Parameters
        ----------
        mean: float
            Normal distribution mean.
        sigma: float
            Normal distribution standard deviation.
        clip: float
            The range of Normal distribution is clipped to:
                [mean - sigma*clip, mean + sigma*clip]

        Examples
        --------
        Creates object representing normal distribution with mean 1 and
        standard deviation 0.1.

        >>> nd = Normal(1.0, 0.1, clip=5)
        >>> from matplotlib import pyplot as pp
        >>> x = np.linspace(0.5, 1.5, 1000)
        >>> pp.plot(x, nd(x))
        >>>
        '''
        self._sigma = float(sigma)
        self._mean = float(mean)
        self._clip = clip
        self._update()

    def raw_moment(self, n: int) -> float:
        '''
        Computes the n-th raw moment of the distribution p(x) as:
            int(p(x)*x**n) on the interval [range[0], range[1]].

        Parameters
        ----------
        n: int
            The requested raw moment.

        Returns
        -------
        n-th raw moment of the distribution defined on the interval
        on the interval [range[0], range[1]].
        '''
        return quad(lambda x: self(x)*(x**n), *self._range)[0]

    def _update(self):
        self._range = (self._mean - self._sigma*self._clip,
                       self._mean + self._sigma*self._clip,)
        self._cdfoffset = norm.cdf(-self._clip)
        self._k = 1.0/(1.0 - 2*self._cdfoffset)

    def __call__(self, x):
        if np.isscalar(x):
            if x < self._range[0] or x > self._range[1]:
                f = 0.0
            else:
                f = norm.pdf(x, self._mean, self._sigma)*self._k
        else:
            x = np.asarray(x)
            f = norm.pdf(x, self._mean, self._sigma)*self._k
            f[x < self._range[0]] = 0.0
            f[x > self._range[1]] = 0.0

        return f

    def cdf(self, x: float) -> float:
        return np.clip(self._k*norm.cdf(x, self._mean, self._sigma) -
                       self._cdfoffset, 0.0, 1.0)

    def _get_mean(self) -> float:
        return self._mean
    def _set_mean(self, mean: float):
        self._mean = np.float64(mean)
        self._update()
    mean = property(_get_mean, _set_mean, None,
                    'Mean value of the Normal distribution.')

    def _get_sigma(self) -> float:
        return self._sigma
    def _set_sigma(self, sigma: float):
        self._sigma = float(sigma)
        self._update()
    sigma = property(_get_sigma, _set_sigma, None,
                     'Standard deviation of the Normal distribution.')

    def _get_parameters(self) -> dict:
        return {'mean':self._mean, 'sigma':self._sigma, 'clip':self._clip}
    parameters = property(
        _get_parameters, None, None,
        'All Normal distribution parameters as a dict object.')

    def _get_range(self) -> Tuple[float, float]:
        return self._range
    range = property(_get_range, None, None, 'Distribution range.')

    def __repr__(self):
        return 'Normal(mean={}, sigma={}, clip={})'.format(
            self._mean, self._sigma, self._clip)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))


class Fractal(object):
    def __init__(self, alpha: float, range: Tuple[float, float] = (10e-9, 10e-6)):
        '''
        Fractal distribution with one parameter alpha:

            p(d) = A*(1/d)**alpha

        The value of parameter A is computed so as to normalize the integral
        of the density function on the specified range to 1:

            A = (range[0])**(alpha + 1)/(alpha*(1 - (range[0]/range[1])**(alpha + 1)))

        Parameters
        ----------
        alpha: float
            Parameter alpha of the fractal distribution.
        range: Tuple[float, float]
            The nonzero distribution range as a tuple (start, stop), e.g.
            use (10e-9, 10e-6) for fractal distribution of particles
            limited to the finite range from 10 nm to 10 um.

        Examples
        --------
        Creates object representing fractal distribution with alpha=2.4.

        >>> fd = Fractal(2.4)
        >>> from matplotlib import pyplot as pp
        >>> x = np.linspace(10e-9, 10e-6, 1000)
        >>> pp.semilogy(x, fd(x))
        >>>
        '''
        self._alpha = float(alpha)
        self._range = (float(range[0]), float(range[1]),)
        self._update()

    def _update(self):

        d1, d2 = self._range
        self._k = (1 - self._alpha)/(
            (d2**(1 - self._alpha) - d1**(1 - self._alpha)))
        self._cdfoffset = self._k/(1 - self._alpha)*d1**(1 - self._alpha)

    def raw_moment(self, n: int) -> float:
        '''
        Computes the n-th raw moment of the distribution p(x) as:
            int(p(x)*x**n) on the interval [range[0], range[1]].

        Parameters
        ----------
        n: int
            The requested raw moment.

        Returns
        -------
        n-th raw moment of the distribution defined on the interval
        on the interval [range[0], range[1]].
        '''
        d1, d2 = self._range
        return self._k/(1 - self._alpha + n)*(d2**(1 - self._alpha + n) -
                                              d1**(1 - self._alpha + n))

    def __call__(self, x):
        return self._k*(x)**(-self._alpha)

    def cdf(self, x: float) -> float:
        x = np.clip(x, *self._range)
        return self._k/(1 - self._alpha)*x**(1 - self._alpha) - self._cdfoffset

    def _get_alpha(self) -> float:
        return self._alpha
    def _set_alpha(self, alpha: float):
        self._alpha = float(alpha)
        self._update()
    alpha = property(
        _get_alpha, _set_alpha, None,
        'Parameter alpha of the fractal distribution (1/d)**alpha.')

    def _get_range(self) -> Tuple[float, float]:
        return self._range
    def _set_range(self, range: Tuple[float, float]):
        self._range = (float(range[0]), float(range[1]),)
        self._update()
    range = property(_get_range, _set_range, None,
                     'Numerical integration range.')

    def _get_parameters(self) -> dict:
        return {'alpha':self._alpha, 'range':self._range}
    parameters = property(
        _get_parameters, None, None,
        'All fractal distribution parameters as a dict object.')

    def __repr__(self):
        return 'Fractal(alpha={}, range={})'.format(
            self._alpha, self._range)

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))


class Mixture(object):
    def __init__(self, distributions: Tuple[Callable[[float], float]],
                 weights: Tuple[float] or np.ndarray):
        '''
        Create a mixture distribution instance. Distribution objects
        must be callable as obj(x) and return the value of distribution at x.

        Parameters
        ----------
        distributions: Tuple[Callable[[float], float]]
            Individual distributions of the mixture.
        weights: Tuple[float] or np.ndarray
            Weights of the individual distributions. The sum of weights should
            equal 1.

        Examples
        --------
        Creates object representing normal distribution with mean 1 and
        standard deviation 0.1.

        >>> nd = Normal(1.0, 0.1, clip=5)
        >>> from matplotlib import pyplot as pp
        >>> x = np.linspace(0.5, 1.5, 1000)
        >>> pp.plot(x, nd(x))
        >>>
        '''
        if not isinstance(distributions, (tuple, list)):
            distributions = (distributions,)

        if not isinstance(weights, (tuple, list, np.ndarray)):
            weights = (weights,)

        self._distributions = tuple(distributions)
        self._weights = tuple(np.asarray(weights).tolist())
   
    def  weight(self, index: int) -> float:
        '''
        Returns weight of the distribution at the specified index or slice.
        
        Parameters
        ----------
        index: int, slice
            Index or slice of the selected distribution weight.

        Returns
        -------
        weight: float or tuple
            The selected distribution weight(s).
        '''
        return self._weights[index]

    def  distribution(self, index: int) -> Callable[[float], float]:
        '''
        Returns distribution at the specified index or slice.
        
        Parameters
        ----------
        index: int, slice
            Index or slice of the selected distribution.

        Returns
        -------
        pd: Uniform, Normal, Fractal, ... or list of Uniform, Normal, Fractal, ...
            The selected distribution object(s).
        '''
        return self._distributions[index]

    def raw_moment(self, n: int) -> float:
        '''
        Computes the n-th raw moment of the distribution p(x) as:
            int(p(x)*x**n) on the interval [range[0], range[1]].

        Parameters
        ----------
        n: int
            The requested raw moment.

        Returns
        -------
        gn: float
            The n-th raw moment of the distribution defined on the interval
            on the interval [range[0], range[1]].
        '''
        m = 0.0
        for weight, distribution in zip(self._weights, self._distributions):
            m += weight*distribution.raw_moment(n)
        return m

    def __call__(self, x: float) -> float:
        f = None
        for weight, distribution in zip(self._weights, self._distributions):
            if f is None:
                f = weight*distribution(x)
            else:
                f += weight*distribution(x)
        return f

    def cdf(self, x: float) -> float:
        f = None
        for weight, distribution in zip(self._weights, self._distributions):
            if f is None:
                f = weight*distribution.cdf(x)
            else:
                f += weight*distribution.cdf(x)

        return np.clip(f, 0.0, 1.0)

    def _get_range(self) -> Tuple[float, float]:
        ranges = [pd.range for pd in self._distributions]
        return (np.min(ranges), np.max(ranges))   
    range = property(_get_range, None, None, 'Mixture range.')

    def __repr__(self):
        return 'Mixture(distributions={}, weights={})'.format(
            self._distributions, self._weights, id(self))

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))
