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
        Create a uniform size distribution:

        .. math::

            p(d) =
            \\begin{cases}
            \\frac{d}{d_{max} - d_{min}} & d_{min} \\geq d \\leq d_{max} \\\\
            0 & \\text{elsewhere}
            \\end{cases}

        The object is callable as `obj(x)` and returns the value
        of the distribution at x.

        Parameters
        ----------
        lower: float
            Lower bound of the uniform distribution
            (:math:`d_{min}` in :math:`p(d)`).
        upper: float
            Upper bound of the uniform distribution
            (:math:`d_{max}` in :math:`p(d)`).

        Note
        ----
        This class reimplements the equality/inequality operator by
        overloading the `__eq__` method.

        The overloaded `__hash__` method computes the instance hash from the
        parameters of the distribution. All distribution instances with
        exactly the same parameter values return the same hash value. 

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

    def __eq__(self, other: 'Uniform'):
        return self.range[0] == other.range[0] and \
               self.range[1] == other.range[1]

    def __hash__(self) -> int:
        return hash((float(self._range[0]), float(self._range[1])))

    def todict(self) -> dict:
        '''
        Export distribution object to a dict.

        Returns
        -------
        data: dict
            Distribution object exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'lower': self.lower, 'upper': self.upper
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Uniform':
        '''
        Create a new object from dict data.

        Parameters
        ----------
        data: dict
            Data as returned by the  :py:meth`Uniform.todict` method.

        Returns
        -------
        distribution: Uniform
            A new instance of :py:class:`Uniform` class that is initialized
            with data from the input dict.
        '''
        data = dict(data)
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Expected data for type "Uniform" but got "{}"!'.format(T))
        return cls(**data)

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
        Create a Normal distribution instance with mean :math:`\mu` and
        standard deviation :math:`\\sigma`:

        .. math::

            p(d) = \\frac{1}{\\sigma\\sqrt{2\\pi}}e^{-\\frac{1}{2}\\left(\\frac{d - \\mu}{\\sigma}\\right)^2}

        The created object is callable as obj(x) and returns the value
        of Normal distribution at x.

        Parameters
        ----------
        mean: float
            Normal distribution mean (:math:`\\mu` in :math:`p(d)`).
        sigma: float
            Normal distribution standard deviation
            (:math:`\\sigma` in :math:`p(d)`).
        clip: float
            The range of Normal distribution is clipped to
            :math:`[\\mu - \\sigma \\cdot \\text{clip}, \\mu + \\sigma \\cdot \\text{clip}]`.

        Note
        ----
        This class reimplements the equality/inequality operator by
        overloading the `__eq__` method.

        The overloaded `__hash__` method computes the instance hash from the
        parameters of the distribution (including the value of `clip`).
        All distribution instances with exactly the same parameter values
        return the same hash value. 

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

    def __eq__(self, other: 'Normal'):
        return self.sigma == other.sigma and self.mean == other.mean and \
               self.range[0] == other.range[0] and \
               self.range[1] == other.range[1]

    def __hash__(self) -> int:
        return hash((float(self.mean), float(self.sigma),
                     float(self._range[0]), float(self._range[1])))

    def todict(self) -> dict:
        '''
        Export distribution object to a dict.

        Returns
        -------
        data: dict
            Distribution object exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'mean': self._mean, 'sigma': self._sigma,
            'clip': self._clip
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Normal':
        '''
        Create a new object from dict data.

        Parameters
        ----------
        data: dict
            Data as returned by the  :py:meth`Normal.todict` method.

        Returns
        -------
        distribution: Normal
            A new instance of :py:class:`Normal` class that is initialized
            with data from the input dict.
        '''
        data = dict(data)
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Expected data for type "Normal" but got "{}"!'.format(T))
        return cls(**data)

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
        Fractal distribution of size :math:`d` has one parameter
        :math:`\\alpha`:

        .. math::

            p(d) =
            \\begin{cases}
            A \\left(\\frac{1}{d}\\right)^{\\alpha} & d_{min} \\geq d \\leq d_{max} \\\\
            0 & \\text{elsewhere}
            \\end{cases}

        The value of constant :math:`A` is computed so as to normalize the
        integral of :math:`p(d)` over :math:`[d_{min}, d_{max}]` to 1:

        .. math::

            A =\\frac{r_{0}^{\\alpha + 1}}
            {\\alpha \\left(1 - \\left(\\frac{r_0}{r_1}\\right)^{\\alpha + 1}\\right)}

        Parameters
        ----------
        alpha: float
            Parameter alpha of the fractal distribution
            (:math:`\\alpha` in :math:`p(d)`).
        range: Tuple[float, float]
            The nonzero distribution range as a tuple (start, stop), e.g.
            use (10e-9, 10e-6) for fractal distribution of particles
            limited to the finite range from 10 nm to 10 um.
            (:math:`\\left(d_{min}, d_{max})\\right)` in :math:`p(d)`).

        Note
        ----
        This class reimplements the equality/inequality operator by
        overloading the `__eq__` method.

        The overloaded `__hash__` method computes the instance hash from the
        parameters of the distribution. All distribution instances with
        exactly the same parameter values return the same hash value. 

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

    def __eq__(self, other: 'Fractal'):
        return self.alpha == other.alpha and \
               self.range[0] == other.range[0] and \
               self.range[1] == other.range[1]

    def __hash__(self) -> int:
        return hash((float(self.alpha),
                     float(self._range[0]), float(self._range[1])))

    def todict(self) -> dict:
        '''
        Export distribution object to a dict.

        Returns
        -------
        data: dict
            Distribution object exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'alpha': self._alpha, 'range': self._range
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Fractal':
        '''
        Create a new object from dict data.

        Parameters
        ----------
        data: dict
            Data as returned by the  :py:meth`Fractal.todict` method.

        Returns
        -------
        distribution: Fractal
            A new instance of :py:class:`Fractal` class that is initialized
            with data from the input dict.
        '''
        data = dict(data)
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Expected data for type "Fractal" but got "{}"!'.format(T))
        return cls(**data)

    def _update(self):
        '''
        Internal method that is used to update precalculated values.
        '''
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

        Note
        ----
        This class reimplements the equality/inequality operator by
        overloading the `__eq__` method.

        The overloaded `__hash__` method computes the instance hash from the
        parameters and weights of all the distributions.
        All distribution instances with exactly the same parameters and weights
        return the same hash value. 

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

        if len(distributions) != len(weights):
            raise ValueError('The number of distributions must be the same '
                             'as the number of weights!')

        self._distributions = tuple(distributions)
        self._weights = tuple(np.asarray(weights).tolist())

    def __eq__(self, other: 'Mixture'):
        if len(self) == len(other):
            for dist1, dist2 in zip(self.distributions, other.distributions):
                if dist1 != dist2:
                    return False

            for w1, w2 in zip(self.weights, other.weights):
                if w1 != w2:
                    return False

        return True

    def __hash__(self) -> int:
        return hash(
            tuple((p, w) for p, w in zip(self.distributions, self.weights))
        )

    def todict(self) -> dict:
        '''
        Export distribution object to a dict.

        Returns
        -------
        data: dict
            Distribution object exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'distributions': [item.todict() for item in self._distributions],
            'weights': [float(item) for item in self._weights],
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Mixture':
        '''
        Create a new object from dict data.

        Parameters
        ----------
        data: dict
            Data as returned by the  :py:meth`Mixture.todict` method.

        Returns
        -------
        distribution: Mixture
            A new instance of :py:class:`Mixture` class that is initialized
            with data from the input dict.
        '''
        data = dict(data)
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Expected data for type "Mixture" but got "{}"!'.format(T))
        distributions = []
        for item in data.pop('distributions'):
            Ts = item.pop('type')
            T = globals().get(Ts)
            if T is None:
                T = locals().get(Ts)
            if T is None:
                raise TypeError('Unsupported distribution "{}"!'.format(Ts))
            distributions.append(T(**item))

        return cls(distributions, **data)

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

    def _get_weights(self) -> Tuple[float]:
        return self._weights
    weights = property(_get_weights, None, None, 'Distribution weights.')

    def _get_distributions(self) -> Tuple[Callable[[float], float]]:
        return self._distributions
    distributions = property(_get_distributions, None, None,
                             'Mixture distributions')

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

    def __len__(self):
        return len(self._distributions)

    def __repr__(self):
        return 'Mixture(distributions={}, weights={})'.format(
            self._distributions, self._weights, id(self))

    def __str__(self):
        return '{} # object @{}'.format(self.__repr__(), id(self))
