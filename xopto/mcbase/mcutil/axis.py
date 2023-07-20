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


class Axis:
    '''
    An accumulator/detector axis with linearly or logarithmically spaced
    points/bins.
    '''
    def __init__(self, start: float or 'Axis' = 0.0, stop: float = 1.0, 
                 n: int = 1, logscale: bool=False):
        '''
        Creates a linearly or logarithmically spaced accumulator/detector axis.
        The parameters and data points of the axis can be accessed through
        several properties.

        Parameters
        ----------
        start: float or Axis
            Start coordinate (left edge) of the accumulator. If an instance of
            Axis, an identical copy of the instance is created.
        stop: float
            Stop coordinate (right edge) of the accumulator.
        n: int
            Number of bins in the accumulator.
        logscale: bool
            Scale of the accumulator. Linear if False (default) or logarithmic
            if True. Note that values of the start and stop parameters are
            threated as unscaled even if logscale is True.

        Note
        ----
        Note that only the r-axis of RT accumulators supports log-scale.

        Detailed description of the class properties:

        - start: float
            Left edge of the first bin.

        - scaled_start: float
            Logarithm of the start point if bins are distributed
            logarithmically, else equals start.

        - stop: float
            Right edge of the final bin.

        - scaled_stop: float
            Logarithm of the stop point if bins are distributed
            logarithmically, else equals stop.

        - span: tuple of two float
            Full range of the accumulator axis as a tuple (start, stop).

        - scaled_span: tuple of two float
            Full scaled range of the accumulator axis as a tuple
            (scaled_start, scaled_stop).

        - n: int
            Number of bins along the axis.

        - logscale: bool
            True, if bins are distributed logarithmically.

        - edges: np.ndarray vector
            Numpy vector of bin edges (number of points is n + 1).

        - centers: np.ndarray vector
            Central positions of the bins.
        '''
        if isinstance(start, Axis):
            axis = start
            start = axis.start
            stop = axis.stop
            n = axis.n
            logscale = axis.logscale
        else:
            start = float(start)
            stop = float(stop)
            n = int(n)
            logscale = bool(logscale)

        self._logscale = logscale
        self._n = n

        if self._logscale:
            start = max(start, float(np.finfo(np.float64).eps))
            self._edges = np.logspace(np.log(start), np.log(stop), n + 1,
                                      base=np.e)
            self._step = (np.log(stop) - np.log(start))/self._n
        else:
            self._edges = np.linspace(start, stop, n + 1)
            self._step = (stop - start)/self._n
        self._span = np.array((start, stop), dtype=np.float64)

        if self._logscale:
            self._scaled_span = np.log(self._span)
        else:
            self._scaled_span = self._span

        self._centers = 0.5*(self._edges[:-1] + self._edges[1:])

    def todict(self):
        '''
        Export object toa dict.
        '''
        return {'start':self._span[0], 'stop':self._span[1], 'n':self._n,
                'logscale':self._logscale, 'type':'Axis'}

    @classmethod
    def fromdict(cls, data: dict) -> 'Axis':
        '''
        Create a new object from dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        data_ = dict(data)
        type_name = data_.pop('type')
        if type_name != cls.__name__:
            raise TypeError('Expected "{}" type bot got "{}"!'.format(
                cls.__name__, type_name))
        return cls(**data_)

    def _get_logscale(self) -> bool:
        return self._logscale
    logscale = property(_get_logscale, None, None, 'Accumulator scale.')

    def _get_step(self) -> float:
        return self._step
    step = property(_get_step, None, None, 'Accumulator step.')

    def _get_n(self) -> int:
        return self._n
    n = property(_get_n, None, None, 'Number of accumulators along the axis.')

    def _get_span(self):
        return self._span
    span = property(_get_span, None, None, 'Accumulator span as [min, max].')

    def _get_start(self) -> float:
        return self._span[0]
    start = property(_get_start, None, None, 'Accumulator start coordinate.')

    def _get_stop(self) -> float:
        return self._span[1]
    stop = property(_get_stop, None, None, 'Accumulator stop coordinate.')

    def _get_scaled_span(self) -> Tuple[float, float]:
        return self._scaled_span
    scaled_span = property(_get_scaled_span, None, None,
                           'Scaled accumulator span as [min, max].')

    def _get_scaled_start(self) -> float:
        return self._scaled_span[0]
    scaled_start = property(_get_scaled_start, None, None,
                            'Scaled accumulator start coordinate.')

    def _get_scaled_stop(self) -> float:
        return self._scaled_span[1]
    scaled_stop = property(_get_scaled_stop, None, None,
                           'Scaled accumulator stop coordinate.')

    def _get_edges(self) -> np.ndarray:
        return self._edges
    edges = property(_get_edges, None, None, 'Edges of the accumulators.')

    def _get_centers(self) -> np.ndarray:
        return self._centers
    centers = property(_get_centers, None, None,
                       'Center points of the accumulators.')

    def __str__(self):
        return 'Axis(start={}, stop={}, n={}, logscale={})'.format(
            self._span[0], self._span[1], self._n, self._logscale)

    def __repr__(self):
        return self.__str__() + \
            ' # object at 0x{:>08X}.'.format(id(self))


class RadialAxis(Axis):
    '''
    A radial accumulator/detector axis with linearly or logarithmically spaced
    points/bins and corrected accumulator centers.
    '''
    def __init__(self, start: float=0.0, stop: float=1.0, n: int=1,
                 logscale: bool=False):
        '''
        Radial axis with corrected bin centers.

        The central point of each bin is computed as:
            :math:`2/3(r_{i}^2 + r_{i}*r_{i + 1} + r_{i + 1}^2)/(r_{i} + r_{i + 1})`

        For a detailed description of parameters see the Axis documentation.

        Note
        ----
        Correction is based on the assumption that for linear functions of
        reflectance :math:`f(r) = a + b r` the integral of
        :math:`f(r) 2 \\pi r dr` from :math:`r_i` to :math:`r_{i + 1}` is equal
        to the integral of :math:`f(r_center) 2 \pi r dr` from
        :math:`r_i` to :math:`r_{i + 1}`, where
        :math:`r_i < r_center < r_{i + 1}`.
        Computing the two integrals and equating the terms at coefficients
        :math:`a` and :math:`b` yields the "equivalent" center point of the bin
        :math:`[r_i, r_{i + 1})`.
        Note that the two integrals represent surface integral of reflectace.

        Warning
        -------
        The spacing of points in the RadialAxis is uneven even if logscale is
        set to False, since correction of the bin centers depends on the
        absolute position of the bin.
        If using data from RadialAxis with the simpson integration method make
        sure to use a version that does not require a fixed step (specify x
        instead of dx in a call to scipy.integrate.simps).
        If using the :py:meth:`xopto.util.hankel.discrete_simpson` method make
        sure to set the value of the uneven parameter to True.
        '''
        Axis.__init__(self, start, stop, n, logscale)
        edges = self.edges
        self._true_centers = (2.0/3.0)*(edges[:-1]**2 + edges[:-1]*edges[1:] +
                                        edges[1:]**2)/(edges[:-1] + edges[1:])

    def _get_centers(self) -> np.ndarray:
        return self._true_centers
    centers = property(_get_centers, None, None,
                       'Corrected centers of the accumulator bins.')

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        data = Axis.todict(self)
        data['type'] = 'RadialAxis'

    @classmethod
    def fromdict(cls, data:dict) -> 'RadialAxis':
        '''
        Create a new object from dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        return cls(**data)


class SymmetricAxis:
    '''
    An accumulator/detector axis with linearly or logarithmically spaced
    points/bins placed symmetrically around a center.
    '''

    def __init__(self, center: float=0.0, range: float=1.0, n_half: int=1000,
                 logscale: bool=False):
        '''
        Creates a linearly or logarithmically spaced symmetric accumulator axis.
        The parameters and data points of the axis can be accessed through
        several properties.

        Parameters
        ----------
        center: float
            Position of the axis/accumulator origin.
        range: float
            Extent of the axis as [center - range, center + range].
        n_half: int
            Number of bins in the interval [center, center + range]. The total
            number of bins equals 2*n_half.
        logscale: bool
            Scale of the accumulator. Linear if False (default) or logarithmic
            if True. Note that the values of the center and range parameters
            are threated as unscaled even if logscale is set to True.

        Note
        ----
        Note that only some RT classes support the SymmetricAxis class.

        Detailed description of the class properties:

        - center: float
            Coordinate of the center/origin of the symmetric accumulator.

        - range: float
            Range of the symmetric accumulator such that the full span of
            the symmetric accumulator is [center - range, center + range].

        - start: float
            Left edge of the first bin of the positive half interval. Note that
            for logarithmically spaced accumulators the start does not equal
            center but is instead slightly (eps) larger than the center.

        - stop: float
            Right edge of the final bin of the positive interval.

        - offset: float
            Offset of the left edge of the first bin along the
            positive interval. This value is 0 for linearly spaced
            accumulators but is a small nonzero value (eps) for
            logarithmically spaced values.

        - scaled_offset: float
            Logarithm of the offset of the left edge of the first bin along the
            positive interval. This value is 0 for linearly spaced
            accumulators but is a logarithm of a small nonzero value (ln(eps))
            for logarithmically spaced accumulators.

        - span: tuple of two float
            Full range of the accumulators (start, stop).

        - n_half:int
            Number of bins along the positive axis [center, center + range].

        - n: int
            Number of bins (2*n_half) along the full axis range
            [center - range, center + range].

        - logscale: bool
            True, if bins are distributed logarithmically.

        - edges: np.ndarray vector
            Numpy vector of bin edges (number of points is n + 1) in the
            full axis interval [center - range, center + range].

        - centers: np.ndarray vector
            Central positions of the bins in the
            full axis interval [center - range, center + range].
        '''
        if isinstance(center, SymmetricAxis):
            axis = center
            range = axis.range
            center = axis.center
            n_half = axis.n_half
            logscale = axis.logscale
        else:
            range = float(range)
            center = float(center)
            n_half = int(n_half)
            logscale = bool(logscale)

        self._range = range
        self._center = center
        self._logscale = logscale
        self._n_half = n_half
        self._n = 2*n_half

        if self._logscale:
            eps = np.finfo(np.float64).eps
            tmp = center + np.logspace(
                np.log(eps), np.log(range), n_half + 1, base=np.e)
            edges_right = tmp[1:]
            edges_left = -edges_right[::-1]

            self._offset = eps
            self._edges = np.hstack((edges_left, [center], edges_right))
            self._step = (np.log(range) - np.log(eps))/n_half

        else:
            self._offset = 0.0
            self._edges = np.linspace(
                center - range, center + range, 2*n_half + 1)
            self._step = (self._edges[-1] - self._edges[0])/self._n


        self._span = np.array(
            (self._center + self._offset, self._center + self._range),
            dtype=np.float64)
        if self._logscale:
            self._scaled_offset = np.log(self._offset)
        else:
            self._scaled_offset = self._offset

        self._centers = 0.5*(self._edges[:-1] + self._edges[1:])

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'range':self._range, 'center':self._center,
                'n_half':self._n_half, 'logscale':self._logscale,
                'type':'SymmetricAxis'}

    @classmethod
    def from_dict(cls, data: dict) -> 'SymmetricAxis':
        '''
        Create a new object from dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        return cls(**data)

    def _get_logscale(self) -> bool:
        return self._logscale
    logscale = property(_get_logscale, None, None,
                        'True if the accumulator is defined in log scale.')

    def _get_center(self) -> np.ndarray:
        return self._center
    center = property(_get_center, None, None, 'Axis/accumulator origin.')

    def _get_range(self) -> float:
        return self._range
    range = property(_get_range, None, None, 'Axis/accumulator range.')

    def _get_step(self) -> float:
        return self._step
    step = property(_get_step, None, None, 'Accumulator step/bin size.')

    def _get_n_half(self) -> int:
        return self._n_half
    n_half = property(_get_n_half, None, None,
                      'Total number of accumulator bins along '\
                      'the positive half axis.')

    def _get_n(self) -> int:
        return self._n
    n = property(_get_n, None, None,
                 'Total number of accumulator bins along the axis.')

    def _get_span(self) -> Tuple[float, float]:
        return self._span
    span = property(_get_span, None, None,
                    'Accumulator positive half axis span as [min, max].')

    def _get_start(self) -> float:
        return self._span[0]
    start = property(_get_start, None, None,
                     'Accumulator positive half interval start coordinate.')

    def _get_stop(self) -> float:
        return self._span[1]
    stop = property(_get_stop, None, None,
                    'Accumulator positive half interval stop coordinate.')

    def _get_scaled_offset(self) -> float:
        return self._scaled_offset
    scaled_offset = property(_get_scaled_offset, None, None,
                             'Scaled offset of the first bin ' \
                             'of the positive half interval '\
                             'from the center of the accumulator.')

    def _get_edges(self) -> np.ndarray:
        return self._edges
    edges = property(_get_edges, None, None,
                     'Edges of the accumulator bins.')

    def _get_centers(self) -> np.ndarray:
        return self._centers
    centers = property(_get_centers, None, None,
                       'Central points of the accumulator bins.')

    def __str__(self):
        return 'SymmetricAxis(range={}, center={}, n_half={}, '\
               'logscale={})'.format(
                   self._range, self._center, self._n_half, self._logscale)

    def __repr__(self):
        return self.__str__() + \
            ' # object at 0x{:>08X}.'.format(id(self))
