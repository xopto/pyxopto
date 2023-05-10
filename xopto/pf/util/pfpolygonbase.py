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
import pickle
import os.path
import sys

import numpy as np
from shapely import geometry

from . import helpers
from xopto import USER_DATA_PATH, PICKLE_PROTOCOL


class GammaDeltaPolygonBase:
    # Default location of precalculated boundary/polygon points
    DEFAULT_POLYGON_PATH = os.path.join(USER_DATA_PATH, 'pf')

    # Default filename for the precalculated boundary/polygon points.
    # Overload this static member in derived classes.
    DEFAULT_POLYGON_FILE = None

    # Default number of boundary points along a boundary segment
    # Overload this static member in derived classes.
    NUM_POINTS = 10000

    @classmethod
    def default_data_file(cls):
        '''
        Creates a full filename for the default data file.

        Returns
        -------
        filename: str
            Full filename of the default data file.
        '''
        return os.path.join(cls.DEFAULT_POLYGON_PATH, cls.DEFAULT_POLYGON_FILE)

    @classmethod
    def fromfile(cls, filename: str = None) -> 'GammaDeltaPolygonBase':
        '''
        Loads a domain polygon of GK scattering phase function from a file.

        Parameters
        ----------
        filename: str
            File with the parameter map data. If None, attempts to load a
            default map file from data/pf folder.

        Returns
        -------
        polygon: GkPolygon
            A GkPolygon instance
        '''
        if filename is None:
            filename = cls.default_data_file()
            if not os.path.isfile(filename):
                print(
                    '\n'
                    'Precalculating default boundary polygon {} on the'
                    'first use.\n'
                    'This one-time process can take a couple of minutes to '
                    'complete.\n'
                    'The results will be saved to "{}" for future '
                    'use.'.format(
                        cls.__name__, filename
                    )
                )
                cls.precalculate(filename=filename)

        return cls(filename=filename)

    @classmethod
    def precalculate(cls, n: int = None, filename: str = None,
                     verbose: bool = True):
        '''
        Precalculate a scattering phase function domain boundary points/polygon
        and save the data to the given file.

        Parameters
        ----------
        n: int
            Number of steps along the scattering phase function parameters.
            If None, a default value from NUM_POINTS static member variable is
            used.
        filename: str
            Output file or None to save as the default boundary.
        verbose: bool
            Turn on verbose progress report.
        '''
        if verbose:
            print('\nCreating {:s} gamma-delta polygon:'.format(cls.__name__))

        if n is None:
            n = cls.NUM_POINTS

        if filename is None:
            filename = os.path.join(
                cls.DEFAULT_POLYGON_PATH, cls.DEFAULT_POLYGON_FILE)

        target_dir = os.path.dirname(os.path.abspath(filename))
        try:
            os.makedirs(target_dir)
        except OSError:
            pass

        if not os.path.isdir(target_dir):
            raise OSError(
                'Failed to create user data directory for the precalculated '
                'scattering phase function boundary!'
            )

        cls(n=n, verbose=True).save(filename)

    @staticmethod
    def _contains_one_point(validator, x: float, y: float) -> bool:
        '''
        A helper function for checkin the location of a single point.
        '''
        return validator(geometry.Point(x, y))

    def __init__(self, gamma: np.ndarray, delta: np.ndarray,
                 filename: str = None):
        '''
        Constructor of a domain boundary polygon.

        Parameters
        ----------
        gamma: np.ndarray vector
            A vector of gamma coordinates that lie on the boundary.
        delta: np.ndarray vector
            A corresponding vector of delta coordinates that lie on the boundary.
        filename: str
            Load boundary from a file.
        '''
        if filename is not None:
            data = self._load(filename)
            gamma, delta = data['gamma'], data['delta']

        self._data = {'gamma':gamma, 'delta':delta}

        self._gamma_delta_polygon = self._make_polygon(
            self._data['gamma'], self._data['delta'])

        self._gamma_range = (gamma.min(), gamma.max())
        self._delta_range = (delta.min(), delta.max())

        self._contains_method = np.vectorize(
            GammaDeltaPolygonBase._contains_one_point, otypes=(np.bool,))

    def _make_polygon(self, x: np.ndarray, y: np.ndarray) -> geometry.Polygon:
        return geometry.Polygon(np.vstack((x, y,)).T)

    def boundary(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns points along the boundary of the domain as two ndarray vectors
        gamma, delta.

        Returns
        -------
        gamma: np.ndarray vector
            Gamma coordinates of the points along the boundary of the domain.

        delta: np.ndarray vector
            Delta coordinates of the points along the boundary of the domain.
        '''
        return self._data['gamma'], self._data['delta']

    def delta(self, gamma: float) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns two nearest points on boundary of the domain for the
        given gamma value or None if the given gamma value is out of range.

        Parameters
        ----------
        gamma: float
            Gamma value at which to find two nearest points on the boundary
            of the domain.

        Returns
        -------
        gamma: np.ndarray vector of 2 floats
            Gamma coordinates of the nearest two points on the
            boundary of the domain.
        delta: np.ndarray vector of 2 floats
            Delta coordinates of the nearest two points on the
            boundary of the domain.
        '''
        ref = self._data['gamma']
        if self._gamma_range[0] <= gamma <= self._gamma_range[1]:
            distance = ref - gamma
            abs_distance = np.abs(distance)

            left, = np.nonzero(distance[:-1]*distance[1:] <= 0)
            right = left + 1

            distance_left = abs_distance[left]
            distance_right = abs_distance[right]

            left_mask = distance_left < distance_right
            right_mask = np.logical_not(left_mask)

            delta_1 = self._data['delta'][left][left_mask]
            delta_2 = self._data['delta'][right][right_mask]

            gamma_1 = self._data['gamma'][left][left_mask]
            gamma_2 = self._data['gamma'][right][right_mask]

            delta_candidates = np.hstack((delta_1, delta_2))
            gamma_candidates = np.hstack((gamma_1, gamma_2))
            n = min(2, delta_candidates.size)

            return gamma_candidates[:n], delta_candidates[:n]

        return None, None

    def gamma(self, delta: float) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns two nearest points on the domain boundary for the
        given delta value or None if the given delta value is out of range.

        Parameters
        ----------
        gamma: float
            Gamma value at which to find two nearest points on domain boundary.

        Returns
        -------
        gamma: np.ndarray vector of 2 floats
            Gamma coordinates of the nearest two points on the
            boundary of the domain.
        delta: np.ndarray vector of 2 floats
            Delta coordinates of the nearest two points on the
            boundary of the domain.
        '''
        ref = self._data['delta']
        if self._delta_range[0] <= delta <= self._delta_range[1]:
            distance = ref - delta
            abs_distance = np.abs(distance)

            left, = np.nonzero(distance[:-1]*distance[1:] <= 0)
            right = left + 1

            distance_left = abs_distance[left]
            distance_right = abs_distance[right]

            left_mask = distance_left < distance_right
            right_mask = np.logical_not(left_mask)

            delta_1 = self._data['delta'][left][left_mask]
            delta_2 = self._data['delta'][right][right_mask]

            gamma_1 = self._data['gamma'][left][left_mask]
            gamma_2 = self._data['gamma'][right][right_mask]

            delta_candidates = np.hstack((delta_1, delta_2))
            gamma_candidates = np.hstack((gamma_1, gamma_2))
            n = min(2, delta_candidates.size)

            return gamma_candidates[:n], delta_candidates[:n]

        return None, None

    def bounding_box(self,
                    gamma: Tuple[float, float] or None = None,
                    delta: Tuple[float, float] or None = None) -> \
                        Tuple[Tuple[float, float], Tuple[float, float]]:
        '''
        Compute valid bounding box in the gamma-delta plane.

        Parameters
        ----------
        gamma: Tuple[float, float] or None
            Range of gamma values as a tuple (low, high). If None, the
            range is derived from the delta parameter. If gamma and
            delta are None, the full range of values is included in
            in the bounding box.
        delta: Tuple[float, float] or None
            Range of delta values as a tuple (low, high). If gamma and
            delta are None, the full range of values is included in
            in the bounding box.

        Returns
        -------
        gamma_range: Tuple[float, float]
            Valid range of gamma values as a tuple (low, high).
        delta_range: Tuple[float, float]
            Valid range of delta values as a tuple (low, high).
        '''
        gamma_boundary, delta_boundary = self.boundary()
        gamma_min, gamma_max = gamma_boundary.min(), gamma_boundary.max()
        delta_min, delta_max =  delta_boundary.min(), delta_boundary.max()

        if gamma is None and delta is None:
            gamma = (gamma_min, gamma_max)
            delta = (delta_min, delta_max)
        else:
            if gamma is not None:
                if gamma[0] < gamma_min and gamma[1] < gamma_min:
                    raise ValueError('The given range of gamma values does not '
                                     'intersect the domain!')
                if gamma[0] > gamma_max and gamma[1] > gamma_max:
                    raise ValueError('The given range of gamma values does not '
                                    'intersect the domain!')

                gamma = np.clip(gamma, gamma_min, gamma_max)
            if delta is not None:
                if delta[0] < delta_min and delta[1] < delta_min:
                    raise ValueError('The given range of delta values does not '
                                     'intersect the domain!')
                if delta[0] > delta_max and delta[1] > delta_max:
                    raise ValueError('The given range of delta values does not '
                                     'intersect the domain!')
                delta = np.clip(delta, delta_min, delta_max)

        if delta is None:
            _, delta_1 = self.delta(gamma[0])
            _, delta_2 = self.delta(gamma[1])

            delta = (
                float(np.min([np.min(delta_1), np.min(delta_2)])),
                float(np.max([np.max(delta_1), np.max(delta_2)]))
            )
        else:
            delta = float(np.min(delta)), float(np.max(delta))

        if gamma is None:
            gamma_1, _ = self.gamma(delta[0])
            gamma_2, _ = self.gamma(delta[1])

            gamma = (
                float(np.min(np.min(gamma_1), np.min(gamma_2))),
                float(np.max(np.max(gamma_1), np.max(gamma_2)))
            )
        else:
            gamma = float(np.min(gamma)), float(np.max(gamma))

        return gamma, delta

    def contains(self, gamma: float or np.ndarray, delta: float or np.ndarray,
                 boundary: bool = True) -> bool or np.ndarray:
        '''
        Fast polygon-based validation of the scattering phase function
        gamma-delta domain.

        Parameters
        ----------
        gamma: float or np.ndarray
            Parameter gamma
        delta: float or np.ndarray
            Parameter delta
        boundary: bool
            If True the region includes the boundary, else not.

        Returns
        -------
        valid: bool, ndarray
            True if valid (gamma and delta lie within the domain boundary),
            else False.
        '''
        gamma = np.asarray(gamma)
        delta = np.asarray(delta)

        #return mask
        if boundary:
            return self._contains_method(
                self._gamma_delta_polygon.intersects, gamma, delta)
        else:
            return self._contains_method(
                self._gamma_delta_polygon.contains, gamma, delta)

    def _load(self, filename: str) -> dict:
        '''
        Load object data from a file (compressed numpy file).

        Parameters
        ----------
        filename: str
            Load data from this file.

        Returns
        -------
        data: dict
            Object data with entries "gamma" and "delta".
        '''
        filename = os.path.abspath(filename)
        data = np.load(filename)

        for field in ('type', 'gamma', 'delta'):
            if field not in data.keys():
                raise ValueError(
                    'File \"{}\" is not a valid polygon! '
                    'Missing field {}!'.format(filename, field))

        if data.get('type') != self.__class__.__name__:
            raise ValueError(
                'Data in file \"{}\" bellong to a wrong class \"{}\" and '
                'cannot be used to initialize the target class {}!'.format(
                    filename, data['type'], self.__class__.__name__))

        return {'gamma': data['gamma'], 'delta': data['delta']}

    def save(self, filename: str):
        '''
        Save polygon data to a file.
        Parameters
        ----------
        filename: str
            Name of the file. Extension is not required.
        '''
        filename = os.path.abspath(filename)
        filename = os.path.splitext(filename)[0]

        self._data['type'] = self.__class__.__name__
        np.savez_compressed(filename + '.npz', **self._data)

    def show(self, marker='-', step=1, pfname=''):
        '''
        Show the domain boundary using matlotlib.
        '''
        from matplotlib import pyplot as pp
        pp.plot(self._data['gamma'].flat[::step],
                self._data['delta'].flat[::step],
                marker, label=pfname + ' gamma-delta')
        pp.xlabel('Gamma')
        pp.ylabel('Delta')
        pp.legend()
