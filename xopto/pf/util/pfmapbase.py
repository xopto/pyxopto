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
import os.path
import pickle
import sys

import numpy as np
import scipy.io
from scipy.optimize import fmin, fmin_l_bfgs_b

from .helpers import sigma
from xopto import PICKLE_PROTOCOL, USER_DATA_PATH
from xopto.pf import PfBase


def load_map(filename: str) -> 'PfMap2DBase':
    '''
    Loads any of the scattering phase function maps from a file.
    For details see the documentation of MieFractalPolystyreneMap,
    MiePolystyreneMap, MHgMap, GkMap and MPcMap.
    '''
    ext = os.path.splitext(filename)[1]
    if ext == '.mat':
        data = scipy.io.loadmat(filename)
    else:
        with open(filename, 'rb') as fid:
            data = pickle.load(fid)
    type_name = data.get('type')
    if type_name is None:
        raise ValueError('File {} does not seem to contain a valid scattering '
                         'phase function parameter map! Field "type" is not '\
                         'defined!'.format(filename))
    if type_name in globals():
        return globals()[type_name].fromfile(filename)
    else:
        raise TypeError('A scattering phase function map type '\
                        '{} is not supported!'.format(type_name))


class PfMap2DBase:
    DEFAULT_MAP_FILE = None
    '''
    A filename for the default map.
    Overload this static member in derived classes.
    '''

    DEFAULT_MAP_PATH = os.path.join(USER_DATA_PATH, 'pf')
    '''
    A directory for the precalculated default map.
    Overload this static member in derived classes.
    '''

    XLABEL = '$first parameter$'
    YLABEL = '$second parameter$'
    PLOTSCALEFACTORX = 1.0
    PLOTSCALEFACTORY = 1.0

    '''
    Base class for computing maps of scattering phase function quantifiers
    as a function ot the scattering phase functions parameters - only for
    two-parametric scattering phase functions.
    '''

    @classmethod
    def default_data_file(cls) -> str:
        '''
        Creates a full filename for the default data file.

        Returns
        -------
        filename: str
            Full filename of the default data file.
        '''
        return os.path.join(cls.DEFAULT_MAP_PATH, cls.DEFAULT_MAP_FILE)

    @classmethod
    def fromfile(cls, filename: str = None) -> 'PfMap2DBase':
        '''
        Loads a scattering phase function map from a file.

        Parameters
        ----------
        filename: str
            Source file. If None load the map from the default file or
            create and save a new default map if a default file does not exist.
            Location of the default map file is defined by the
            :py:attr:`DEFAULT_MAP_PATH` and :py:attr:`DEFAULT_MAP_FILE`
            static attributes.
        '''
        if filename is None:
            filename = cls.default_data_file()

            if not os.path.isfile(filename):
                print(
                    '\n'
                    'Precalculating default scattering phase function map {} '
                    'on the first use.\n'
                    'This one-time process can take a couple of minutes to '
                    'complete.\n'
                    'The results will be saved to "{}" for future '
                    'use.'.format(
                        cls.__name__, filename
                    )
                )
                cls.precalculate(filename=filename)

        return cls(filename=filename)

    def __init__(self, param1: np.ndarray = None, param2: np.ndarray = None,
                 ng: int = None, pf: PfBase = None, filename: str = None,
                 ncostheta: int = None, **kwargs):
        '''
        Base class constructor.

        Parameters
        ----------
        param1: np.ndarray vector
            Sampling points of the first parameter of the scattering phase
            function.

        param2: np.ndarray vector
            Sampling points of the second parameter of the scattering phase
            function.

        pf: xopto.pf.PfBase
            Scattering phase function instance used with this parameter map.

        ncostheta: int
            Number of nodes used to compute the Legendre moments. Use a
            large number (> 1000) for accurate results. If None, adaptive
            step integration is used, which is accurate but can become slow.

        kwargs: dict
            Additional parameters passed to the constructor of the
            scattering phase function.
        '''
        self._kwargs = kwargs

        self._pf = pf

        if ncostheta is not None:
            ncostheta = int(ncostheta)
            if ncostheta < 1000:
                print('Warning: the number of nodes (ncostheta parameter) '\
                      'used for fast computation of Legendre moments by a '\
                      'fixed step Simpson quadrature is low! Use at least '\
                      '1000 nodes! If higher order Legendre moments are '\
                      'computed (parameter ng is > 3) use at least 5000 nodes!')
        self._ncostheta = ncostheta

        self._param1 = None
        self._param1_bounds = None

        self._param2 = None
        self._param2_bounds = None

        self._param1_grid = None
        self._param2_grid = None

        # map of Legendre moments
        self._g_map = None

        # parameters of derived scattering phase function parameters
        self._gamma_map = None
        self._delta_map = None
        self._sigma_map = None

        self._ng = None

        if filename is not None:
            self._load(filename)
            self._param1_bounds = (self._param1.min(), self._param1.max())
            self._param2_bounds = (self._param2.min(), self._param2.max())
        else:
            self._param1 = np.asarray(param1, dtype=np.float64)
            self._param1_bounds = (self._param1.min(), self._param1.max())

            self._param2 = np.asarray(param2, dtype=np.float64)
            self._param2_bounds = (self._param2.min(), self._param2.max())

            self._param1_grid, self._param2_grid = np.meshgrid(
                self._param1, self._param2, indexing='ij')
            self._ng = int(ng)

            self._make_maps()

    def _load(self, filename: str):
        filename = os.path.abspath(filename)
        data = np.load(filename)

        for field in ('type', 'param1', 'param2', 'grid1', 'grid2', 'g_map',
                      'gamma_map', 'delta_map', 'sigma_map', 'ng'):
            if field not in data.keys():
                raise ValueError(
                    'File \"{}\" is not a valid map! Missing field {}!'.format(
                        filename, field))

        if data.get('type') != self.__class__.__name__:
            raise ValueError(
                'Data in file \"{}\" bellong to a wrong class \"{}\" and '
                'cannot be used to initialize the target class {}!'.format(
                    filename, data['type'], self.__class__.__name__))

        self._param1 = data.get('param1')
        self._param2 = data.get('param2')
        self._param1_grid = data.get('grid1')
        self._param2_grid = data.get('grid2')
        self._g_map = data.get('g_map')
        self._delta_map = data.get('delta_map')
        self._gamma_map = data.get('gamma_map')
        self._sigma_map = data.get('sigma_map')
        self._ng = int(data.get('ng'))

    def save(self, filename: str):
        '''
        Save instance to a file.

        Parameters
        ----------
        filename: str
            The target file name.
        '''
        filename = os.path.abspath(filename)
        filename = os.path.splitext(filename)[0]

        target_dir = os.path.dirname(os.path.abspath(filename))
        try:
            os.makedirs(target_dir)
        except OSError:
            pass

        if not os.path.isdir(target_dir):
            raise OSError(
                'Failed to create a destination directory for the scattering '
                'phase function map!'
            )

        data = {'type':self.__class__.__name__,
                'param1': self._param1,
                'param2': self._param2,
                'grid1': self._param1_grid,
                'grid2': self._param2_grid,
                'g_map': self._g_map,
                'gamma_map': self._gamma_map,
                'delta_map': self._delta_map,
                'sigma_map': self._sigma_map,
                'ng': self._ng,
                'pftype': self.pf_type_name()}

        np.savez_compressed(filename + '.npz', **data)

    def _make_maps(self):
        if self._ncostheta is not None:
            print('Fast estimation of Legendre moments using {} '\
                  'nodes in the direction of cos(theta).'.format(
                      self._ncostheta))

        name = self.pf_type_name()
        self._g_map = np.zeros((self.n1(), self.n2(), self._ng + 1))

        for i in range(self.n1()):
            for j in range(self.n2()):

                sys.stdout.write(
                    '\rComputing parameter map for the ' \
                    '{} scattering phase function ... {:.1f}% '.format(
                        name,
                        (j + i*self.n2() + 1)/
                        self._param1_grid.size*100)
                )
                sys.stdout.flush()

                if self._ncostheta is None:
                    self._g_map[i, j, :] = self._pf(
                        self._param1[i], self._param2[j],
                        **self._kwargs).gs(self._ng)
                else:
                    self._g_map[i, j, :] = self._pf(
                        self._param1[i], self._param2[j],
                        **self._kwargs).fastgs(self._ng, npts=self._ncostheta)

        if self._ng >= 2:
            self._gamma_map = (1.0 - self._g_map[:, :, 2]) /\
                              (1.0 - self._g_map[:, :, 1])

        if self._ng >= 3:
            tmp1 = 1.0/(1.0 - self._g_map[:, :, 1])
            self._delta_map = (1.0 - self._g_map[:, :, 3])*tmp1

            i = np.arange(1, self._ng)
            k = (-0.5)**(i + 1 - 2)
            k.shape = [1, 1, k.size]
            tmp1.shape = (self._g_map.shape[0], self._g_map.shape[1], 1)
            self._sigma_map = (k*(1.0 - self._g_map[:, :, 2:])*tmp1).sum(-1)

    def pf_type_name(self) -> str:
        '''
        Type name of the scattering phase function that is used in this map.
        '''
        return self._pf.__name__

    def pf(self) -> PfBase:
        '''
        Type/class of the scattering phase function that is used in this map.
        '''
        return self._pf

    def param1(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of
        the first scattering phase function parameter.
        '''
        return self._param1

    def n1(self) -> int:
        '''
        Returns the number of points defining the grid of
        the first scattering phase function parameter.
        '''
        return self._param1.size

    def param2(self) -> np.ndarray:
        '''
        Returns a vector of points defining the grid of
        the second scattering phase function parameter.
        '''
        return self._param2

    def n2(self) -> int:
        '''
        Returns the number of points defining the grid of
        the second scattering phase function parameter.
        '''
        return self._param2.size

    def grid1(self) -> np.ndarray:
        '''
        Returns a 2D array (meshgrid) of points defining the grid of
        the first scattering phase function parameter.
        '''
        return self._param1_grid

    def grid2(self) -> np.ndarray:
        '''
        Returns a 2D array (meshgrid) of points defining the grid of
        the second scattering phase function parameter.
        '''
        return self._param2_grid

    def ng(self) -> int:
        '''
        Returns the number of Legendre moments.
        '''
        return self._ng

    def g(self, param1: float or np.ndarray = None,
          param2: float or np.ndarray = None, p: int = 1)\
                -> float or np.ndarray:
        '''
        Returns the value of the p-the Legendre moment for the scattering
        phase function using the specified values of the scattering phase
        function parameters.
        If one of the scattering phase function parameters is None,
        the full map is returned.

        Parameters
        ----------
        param1: float or np.ndarray
            Value(s) of the first phase function parameter.
        param2: float  or np.ndarray
            Value(s) of the second phase function parameter.
        p: int
            Requested Legendre moment (p=1 ... anisotropy factor g).

        Returns
        -------
        g_p: float or np.ndarray
            The value(s) of the requested Legendre moments.
        '''
        if param1 is None or param2 is None:
            return self._g_map[:, :, p]
        else:
            return self._pf(param1, param2).g(p)

    def gamma(self, param1: float or np.ndarray = None,
              param2: float or np.ndarray = None) -> float or np.ndarray:
        '''
        Returns the value of :math:`\\gamma=(1 - g_2)/(1 - g_1)` for the
        scattering phase function ussing the specified values of the scattering
        phase function parameters.
        If one of the scattering phase function parameters is None,
        the full map is returned.

        Parameters
        ----------
        param1: float or np.ndarray
            Value(s) of the first parameter of the scattering phase function.
        param2: float or np.ndarray
            Value(s) of the second parameter of the scattering phase function.

        Returns
        -------
        g_p: float or np.ndarray
            The value(s) of gamma.
        '''
        if param1 is None or param2 is None:
            return self._gamma_map
        else:
            gs = self._pf(param1, param2).gs(2)
            return (1.0 - gs[2])/(1.0 - gs[1])

    def delta(self, param1: float or np.ndarray = None,
              param2: float or np.ndarray = None) -> float or np.ndarray:
        '''
        Returns the values of delta=(1 - g_3)/(1 - g_1) for the scattering
        phase function using the specified values of the scattering phase
        function parameters.
        If one of the scattering phase function parameters is None,
        the full map is returned.

        Parameters
        ----------
        param1: float, ndarray
            Value(s) of the first parameter of the scattering phase function.
        param2: float, ndarray
            Value(s) of the second parameter of the scattering phase function.

        Returns
        -------
        g_p: float, ndarray
            The value(s) of delta.
        '''
        if param1 is None or param2 is None:
            return self._delta_map
        else:
            gs = self._pf(param1, param2).gs(3)
            return (1.0 - gs[3]) / (1.0 - gs[1])

    def gammadelta(self, param1: float or np.ndarray = None,
                   param2: float or np.ndarray = None) -> float or np.ndarray:
        '''
        Returns the value of gamma=(1 - g_2)/(1 - g_1) for the scattering phase
        function using the specified values of the scattering phase function
        parameters.
        If one of the scattering phase function parameters is None,
        the full map is returned.

        Parameters
        ----------
        param1: float, ndarray
            Value(s) of the first parameter of the scattering phase function.
        param2: float, ndarray
            Value(s) of the second parameter of the scattering phase function.

        Returns
        -------
        g_p: float, ndarray
            The value(s) of gammas.
        '''
        if param1 is None or param2 is None:
            return self._gamma_map, self._delta_map
        else:
            gs = self._pf(param1, param2).gs(3)
            return (1.0 - gs[2])/(1.0 - gs[1]), (1.0 - gs[3])/(1.0 - gs[1])

    def sigma(self, param1: float or np.ndarray = None,
              param2: float or np.ndarray = None, ng: int = 15) \
                  -> float or np.ndarray:
        '''
        Returns the value of sigma for the scattering phase
        function at the specified values of the scattering phase function
        parameters.
        If one of the scattering phase function parameters is None,
        the full map is returned.

        Parameters
        ----------
        param1: float, ndarray
            Value(s) of the GK phase function parameter.
        param2: float, ndarray
            Value(s) of the GK phase function parameter.
        ng: int
            The number of Legendre moments that is used for sigma calculation.

        Returns
        -------
        g_p: float, ndarray
            The value(s) of sigma.
        '''

        if ng > self._g_map.shape[2] - 1:
            print(
                'The specified number of Legendre moments ({}) exceeds '
                'the accumulator size. Only the first {} Legendre moments '
                'will be used for sigma calculation!'.format(
                    ng, self._g_map.shape[2] - 1)
            )
            ng = self._g_map.shape[2] - 1

        if param1 is None or param2 is None:
            return self._sigma_map
        else:
            gs = self._pf(param1, param2).gs(ng)
            return sigma(gs)

    def _invgammadelta_init(self, gamma: float, delta: float) -> Tuple[float, float]:
        '''
        Estimates the scattering phase function parameters from the
        given values of gamma and delta. The best fit from the
        precalculated lookup table is returned.

        Parameters
        ----------
        gamma: float
            Target value of parameter gamma.

        delta: float
            Target value of parameter delta.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''
        err = (self._gamma_map - gamma)**2 + (self._delta_map - delta)**2
        i, j = np.unravel_index(np.argmin(err), err.shape)
        return self._param1[i], self._param2[j]

    def _invgammadelta_interp(self, gamma: float, delta: float) -> float:
        '''
        Estimates the scattering phase function parameters alpha and gg
        from the given values of gamma and delta. The best fit from the
        precalculated lookup table is refined by interpolation/optimization.

        Parameters
        ----------
        gamma: float
            Target value of parameter gamma.
        delta: float
            Target value of parameter delta.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''
        x0 = self._invgammadelta_init(gamma, delta)

        def _kfungd(self, x, gamma, delta):
            return ((self.gamma(x[0], x[1]) - gamma)**2 + \
                    (self.delta(x[0], x[1]) - delta)**2)**0.5
        return fmin(lambda x: _kfungd(self, x, gamma, delta), x0, disp=False)

    def invgammadelta(self, gamma: float , delta: float, **kwargs) \
            -> Tuple[float, float]:
        '''
        Legendre moments, and from there gamma and delta, are calculated on
        the fly, using the scattering phase function model.
        Calculating the scattering phase function moments on the fly is more
        accurate then using a precalculated map.

        Parameters
        ----------
        gamma: float
            Target value of parameter gamma.
        delta: float
            Target value of parameter delta.
        kwargs: dict
            Keyword arguments passed to the fmin_l_bfgs_b optimization
            function.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''

        # x0 = self._invgammadelta_init(gamma, delta)
        x0 = self._invgammadelta_interp(gamma, delta)
        bounds = (self._param1_bounds, self._param2_bounds)

        def _kfunslowgg(x, gammaTarget, deltaTarget):
            #gammaEst = self.gamma(x[0], x[1])
            #deltaEst = self.delta(x[0], x[1])
            gammaEst, deltaEst = self.gammadelta(x[0], x[1])
            return ((deltaEst - deltaTarget)**2 +
                    (gammaEst - gammaTarget)**2)**0.5

        if 'bounds' not in kwargs:
            kwargs['bounds'] = bounds

        #return fmin_tnc(lambda x: _kfunslowgg(x, gamma, delta), x0,
        #                approx_grad=True, **kwargs)[0]
        return fmin_l_bfgs_b(lambda x: _kfunslowgg(x, gamma, delta), x0,
                             approx_grad=True, **kwargs)[0]

    def _invgamma_init(self, g: float, gamma: float) -> Tuple[float, float]:
        '''
        Estimates the scattering phase function parameters from the
        given values of gamma and g (the first Legendre moment).
        The best fit from the precalculated lookup table is returned.

        Parameters
        ----------
        g: float
            Target value of the first Legendre moment.
        gamma: float
            Target value of parameter gamma.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''
        err = (self._g_map[:, :, 1] - g)**2 + (self._gamma_map - gamma)**2
        i, j = np.unravel_index(np.argmin(err), err.shape)
        return self._param1[i], self._param2[j]

    def _invgamma_interp(self, g: float, gamma: float) -> Tuple[float, float]:
        '''
        Estimates the scattering phase function parameters from the
        given values of gamma and g (the first Legendre moment).
        The best fit from the precalculated lookup table is refined
        by interpolation/optimization.

        Parameters
        ----------
        g: float
            Target value of the first Legendre moment.
        gamma: float
            Target value of parameter gamma.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''
        x0 = self._invgamma_init(g, gamma)

        def _kfungg(self, x, g, gamma):
            return (self.g(x[0], x[1]) - g)**2 + \
                   (self.gamma(x[0], x[1]) - gamma)**2

        return fmin(lambda x: _kfungg(self, x, g, gamma), x0, disp=False)

    def invgamma(self, g: float, gamma: float) -> Tuple[float, float]:
        '''
        Legendre moments, and from there gamma, are calculated on
        the fly, using the scattering phase function model.
        Calculating the scattering phase function moments on the fly is
        more accurate then using a precalculated map.

        Parameters
        ----------
        g: float
            Target value of the first Legendre moment.
        gamma: float
            Target value of parameter gamma.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''

        # x0 = self._invgamma_init(g, gamma)
        x0 = self._invgamma_interp(g, gamma)
        bounds = (self._param1_bounds, self._param2_bounds)

        def _kfunslowgg(x, g1Target, gammaTarget):
            g1Est = self.g(x[0], x[1])
            gammaEst = self.gamma(x[0], x[1])
            return (g1Est - g1Target)**2 + (gammaEst - gammaTarget)**2

        return fmin_l_bfgs_b(lambda x: _kfunslowgg(x, g, gamma), x0,
                             approx_grad=True, bounds=bounds, disp=False)[0]

    def _invsigma_init(self, g: float, sigma: float) -> Tuple[float, float]:
        '''
        Estimates the scattering phase function parameters from the
        given values of sigma and g (the first Legendre moment).
        The best fit from the precalculated lookup table is returned.

        Parameters
        ----------
        g: float
            Target value of the first Legendre moment.
        sigma: float
            Target value of parameter sigma.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''
        err = (self._g_map[:, :, 1] - g)**2 + (self._sigma_map - sigma)**2
        i, j = np.unravel_index(np.argmin(err), err.shape)
        return self._param1[i], self._param2[j]

    def _invsigma_interp(self, g: float, sigma: float) -> Tuple[float, float]:
        '''
        Estimates the scattering phase function parameters from the
        given values of sigma and g (the first Legendre moment).
        The best fit from the precalculated lookup table is refined
        by interpolation/optimization.

        Parameters
        ----------
        g: float
            Target value of the first Legendre moment.
        sigma: float
            Target value of parameter sigma.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''
        x0 = self._invsigma_init(g, sigma)

        def _kfungg(self, x, g, sigma):
            return (self.g(x[0], x[1]) - g)**2 + \
                    (self.sigma(x[0], x[1]) - sigma)**2
        return fmin(lambda x: _kfungg(self, x, g, sigma), x0, disp=False)

    def invsigma(self, g: float, sigma: float) -> Tuple[float, float]:
        '''
        Legendre moments, and from there sigma, are calculated on
        the fly, using the scattering phase function model.
        Calculating the phase function moments on the fly is more accurate
        then using a precalculated map.

        Parameters
        ----------
        g: float
            Target value of the first Legendre moment.
        sigma: float
            Target value of parameter sigma.

        Returns
        -------
        param1, param2: float, float
            The estimated values of the scattering phase function parameters.
        '''
        x0 = self._invsigma_interp(g, sigma)
        bounds = (self._param1_bounds, self._param2_bounds)

        def _kfunslowgg(x, g1Target, sigmaTarget):
            g1Est = self.g(x[0], x[1])
            sigmaEst = self.sigma(x[0], x[1])
            return (g1Est - g1Target)**2 + (sigmaEst - sigmaTarget)**2

        return fmin_l_bfgs_b(lambda x: _kfunslowgg(x, g, sigma), x0,
                             approx_grad=True, bounds=bounds, disp=False)[0]

    def show(self):
        '''
        Show the scattering phase function maps using matplotlib.
        '''
        import matplotlib.pyplot as pp

        extent = np.array((*self._param1_bounds, *self._param2_bounds))
        extent[:2] *= self.PLOTSCALEFACTORX
        extent[2:] *= self.PLOTSCALEFACTORY

        imshow_options = {'cmap': 'hot',
                           'origin': 'lower',
                           'aspect': 'auto',
                           'interpolation': 'None',
                           'extent': extent}

        listMaps = [self.g(), self._gamma_map, self._delta_map, self._sigma_map]
        titles = ['$g_1$', '$\\gamma$', '$\\delta$', '$\\sigma$']

        fig, axs = pp.subplots(nrows=2, ncols=2)
        for i in range(axs.size):
            axs.flat[i].imshow(listMaps[i],
                               vmax=np.percentile(listMaps[i], 99),
                               **imshow_options)
            axs.flat[i].set_title(titles[i], fontsize=16)
            axs.flat[i].set_xlabel(self.XLABEL, fontsize=16)
            axs.flat[i].set_ylabel(self.YLABEL, fontsize=16)
        pp.suptitle(self.pf_type_name(), fontsize=20)

        # divide phase function space to four quadrants
        param1avg = np.mean(self._param1_bounds)
        param2avg = np.mean(self._param2_bounds)
        quad1 = np.where((self._param1_grid <= param1avg) &
                         (self._param2_grid <= param2avg))
        quad2 = np.where((self._param1_grid <= param1avg) &
                         (self._param2_grid > param2avg))
        quad3 = np.where((self._param1_grid > param1avg) &
                         (self._param2_grid <= param2avg))
        quad4 = np.where((self._param1_grid > param1avg) &
                         (self._param2_grid > param2avg))
        listQuads = [quad1, quad2, quad3, quad4]
        colors = ('r', 'b', 'g', 'k')

        pp.figure()
        ax1 = pp.subplot(121)
        ax2 = pp.subplot(122)
        for i in range(len(listQuads)):
            ax1.plot(self._param1_grid[listQuads[i]] * self.PLOTSCALEFACTORX,
                     self._param2_grid[listQuads[i]] * self.PLOTSCALEFACTORY,
                     '+', c=colors[i])
            ax2.plot(self._gamma_map[listQuads[i]],
                     self._delta_map[listQuads[i]],
                     '+', c=colors[i])
        ax1.set_xlabel(self.XLABEL, fontsize=16)
        ax1.set_ylabel(self.YLABEL, fontsize=16)
        ax1.set_title('parameter space', fontsize=16)
        ax1.set_xlim(imshow_options['extent'][:2])
        ax1.set_ylim(imshow_options['extent'][2:])
        ax2.set_xlabel('$\\gamma$', fontsize=16)
        ax2.set_ylabel('$\\delta$', fontsize=16)
        ax2.set_title('$\\gamma - \\delta$ space', fontsize=16)
        pp.suptitle(self.pf_type_name(), fontsize=20)
