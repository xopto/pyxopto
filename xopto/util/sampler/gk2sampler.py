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

from xopto.util.sampler.base import UniformSampler, PfSampler
from xopto.mcbase import mcpf
from xopto.pf.util import g2gamma, g2delta
from xopto.pf import Gk2
from scipy.optimize import fmin_l_bfgs_b


class Gk2Sampler(PfSampler):
    def __init__(self,
                 gamma: Tuple[float, float] or None = None,
                 delta: Tuple[float, float] or None = None,
                 g: Tuple[float, float] or None = None):
        '''
        Uniformly samples the native parameters of the two GK scattering phase
        functions. Samples with gamma or delta value outside of the specified
        range are rejected.

        Parameters
        ----------
        gamma: Tuple[float, float] or None
            Range of gamma values that will be sampled. If None, the range
            is computed from delta values. If gamma and delta are None,
            all samples are accepted.
        delta: Tuple[float, float] or None:
            Range of delta values that will be sampled. If None, the range
            is computed from the gamma values. If gamma and delta are None,
            all samples are accepted.
        g: Tuple[float, float] or None:
            The anisotropy of the sample must be within the specified range
            or the sample is rejected. If None, samples are accepted regardless
            of their anisotropy.
        '''
        self._g = self._gamma = self._delta =  None
        if gamma is not None:
            self._gamma = (float(gamma[0]), float(gamma[1]))
        if delta is not None:
            self._delta = (float(delta[0]), float(delta[1]))
        if g is not None:
            self._g = (float(g[0]), float(g[1]))

        g1_sampler = UniformSampler(0.0, 0.98)
        a1_sampler = UniformSampler(-0.5, 10.0)

        g2_sampler = UniformSampler(-0.98, 0.0)
        a2_sampler = UniformSampler(-0.5, 10.0)

        b_sampler = UniformSampler(0.0, 1.0)

        super().__init__(mcpf.Gk2, (g1_sampler, a1_sampler,
                               g2_sampler, a2_sampler, b_sampler))

    def __call__(self, freeze=0):
        '''
        Sample the scattering phase function and return a new
        instance. Samples with gamma delta values outside of the domain
        specified in the constructor are rejected.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.
        '''
        sample = None
        while (True):
            sample = super().__call__(freeze=freeze)
            gs = sample.pf().gs(3)
            g = gs[1]
            gamma, delta = g2gamma(g, gs[2]), g2delta(g, gs[3])
            # check g range if required
            if self._g is not None:
                if self._g[0] > g or g > self._g[1]:
                    continue
            # check gamma range if required
            if self._gamma is not None:
                if self._gamma[0] > gamma or gamma > self._gamma[1]:
                    continue
            # check delta range if required
            if self._delta is not None:
                if self._delta[0] > delta or delta > self._delta[1]:
                    continue

            break

        return sample

    def todict(self) -> dict:
        '''
        Export instance to a dict.

        Returns
        -------
        data: dict
            Instance data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'gamma': self._gamma,
            'delta': self._delta,
            'g': self._g
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Gk2Sampler':
        '''
        Create a new instance from a dict.

        Parameters
        ----------
        data: dict
            Instance data that were exported to a dict.
        '''
        data = dict(data)
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Expected data for Gk2Sampler but got {}!'.format(T))

        return cls(**data)

    def __str__(self):
        return 'Gk2Sampler(gamma={}, delta={}, g={})'.format(
            self._gamma, self._delta, self._g)


def _gk2_inv_gamma_delta(gamma: float , delta: float, **kwargs) \
        -> Tuple[float, float]:
    '''
    Legendre moments, and from there gamma and delta, are calculated on
    the fly, using the scattering phase function model.

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

    x0 = [0.5, 0.5, -0.5, 0.5, 0.0]
    bounds = ((0.0, 0.99), (-0.5, 20.0), (-0.99, 0.0), (-0.5, 20), (0.0, 1.0))

    def _kfunslowgg(x, gammaTarget, deltaTarget):
        gs = Gk2(*x).gs(3)
        gammaEst, deltaEst = g2gamma(gs[1], gs[2]), g2delta(gs[1], gs[3])
        err =((deltaEst - deltaTarget)**2 +
                (gammaEst - gammaTarget)**2)**0.5
        return err

    if 'bounds' not in kwargs:
        kwargs['bounds'] = bounds

    x, err = fmin_l_bfgs_b(lambda x: _kfunslowgg(x, gamma, delta), x0,
                           approx_grad=True, **kwargs)[:2]

    return x, err


class Gk2GammaDeltaSampler(PfSampler):
    def __init__(self,
                 gamma: Tuple[float, float] or None = None,
                 delta: Tuple[float, float] or None = None,
                 g: Tuple[float, float] or None = None):
        '''
        Uniformly samples the native parameters of the two GK scattering phase
        functions. Samples with gamma or delta value outside of the specified
        range are rejected.

        Parameters
        ----------
        gamma: Tuple[float, float] or None
            Range of gamma values that will be sampled. If None, the range
            is computed from delta values. If gamma and delta are None,
            all samples are accepted.
        delta: Tuple[float, float] or None:
            Range of delta values that will be sampled. If None, the range
            is computed from the gamma values. If gamma and delta are None,
            all samples are accepted.
        g: Tuple[float, float] or None:
            The anisotropy of the sample must be within the specified range
            or the sample is rejected. If None, samples are accepted regardless
            of their anisotropy.
        '''
        self._g = self._gamma = self._delta =  None
        if gamma is not None:
            self._gamma = (float(gamma[0]), float(gamma[1]))
        if delta is not None:
            self._delta = (float(delta[0]), float(delta[1]))
        if g is not None:
            self._g = (float(g[0]), float(g[1]))

        gamma_sampler = UniformSampler(*gamma)
        delta_sampler = UniformSampler(*delta)

        super().__init__(mcpf.Gk2, (gamma_sampler, delta_sampler))

    def __call__(self, freeze=0):
        '''
        Sample the scattering phase function and return a new
        instance. Samples with gamma delta values outside of the domain
        specified in the constructor are rejected.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.
        '''
        sample = None
        while (True):
            gamma = self._pf_args[0](freeze=freeze)
            delta = self._pf_args[1](freeze=freeze)

            try:
                raw_args, err = _gk2_inv_gamma_delta(gamma, delta)
            except ValueError:
                continue

            if err < 1e-5:
                sample = self._pf_type(*raw_args)
                # check g range if required
                if self._g is not None:
                    g = sample.pf().g(1)
                    if self._g[0] > g or g > self._g[1]:
                        continue
                break
            #print(err)

        return sample


    def todict(self) -> dict:
        '''
        Export instance to a dict.

        Returns
        -------
        data: dict
            Instance data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'gamma': self._gamma,
            'delta': self._delta,
            'g': self._g
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Gk2Sampler':
        '''
        Create a new instance from a dict.

        Parameters
        ----------
        data: dict
            Instance data that were exported to a dict.
        '''
        data = dict(data)
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Expected data for Gk2GammaDeltaSampler but got {}!'.format(T))

        return cls(**data)

    def __str__(self):
        return 'Gk2GammaDeltaSampler(gamma={}, delta={}, g={})'.format(
            self._gamma, self._delta, self._g)


if __name__ == '__main__':
    import matplotlib.pyplot as pp
    import numpy as np

    N = 500
    for s in [Gk2Sampler(gamma=(1.5, 2.5), g=(0.0, 0.99)),
              Gk2GammaDeltaSampler(gamma=(1.5, 2.5), delta=(1.0, 5.0), g=(0, 0.99))]:
        delta, gamma = [], []
        for i in range(N):
            gs = s().pf().gs(3)
            gamma.append(g2gamma(gs[1], gs[2]))
            delta.append(g2delta(gs[1], gs[3]))
            print('Completed: {:5.1f}%'.format((i + 1)/N*100.0), end='\r')

        pp.figure()
        pp.plot(np.array(gamma), np.array(delta), '.')

    pp.show()