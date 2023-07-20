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

from typing import Tuple, List

from xopto.util.sampler.base import PfSampler, UniformSampler
from xopto.mcbase import mcpf
from xopto.pf.util import MPcMap, MPcPolygon, g2gamma, g2delta


class MPcSampler(PfSampler):
    def __init__(self,
                 gamma: Tuple[float, float] or None = None,
                 delta: Tuple[float, float] or None = None,
                 g: Tuple[float, float] or None = None):
        '''
        Uniformly samples the native parameters of the MPC scattering phase
        function and rejects all samples that have the gamma or delta
        value outside of the specified range.

        Parameters
        ----------
        gamma: Tuple[float, float] or None
            Range of gamma values that will be sampled. Samples with gamma
            value outside of the specified interval will be rejected.
            If None, all values of gamma are accepted.
        delta: Tuple[float, float] or None:
            Range of delta values that will be sampled. Samples with delta
            value outside of the specified interval will be rejected.
            If None, all values of delta are accepted.
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

        n_sampler = UniformSampler(0.0, 500.0)
        b_sampler = UniformSampler(0.0, 1.0)

        super().__init__(mcpf.MPc, (n_sampler, b_sampler))

    def __call__(self, freeze=0):
        '''
        Sample the scattering phase function and return a new
        instance. Samples with gamma or delta value outside of the domain
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
    def fromdict(cls, data: dict) -> 'MPcSampler':
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
                'Expected data for MPcSampler but got {}!'.format(T))

        return cls(**data)

    def __str__(self):
        return 'MPcSampler(gamma={}, delta={}, g={})'.format(
            self._gamma, self._delta, self._g)


class MPcGammaDeltaSampler(PfSampler):
    def __init__(self,
                 gamma: Tuple[float, float] or None = None,
                 delta: Tuple[float, float] or None = None,
                 g: Tuple[float, float] or None = None):
        '''
        Uniformly samples the MPC scattering phase function in
        the gamma-delta plane. All samples with gamma or delta value
        outside of the specified range are rejected. 

        Parameters
        ----------
        gamma: Tuple[float, float] or None
            Range of gamma values that will be sampled. If None, the range
            is computed from delta values. If both gamma and delta are None,
            the full range of the MPC scattering phase function is used.
        delta: Tuple[float, float] or None:
            Range of delta values that will be sampled. If None, the range
            is computed from the gamma values. If both gamma and delta are None,
            the full range of the MPC scattering phase function is used.
        g: Tuple[float, float] or None:
            The anisotropy of the sample must be within the specified range
            or the sample is rejected. If None, samples are accepted regardless
            of their anisotropy.
        '''
        self._map = MPcMap.fromfile()
        self._polygon = MPcPolygon.fromfile()

        if gamma is None and delta is None:
            gamma_boundary, delta_boundary = self._polygon.boundary()
            gamma = gamma_boundary.min(), gamma_boundary.max()
            delta = delta_boundary.min(), delta_boundary.max()
        else:
            gamma, delta = self._polygon.bounding_box(gamma, delta)
        if g is not None:
            g = (float(g[0]), float(g[1]))

        self._gamma, self._delta, self._g = gamma, delta, g

        gamma_sampler = UniformSampler(*gamma)
        delta_sampler = UniformSampler(*delta)

        super().__init__(mcpf.MPc, (gamma_sampler, delta_sampler))

    def __call__(self, freeze=0):
        '''
        Sample tha scattering phase function and return a new instance.
        Samples with gamma or delta value outside of the domain
        specified in the constructor are rejected.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.
        '''
        gamma = delta = None
        while True:
            gamma = self._pf_args[0](freeze=freeze)
            delta = self._pf_args[1](freeze=freeze)
            if not bool(self._polygon.contains(gamma=gamma, delta=delta)):
                continue

            raw_args = self._map.invgammadelta(gamma, delta)
            pf_sample = self._pf_type(*raw_args)
            gs = pf_sample.pf().gs(3)
            g = gs[1]
            # check g range if required
            if self._g is not None:
                if self._g[0] > g or g > self._g[1]:
                    continue
            # check fit
            gamma_, delta_ = g2gamma(g, gs[2]), g2delta(g, gs[3])
            err = ((gamma - gamma_)**2 + (delta - delta_)**2)
            if err < 1e-5:
                break

        return pf_sample


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
    def fromdict(cls, data: dict) -> 'MPcGammaDeltaSampler':
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
                'Expected data for MPcGammaDeltaSampler but got {}!'.format(T))

        return cls(**data)

    def __str__(self):
        return 'MPcGammaDeltaSampler(gamma={}, delta={}, g={})'.format(
            self._gamma, self._delta, self._g)


if __name__ == '__main__':
    import matplotlib.pyplot as pp

    for s in [MPcSampler(g=(0.0, 0.99)), MPcGammaDeltaSampler(g=(0.0, 0.99))]:
        delta, gamma = [], []
        for i in range(1000):
            gs = s().pf().gs(3)
            gamma.append(g2gamma(gs[1], gs[2]))
            delta.append(g2delta(gs[1], gs[3]))
            print('Completed: {:5.1f}%'.format((i + 1)/1000.0*100.0), end='\r')
        pp.figure()
        pp.plot(gamma, delta, '.')

    pp.show()