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
from xopto.pf.util import GkMap, GkPolygon, g2gamma, g2delta


class GkSampler(PfSampler):
    def __init__(self,
                 gamma: Tuple[float, float] or None = None,
                 delta: Tuple[float, float] or None = None):
        '''
        Uniformly samples the native parameters of the GK scattering phase
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
        '''
        self._gamma = self._delta =  None
        if gamma is not None:
            self._gamma = (float(gamma[0]), float(gamma[1]))
        if delta is not None:
            self._delta = (float(delta[0]), float(delta[1]))

        g_sampler = UniformSampler(0.0, 0.98)
        a_sampler = UniformSampler(-0.5, 10.0)

        super().__init__(mcpf.Gk, (g_sampler, a_sampler))

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
            gamma, delta = g2gamma(gs[1], gs[2]), g2delta(gs[1], gs[3])
            # check gamma range if required
            if self._gamma is not None:
                if self._gamma[0] > gamma or gamma > self._gamma[1]:
                    if freeze:
                        raise RuntimeError('Cannot freeze a rejected sample!')
                    continue
            # check delta range if required
            if self._delta is not None:
                if self._delta[0] > delta or delta > self._delta[1]:
                    if freeze:
                        raise RuntimeError('Cannot freeze a rejected sample!')
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
            'delta': self._delta
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'GkSampler':
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
            raise TypeError('Expected data for GkSampler but got {}!'.format(T))

        return cls(**data)

    def __str__(self):
        return 'GkSampler(gamma={}, delta={})'.format(self._gamma, self._delta)


class GkGammaDeltaSampler(PfSampler):
    def __init__(self,
                 gamma: Tuple[float, float] or None = None,
                 delta: Tuple[float, float] or None = None):
        '''
        Uniformly samples the GK scattering phase function in
        the gamma-delta plane. All samples with gamma or delta value
        outside of the specified range are rejected. 

        Parameters
        ----------
        gamma: Tuple[float, float] or None
            Range of gamma values that will be sampled. If None, the range
            is computed from delta values. If both gamma and delta are None,
            the full range of the GK scattering phase function is used.
        delta: Tuple[float, float] or None:
            Range of delta values that will be sampled. If None, the range
            is computed from the gamma values. If both gamma and delta are None,
            the full range of the GK scattering phase function is used.
        '''
        self._map = GkMap.fromfile()
        self._polygon = GkPolygon.fromfile()

        if gamma is None and delta is None:
            gamma_boundary, delta_boundary = self._polygon.boundary()
            gamma = gamma_boundary.min(), gamma_boundary.max()
            delta = delta_boundary.min(), delta_boundary.max()
        else:
            gamma, delta = self._polygon.bounding_box(gamma, delta)

        self._gamma, self._delta = gamma, delta
        gamma_sampler = UniformSampler(*gamma)
        delta_sampler = UniformSampler(*delta)

        super().__init__(mcpf.Gk, (gamma_sampler, delta_sampler))

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
            if bool(self._polygon.contains(gamma=gamma, delta=delta)):
                break
            elif freeze:
                raise RuntimeError('Cannot freeze illegal sample value!')

        raw_args = self._map.invgammadelta(gamma, delta)
        pf_sample = self._pf_type(*raw_args)
        gs = pf_sample.pf().gs(3)
        # gamma_, delta_ = g2gamma(gs[1], gs[2]), g2delta(gs[1], gs[3])
        # print(gamma, delta, gamma_, delta_)
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
            'delta': self._delta
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'GkGammaDeltaSampler':
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
                'Expected data for GkGammaDeltaSampler but got {}!'.format(T))

        return cls(**data)

    def __str__(self):
        return 'GkGammaDeltaSampler(gamma={}, delta={})'.format(
            self._gamma, self._delta)
