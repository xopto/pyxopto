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


class Gk2Sampler(PfSampler):
    def __init__(self,
                 gamma: Tuple[float, float] or None = None,
                 delta: Tuple[float, float] or None = None):
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
        '''
        self._gamma = self._delta =  None
        if gamma is not None:
            self._gamma = (float(gamma[0]), float(gamma[1]))
        if delta is not None:
            self._delta = (float(delta[0]), float(delta[1]))

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
            gs = sample.gs(3)
            gamma, delta = g2gamma(gs[1], gs[2]), g2delta(gs[1], gs[3])
            if self._gamma is not None:
                if self._gamma[0] > gamma or gamma > self._gamma[1]:
                    if freeze:
                        raise RuntimeError('Cannot freeze a rejected sample!')
                    continue
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
        return 'Gk2Sampler(gamma={}, delta={})'.format(self._gamma, self._delta)
