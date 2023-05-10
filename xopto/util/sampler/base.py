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

from typing import Tuple, List, Type

import numpy as np

from xopto.mcbase import mcrmax, mcpf, mcmaterial


class Sampler:
    @classmethod
    def fromdict(cls, data: dict) -> 'Sampler':
        '''
        Create a new instance from the data in the dict.

        Parameters
        ----------
        data: dict
            A sampler instance exported to dict.

        Returns
        -------
        sampler: Sampler
            A new sampler instance.
        '''
        data = dict(data)
        type_name = data.pop('type')
        sampler_type = globals().get(type_name)
        return sampler_type(**data)

    def reset(self):
        '''
        Reset the state of the sampler. Reimplement this method for
        custom handling of the sampler state.
        '''
        pass

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))


class SequenceSampler(Sampler):
    def __init__(self, sequence: np.ndarray, start: int = 0):
        '''
        A sequence sampler.

        Parameters
        ----------
        sequence: np.ndarray
            A sequence of vales to sample. The array is sampled from a flat
            view.
        start: int
            Zero-based index of the first sample.

        Note
        ----
        After reaching the end of sequence, the sampling will continue with
        the first element in the sequence.
        '''
        super().__init__()
        self._sequence = np.asarray(sequence)
        self._pos = self._start = int(start)

    def reset(self):
        '''
        Reset the sequence sampler to start from the first item as specified
        by the start parameter of the constructor.
        '''
        self._pos = self._start

    def _get_pos(self) -> int:
        return self._pos
    pos = property(_get_pos, None, None,
                   'Zero-based index of the next sample.')

    def _get_sequence(self) -> int:
        return self._sequence
    sequence = property(_get_sequence, None, None,
                        'The sampled data sequence.')

    def __call__(self, freeze: bool = False) -> float:
        '''
        Return a sample from the sequence.


        Parameters
        ----------
        freeze: bool
            Do not advance the state/position of the sampler if True.

        Returns
        -------
        sample: float or np.ndarray
            The sample.
        '''
        sample = self._sequence.flat[self._pos]

        if not freeze:
            self._pos += 1
            if self._pos >= self._sequence.size:
                self._pos = 0

        return sample

    def __str__(self):
        return 'SequenceSampler(sequence={}, start={})'.format(
            self._sequence, self._start)


class UniformSampler(Sampler):
    def __init__(self, start: float, stop: float, logscale: bool = False):
        '''
        Randomly samples values from interval [start, stop]. Note that the
        values of start and stop can be provided in arbitrary order.

        Parameters
        ----------
        start: float
            Start of the sampling interval.
        stop: float
            Stop/end of the sampling interval.
        logscale: bool
            Use logarithmic sampling on the start-stop interval.
        '''
        super().__init__()
        self._set_interval((start, stop))
        self._set_logscale(logscale)

    def _get_interval(self) -> Tuple[float, float]:
        return self._interval
    def _set_interval(self, interval: Tuple[float, float]):
        self._interval = (float(interval[0]), float(interval[1]))
    interval = property(_get_interval, None, None,
                        'Sampling interval as a tuple (start, stop)')

    def _get_start(self) -> float:
        return self._interval[0]
    def _set_start(self, start: float):
        self._interval[0] = float(start)
    start = property(_get_start, _set_start, None,
                     'Start of the sampling interval.')

    def _get_stop(self) -> float:
        return self._interval[1]
    def _set_stop(self, stop: float):
        self._interval[1] = float(stop)
    stop = property(_get_stop, _set_stop, None,
                    'Stop of the sampling interval.')

    def _get_min(self) -> float:
        return min(self._interval)
    min = property(_get_min, None, None,
                   'Lower bound of the sampling interval.')

    def _get_max(self) -> float:
        return max(self._interval)
    max = property(_get_max, None, None,
                   'Upper bound of the sampling interval.')

    def _get_logscale(self) -> bool:
        return self._logscale
    def _set_logscale(self, state: bool):
        self._logscale = bool(state)
    logscale = property(_get_logscale, _set_logscale, None,
                        'Sampling in logarithmic scale.')

    def __call__(self, n: int or Tuple[int] = 1, freeze: bool = False) \
            -> np.ndarray:
        '''
        Return the requested number of samples.

        Parameters
        ----------
        n: int or Tuple[int]
            The requested number of samples.
        freeze: bool
            Do not advance the state of the sampler if True.

        Returns
        -------
        samples: float or np.ndarray
            The samples.
        '''
        data = np.random.rand(n)
        if self._logscale:
            log_min, log_max = np.log10(self.min), np.log10(self.max)
            data *= (log_max - log_min)
            data += log_min
            data = np.float_power(10, data, out=data)
        else:
            data *= (self.max - self.min)
            data += self.min

        return data

    def todict(self) -> dict:
        '''
        Export object to a dict.

        Returns
        -------
        data: dict
            Instance data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'start': self.start,
            'stop': self.stop,
            'logscale': self.logscale
        }

    def __str__(self):
        return 'UniformSampler(start={:g}, stop={:g}, logscale={})'.format(
            *self._interval, self._logscale)


class NormalSampler(Sampler):
    def __init__(self, mean: float, sigma: float, clip: float = 5):
        '''
        Sampler that follows a Normal distribution with the given
        mean and standard deviation.

        Parameters
        ----------
        mean: float
            Mean value of the distribution.
        sigma: float
            Standard deviation of the distribution.
        clip: float
            Clip the values to :math:`mean ± clip \\cdot sigma`.
        '''
        super().__init__()

        self._mean = 0.0
        self._sigma = 1.0
        self._clip = float(np.inf)

        self._set_mean(mean)
        self._set_sigma(sigma)
        self._set_clip(clip)

    def _get_mean(self) -> float:
        return self._mean
    def _set_mean(self, value: float):
        self._mean = float(value)
    mean = property(_get_mean, _set_mean, None,
                    'Mean value of the distribution.')

    def _get_sigma(self) -> float:
        return self._sigma
    def _set_sigma(self, value: float):
        self._sigma = float(value)
    sigma = property(_get_sigma, _set_sigma, None,
                    'Standard deviation of the distribution.')

    def _get_clip(self) -> float:
        return self._clip
    def _set_clip(self, value: float):
        self._clip = float(value)
    clip = property(_get_clip, _set_clip, None,
                    'Clip value range to mean ± clip*sigma.')

    def __call__(self, n: int = 1, freeze: bool = False) -> np.ndarray:
        '''
        Return the requested number of samples.

        Parameters
        ----------
        n: int or Tuple[int]
            The requested number of samples.
        freeze: bool
            Do not advance the state of the sampler if True.

        Returns
        -------
        samples: float or np.ndarray
            The samples.
        '''
        data = np.random.normal(self._mean, self._sigma, n)
        delta = self._clip*self._sigma
        np.clip(data, self._mean - delta, self._mean + delta)

        return data

    def todict(self) -> dict:
        '''
        Export object to a dict.

        Returns
        -------
        data: dict
            Instance data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'mean': self.mean,
            'sigma': self.sigma,
            'clip': self.clip
        }

    def __str__(self):
        return 'NormalSampler(mean={:g}, sigma={:g}, clip={:g})'.format(
            self._mean, self._sigma, self._clip)


class ConstantSampler(UniformSampler):
    def __init__(self, value: float):
        '''
        A sampler that returns a constant value.

        Parameters
        ----------
        value: float
            Constant value that will be returned by the sampler.
        '''
        super().__init__(value, value)

    def __call__(self, n: int or Tuple[int] = 1, freeze: bool = False) \
            -> np.ndarray:
        '''
        Return the requested number of samples.

        Parameters
        ----------
        n: int or Tuple[int]
            The requested number of samples.
        freeze: bool
            Do not advance the state of the sampler if True.

        Returns
        -------
        samples: float or np.ndarray
            The samples.
        '''
        return np.tile(self.min, n)

    def _get_value(self) -> float:
        return self.start
    def _set_value(self, value: float):
        self.interval = (value, value)
    value = property(_get_value, _set_value, None,
                     'Value of the sampler.')

    def todict(self) -> dict:
        '''
        Export object to a dict.

        Returns
        -------
        data: dict
            Instance data exported to a dict.
        '''
        return {
            'type': self.__class__.__name__,
            'value': self.value,
        }

    def __str__(self):
        return 'ConstantSampler(value={:g})'.format(self.value)


class PfSampler(Sampler):
    def __init__(self, pf_type: str or Type[mcpf.PfBase],
                 pf_args: Tuple[Sampler]):
        '''
        Monte Carlo simulator scattering phase function sampler.

        Parameters
        ----------
        pf_type: mcpf.PfBase
            Scattering phase function type.
        pf_args: Tuple[Sampler]
            Samplers of the scattering phase function arguments.
        '''
        super().__init__()
        if isinstance(pf_type, str):
            pf_type = getattr(mcpf, pf_type)

        self._pf_type = pf_type
        self._pf_args = tuple(pf_args)

    def _get_pf_type(self) -> mcpf.PfBase:
        return self._pf_type
    pf_type = property(_get_pf_type, None, None,
                       'Scattering phase function type used by the sampler.')

    def _get_pf_args(self) -> Tuple[Sampler]:
        return self._pf_args
    pf_args = property(_get_pf_args, None, None,
                       'Samplers of the scattering phase function arguments.')

    def reset(self):
        '''
        Reset the states of samplers that sample the arguments
        of the scattering phase function.
        '''
        for arg in self._pf_args:
            arg.reset()

    def __call__(self, freeze: bool = False) -> mcpf.PfBase:
        '''
        Sample tha scattering phase function and return a new instance.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.
        '''
        return self._pf_type(*[arg(freeze=freeze) for arg in self._pf_args])

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
            'pf_type': self._pf_type.__name__,
            'pf_args': [arg.todict() for arg in self._pf_args]
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'PfSampler':
        '''
        Create a new instance from a dict.

        Parameters
        ----------
        data: dict
            Instance data that were exported to a dict.
        '''
        data = dict(data)
        data.pop('type')
        return cls(
            pf_type=data.pop('pf_type'),
            pf_args=[Sampler.fromdict(arg) for arg in data.pop('pf_args')]
        )

    def __str__(self):
        return 'PfSampler(pf_type={}, pf_args={})'.format(
            self._pf_type.__name__, self._pf_args)
