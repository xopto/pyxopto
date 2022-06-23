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
import xopto.mcml.mc


def detector_direction(incidence: float, tilt: float, design_angle: float,
                       verbose: bool = False) -> Tuple[float, float, float]:
    '''
    Determines the direction vector of the detector by solving a quadratic
    equation that arises from the following system of equations:

    .. math::

        x\\sin(\\Theta_i) + z\\cos(\\Theta_i) = -\\cos(\\Theta_0)

        y\\cos(\\Theta_t) - z\\sin(\\Theta_t) = 0

        x^2 + y^2 + z^2 = 1

    In the above equations :math:`\\Theata_i` is the angle of incidence,
    :math:`\\Theta_t` the detector tilt angle and :math:`\\Theta_0` the
    design angle, i.e. the angle between the source and detector.

    The first equation is the dot product of the detector direction vector
    and the direction of the incident beam. The angle between the two is
    the design angle :math:`180 - \\Theta_0`. The second equation is
    the dot product of the detector direction vector and the normal
    of the tilted plane that contains the detector direction vector. Note
    that the tilted plane is the x-z plane rotated around the x axis.

    The solution with negative z component that has a polar angle
    in the x-z plane smaller than the polar angle of the incident source
    is taken (detector always rotates from source towards the x axis).

    Parameters
    ----------
    incidence: float
        Source incidence angle (radians) in the x-z plane, measured from
        the z axis.
    tilt: float
        Detector tilt angle (radians) measured as the angle between the
        tilted plane (x-z plane rotated around the x axis) and the z axis. 
    design_angle: float
        Design angle between the source and detector (radians).
    verbose: bool
        Enables verbose report.

    Returns
    -------
    dir: np.ndarray
        Direction vector of the detector.
    '''

    src_direction = np.array([np.sin(incidence), 0.0, np.cos(incidence)])

    a = 1.0/np.tan(incidence)**2 + np.tan(tilt)**2 + 1.0
    b = 2.0*np.cos(design_angle)*np.cos(incidence)/np.sin(incidence)**2
    c = np.cos(design_angle)**2/np.sin(incidence)**2 - 1.0
    d = (b**2 - 4.0*a*c)**0.5

    z1, z2 = (-b + d)/(2.0*a), (-b - d)/(2.0*a)
    y1, y2 = z1*np.tan(tilt), z2*np.tan(tilt)
    x1 = (-np.cos(design_angle) - z1*np.cos(incidence))/np.sin(incidence)
    x2 = (-np.cos(design_angle) - z2*np.cos(incidence))/np.sin(incidence)

    if verbose:
        print('dir 1             :', x1, y1, z1,)
        print('  src-detector (°):', np.rad2deg(
            np.arccos(-np.dot([x1, y1, z1], src_direction))))
        print('  polar in x-z (°):', np.rad2deg(np.arccos(x1)))
        print('dir 2             :', x2, y2, z2)
        print('  src-detector (°):', np.rad2deg(
            np.arccos(-np.dot([x2, y2, z2], src_direction))))
        print('  polar in x-z (°):', np.rad2deg(np.arccos(x2)))

    if z1 >= 0.0 and z2 >= 0.0:
        raise ValueError('Detector direction could not be resolved!')

    if x1 > np.cos(np.pi*0.5 + incidence):
        direction = np.array([float(x1), float(y1), float(z1)])
    else:
        direction = np.array([float(x2), float(y2), float(z2)])

    return direction


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


class MaterialSampler(Sampler):
    def __init__(self, mua: Sampler, musr: Sampler,
                 n: Sampler, pf: PfSampler):
        '''
        Samples a material.
        
        Parameters
        ----------
        mua: Sampler
            Sampler of the absorption coefficient (1/m).
        musr: Sampler
            Sampler of the reduced scattering coefficient (1/m).
        n: Sampler
            Sampler of the refractive index.
        pf: PfSampler
            Sampler Of the scattering phase function.
        '''
        super().__init__()
        self._mua = mua
        self._musr = musr
        self._n = n
        self._pf = pf

    def reset(self):
        '''
        Reset the states of all the enclosed samplers.
        '''
        self._mua.reset()
        self._musr.reset()
        self._n.reset()
        self._pf.reset()

    def __call__(self, freeze: bool = False) -> dict:
        '''
        Sample the layer parameters.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.

        Returns
        -------
        sample: dict
            A dict object with keys "mua" (float), "mus" (float),
            "n" (float), "d" (float) and "pf" (:py:class:`mcpf.PfBase`).
        '''
        mcpf_obj = self._pf(freeze=freeze)
        pf_obj = mcpf_obj.pf()

        musr = self._musr(freeze=freeze)
        g = pf_obj.g(1)

        return {
            'mua': self._mua(freeze=freeze),
            'musr': musr,
            'mus': musr/(1.0 - g),
            'g': g,
            'n':self._n(freeze=freeze),
            'pf': mcpf_obj
        }

    def _get_pf(self) -> Sampler:
        return self._pf
    pf = property(_get_pf, None, None,
                  'Scattering phase function sampler.')

    def _get_mua(self) -> Sampler:
        return self._mua
    mua = property(_get_mua, None, None,
                   'Absorption coefficient sampler.')

    def _get_musr(self) -> Sampler:
        return self._musr
    musr = property(_get_musr, None, None,
                    'Reduced scattering coefficient sampler.')

    def _get_n(self) -> Sampler:
        return self._n
    n = property(_get_n, None, None,
                 'Refractive index sampler.')

    def _get_d(self) -> Sampler:
        return self._d
    d = property(_get_n, None, None,
                 'Layer thickness sampler.')


    def update(self, layer: mcmaterial.Material, sample: dict = None) \
            -> mcmaterial.Material:
        '''
        Update the material with a new sample.

        Parameters
        ----------
        layer: mcmaterial.Material
            Material to update.
        sample: dict
            Optional sample. If None, a new sample is generated.

        Returns
        -------
        layer: mcmaterial.Material
            Updated material.
        '''
        if sample is None:
            sample = self()

        layer.mua = float(sample['mua'])
        layer.mus = float(sample['mus'])
        layer.n = float(sample['n'])
        layer.pf = sample['pf']

        return layer

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
            'mua': self._mua.todict(),
            'musr': self._musr.todict(),
            'n': self._n.todict(),
            'pf': self._pf.todict()
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'MaterialSampler':
        '''
        Create a new instance of :py:class:`MaterialSampler` from a dict.

        Parameters
        ----------
        data: dict
            Data that were exported by the :py:meth:`MaterialSampler.todict`
            method.

        Returns
        -------
        sampler: LayerSampler
            A new sampler instance.
        '''
        data = dict(data)
        data.pop('type')

        return cls(
            musr=Sampler.fromdict(data.pop('musr')),
            mua=Sampler.fromdict(data.pop('mua')),
            n=Sampler.fromdict(data.pop('n')),
            pf=PfSampler.fromdict(data.pop('pf'))
        )

    def __str__(self):
        return 'MaterialSampler(mua={}, musr={}, n={}, pf={})'.format(
            self._mua, self._mus, self._n, self._pf)

    
class LayerSampler(Sampler):
    def __init__(self, mua: Sampler, musr: Sampler,
                 d: Sampler, n: Sampler, pf: PfSampler):
        '''
        Samples layer parameters.
        
        Parameters
        ----------
        mua: Sampler
            Sampler of the absorption coefficient (1/m).
        musr: Sampler
            Sampler of the reduced scattering coefficient (1/m).
        d: Sampler
            Sampler of the layer thickness.
        n: Sampler
            Sampler of the refractive index.
        pf: PfSampler
            Sampler Of the scattering phase function.
        '''
        super().__init__()
        self._mua = mua
        self._musr = musr
        self._n = n
        self._d = d
        self._pf = pf

    def reset(self):
        '''
        Reset the states of all the enclosed samplers.
        '''
        self._mua.reset()
        self._musr.reset()
        self._n.reset()
        self._d.reset()
        self._pf.reset()

    def __call__(self, freeze: bool = False) -> dict:
        '''
        Sample the layer parameters.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.

        Returns
        -------
        sample: dict
            A dict object with keys "mua" (float), "mus" (float),
            "n" (float), "d" (float) and "pf" (:py:class:`mcpf.PfBase`).
        '''
        mcpf_obj = self._pf(freeze=freeze)
        pf_obj = mcpf_obj.pf()

        musr = self._musr(freeze=freeze)
        g = pf_obj.g(1)

        return {
            'mua': self._mua(freeze=freeze),
            'musr': musr,
            'mus': musr/(1.0 - g),
            'g': g,
            'n':self._n(freeze=freeze),
            'd': self._d(freeze=freeze),
            'pf': mcpf_obj
        }

    def _get_pf(self) -> Sampler:
        return self._pf
    pf = property(_get_pf, None, None,
                  'Scattering phase function sampler.')

    def _get_mua(self) -> Sampler:
        return self._mua
    mua = property(_get_mua, None, None,
                   'Absorption coefficient sampler.')

    def _get_musr(self) -> Sampler:
        return self._musr
    musr = property(_get_musr, None, None,
                    'Reduced scattering coefficient sampler.')

    def _get_n(self) -> Sampler:
        return self._n
    n = property(_get_n, None, None,
                 'Refractive index sampler.')

    def _get_d(self) -> Sampler:
        return self._d
    d = property(_get_n, None, None,
                 'Layer thickness sampler.')


    def update(self, layer: xopto.mcml.mc.mclayer.Layer, sample: dict = None) \
            -> xopto.mcml.mc.mclayer.Layer:
        '''
        Update the layer with a new sample.

        Parameters
        ----------
        layer: xopto.mcml.mc.mclayer.Layer
            Layer to update.
        sample: dict
            Optional sample. If None, a new sample is generated.

        Returns
        -------
        layer: xopto.mcml.mc.mclayer.Layer
            Updated input layer.
        '''
        if sample is None:
            sample = self()

        layer.mua = float(sample['mua'])
        layer.mus = float(sample['mus'])
        layer.d = float(sample['d'])
        layer.n = float(sample['n'])
        layer.pf = sample['pf']

        return layer

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
            'mua': self._mua.todict(),
            'musr': self._musr.todict(),
            'n': self._n.todict(),
            'd': self._d.todict(),
            'pf': self._pf.todict()
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'LayerSampler':
        '''
        Create a new instance of :py:class:`LayerSampler` from a dict.

        Parameters
        ----------
        data: dict
            Data that were exported by the :py:meth:`LayerSampler.todict`
            method.

        Returns
        -------
        sampler: LayerSampler
            A new sampler instance.
        '''
        data = dict(data)
        data.pop('type')

        return cls(
            musr=Sampler.fromdict(data.pop('musr')),
            mua=Sampler.fromdict(data.pop('mua')),
            n=Sampler.fromdict(data.pop('n')),
            d=Sampler.fromdict(data.pop('d')),
            pf=PfSampler.fromdict(data.pop('pf'))
        )

    def __str__(self):
        return 'LayerSampler(mua={}, musr={}, d={}, n={}, pf={})'.format(
            self._mua, self._musr, self._d, self._n, self._pf)


class MultilayerSampler(Sampler):
    def __init__(self, samplers = Tuple[LayerSampler]):
        super().__init__()
        self._samplers = samplers

        pf_type = self._samplers[0].pf.pf_type
        for layer_sampler in self._samplers[1:]:
            if layer_sampler.pf.pf_type != pf_type:
                raise ValueError('All layers must use the same scattering '
                                 'phase function type!')

    def __call__(self, freeze: bool = False) -> List[dict]:
        '''
        Sample the parameters of the layer stack.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.

        Returns
        -------
        samples: List[dict]
            A list of layer samples.
        '''
        return [layer_sampler(freeze=freeze)
                for layer_sampler in self._samplers]

    def reset(self):
        '''
        Reset the states of all the layer samplers.
        '''
        for sampler in self._samplers:
            sampler.reset()

    def update(self, layers: xopto.mcml.mc.mclayer.Layers,
               sample: dict = None) -> xopto.mcml.mc.mclayer.Layers:
        '''
        Update the layer stack with a new sample.

        Parameters
        ----------
        layers: xopto.mcml.mc.mclayer.Layers
            Layer stack to update.
        sample: List[dict]
            Optional sample. If None, a new sample is generated.

        Returns
        -------
        layers: xopto.mcml.mc.mclayer.Layers
            Updated layer stack.
        '''
        if sample is None:
            sample = self()

        for layer_sample, layer in zip(sample, layers):
            layer.mua = float(layer_sample['mua'])
            layer.mus = float(layer_sample['mus'])
            layer.d = float(layer_sample['d'])
            layer.n = float(layer_sample['n'])
            layer.pf = layer_sample['pf']

        return layers

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
            'samplers': [sampler.todict() for sampler in self._samplers]
        }

    def __getitem__(self, index):
        return self._samplers[index]

    @classmethod
    def fromdict(cls, data: dict) -> 'MultilayerSampler':
        '''
        Create a new instance of :py:class:`MultilayerSampler` from a dict.

        Parameters
        ----------
        data: dict
            Data that were exported by the :py:meth:`MultilayerSampler.todict`
            method.

        Returns
        -------
        sampler: MultilayerSampler
            A new sampler instance.
        '''
        data = dict(data)
        data.pop('type')

        return cls(
            samplers=[
                LayerSampler.fromdict(sampler)
                    for sampler in data.pop('samplers')
            ]
        )

    def __len__(self):
        return len(self._samplers)

    def __str__(self):
        return 'MultilayerSampler(samplers={})'.format(
            str(self._samplers))


class IncidenceTiltSampler(Sampler):
    def __init__(self, incidence: Sampler, tilt: Sampler,
                 design_angle: float):
        '''
        Samples the incidence angle of a collimated source and tilt of
        a surface detector from the given distribution. The angle between the
        collimated incident beam and the reference direction of the detector
        (design angle) is kept constant.

        Parameters
        ----------
        incidence: Sampler
            Sampler of the angle of incidence (rad). The angle is measured
            in the incidence plane x-z from the z axis.
            The sampled direction vector of the collimated beam lies
            in the x-z plane.
        tilt: Sampler
            Sampler of the detector tilt angle (rad). The detector direction
            lies in the x-z plane that is tilted/rotated around the x axis.
            The angle of tilt is defined as the angle between the tilted
            plane and the x-z plane.
        design_angle: float
            Angle between the source and detector (rad) that is kept constant.
        '''
        super().__init__()

        self._design_angle = float(design_angle)
        self._incidence = incidence
        self._tilt = tilt

    def reset(self):
        '''
        Reset the states of all the enclosed samplers.
        '''
        self._incidence.reset()
        self._tilt.reset()

    def __call__(self, freeze: bool = False) -> dict:
        '''
        Samples direction of a collimated sources and reference direction
        of the detector.

        Parameters
        ----------
        freeze: bool
            Do not advance the state of the sampler if True.

        Returns
        -------
        sample: dict
            A new sample with the following fields: "source_direction",
            "detector_direction", "incidence" and "tilt".
        '''
        incidence = float(self._incidence(freeze=freeze))
        tilt = float(self._tilt(freeze=freeze))
        detector_dir = detector_direction(
            incidence, tilt, self._design_angle)
        source_dir = np.array([np.sin(incidence), 0.0, np.cos(incidence)])

        return {
            'incidence': incidence, 'tilt': tilt,
            'detector_direction': detector_dir,
            'source_direction': source_dir
        }

    def update(self, src: xopto.mcml.mc.mcsource.Line,
               detector: xopto.mcml.mc.mcdetector.SymmetricX,
               sample: dict = None):
        '''
        Update the Monte Carlo source and detector with the sampled propagation
        directions.

        Parameters
        ----------
        src: xopto.mcml.mc.mcsource.Line
            Source that will be updated with a new direction.
        detector: xopto.mcml.mc.mcdetector.Detector
            Detector that will be updated with the new reference direction.

        Returns
        -------
        src: xopto.mcml.mc.mcsource.Line
            The updated source.
        detector: xopto.mcml.mc.mcdetector.SymmetricX
            The updated detector.
        '''
        if sample is None:
            sample = self()

        src.direction = sample['source_direction']
        detector.direction = sample['detector_direction']

        return src, detector

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
            'incidence': self._incidence.todict(),
            'tilt': self._tilt.todict(),
            'design_angle': self._design_angle
        }
    
    def _get_design_angle(self) -> float:
        return self._design_angle
    design_angle = property(_get_design_angle, None, None,
                            'Design angle (rad).')

    @classmethod
    def fromdict(cls, data: dict) -> 'IncidenceTiltSampler':
        '''
        Create a new instance of :py:class:`IncidenceTiltSampler`
        from a dict.

        Parameters
        ----------
        data: dict
            Data that were exported by the
            :py:meth:`IncidenceTiltSampler.todict` method.

        Returns
        -------
        sampler: IncidenceTiltSampler
            A new sampler instance.
        '''
        data = dict(data)
        data.pop('type')

        return cls(
            incidence=Sampler.fromdict(data.pop('incidence')),
            tilt=Sampler.fromdict(data.pop('tilt')),
            design_angle=data['design_angle']
        )


class MultilayerSfdi:
    def __init__(self, layers: MultilayerSampler,
                 source_detector: IncidenceTiltSampler,
                 sample_count: int = 0):
        '''
        A sampler for multilayer SFDI Monte Carlo simulations.

        Parameters
        ----------
        layers: MultilayerSampler
            Sampler for the parameters of the layer stack.
        source_detector: SfdiSourceDetectorSampler
            Sampler for the parameters of the photon packet source and
            top surface detector.
        sample_count: int
            Zero-based index of the first sample.
        '''
        self._layers = layers
        self._source_detector = source_detector
        self._sample_number = int(sample_count)
        self._rmax_cache = {}

    def reset(self):
        '''
        Reset the states of all the enclosed samplers and set the sample
        number to 0.
        '''
        self._layers.reset()
        self._source_detector.reset()
        self._sample_number = 0

    def update(self, mcobj: xopto.mcml.mc.Mc) -> dict:
        '''
        Update the Monte Carlo simulation parameters with a new
        sample. This call will update the sample index by 1.

        Parameters
        ----------
        mcobj: xopto.mcml.mc.Mc
            Monte Carlo simulator instance.

        Returns
        -------
        sample: dict
            The sample that was used to update the simulator as a dict
            with keys "layers" (a list of samples that were used to update
            the sample layers, "source_detector" (sample that was used to
            update the photon packet source and top surface detector).
        '''
        layers_sample = self._layers()
        self._layers.update(mcobj.layers, layers_sample)

        source_detector_sample = self._source_detector()
        self._source_detector.update(
            mcobj.source, mcobj.detectors.top, source_detector_sample)

        self._sample_number += 1

        return {'layers': layers_sample,
                'source_detector': source_detector_sample}

    def _get__sample_number(self) -> int:
        return self._sample_number
    sample_index = property(_get__sample_number, None, None,
                            'Number of the current sample '
                            '(1 for the first sample)')

    def mc_obj(self, acceptance: float = np.pi*0.5, **kwargs) \
            -> xopto.mcml.mc.Mc:
        '''
        Create and initialize a new Monte Carlo simulation model.

        Parameters
        ----------
        acceptance: float
            Acceptance angle of the detector.
        kwargs: dict
            Optional keyword arguments passed to the
            :py:meth:`~xopto.mcml.mc.Mc.__init__` constructor.

        Returns
        -------
        mcobj: mc.Mc
            A new Monte Carlo simulator instance.
        '''
        # Get a scattering phase function initializer without advancing
        # the sampler state.
        pf_obj = self._layers[0].pf(freeze=True)
        num_layers = len(self._layers)

        layers = []
        for _ in range(num_layers):
            layers.append(xopto.mcml.mc.mclayer.Layer(
                mua=0.0, mus=0.0, d=np.inf, n=1.0, pf=pf_obj)
            )

        mc_layers = xopto.mcml.mc.mclayer.Layers(layers)
        mc_source = xopto.mcml.mc.mcsource.Line()
        mc_top_detector = xopto.mcml.mc.mcdetector.SymmetricX(
            xopto.mcml.mc.mcdetector.SymmetricAxis(
                0.0, 150e-3, 4000, logscale=True),
            cosmin=np.cos(acceptance)
        )
        mc_detectors = xopto.mcml.mc.mcdetector.Detectors(top=mc_top_detector)

        mcobj = xopto.mcml.mc.Mc(mc_layers, mc_source, mc_detectors, **kwargs)

        return mcobj

    def rmax(self, mcobj: xopto.mcml.mc.Mc) -> np.ndarray:
        '''
        Estimate simulation radius for the given simulator instance. 

        Parameters
        ----------
        mcobj: xopto.mcml.mc.Mc
            Monte Carlo simulator instance.

        Returns
        -------
        rmax: np.ndarray
            A numpy vector of simulation radius for all the sample
            layers (excluding the topmost and bottommost layers).
        '''
        rmax_args = (float(mcobj.layers[1].n), float(mcobj.layers[0].n),
                     float(np.arccos(mcobj.detectors.top.cosmin)),
                     25.0e-3)
        rmax_estimator = self._rmax_cache.get(rmax_args)
        if rmax_estimator is None:
            rmax_estimator = mcrmax.RmaxDiffusion(*rmax_args) 
            self._rmax_cache[rmax_args] = rmax_estimator

        num_intrnale_layers = len(mcobj.layers) - 2
        rmax = np.zeros((num_intrnale_layers,))
        for index, layer in enumerate(mcobj.layers[1:-1]):
            g = layer.pf.pf().g(1)
            rmax[index] = rmax_estimator(layer.mua, layer.mus*(1.0 - g))

        return rmax
