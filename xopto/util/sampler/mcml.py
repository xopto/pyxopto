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

from .base import Sampler, PfSampler
from .common import detector_direction
from xopto.mcbase import mcrmax
import xopto.mcml.mc


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
    d = property(_get_d, None, None,
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
        sample: dict
            Sample that was used to update the layer stack.
        '''
        if sample is None:
            sample = self()

        for layer_sample, layer in zip(sample, layers):
            layer.mua = float(layer_sample['mua'])
            layer.mus = float(layer_sample['mus'])
            layer.d = float(layer_sample['d'])
            layer.n = float(layer_sample['n'])
            layer.pf = layer_sample['pf']

        return layers, sample

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
        source_detector: IncidenceTiltSampler
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
        mcobj: xopto.mcml.mc.Mc
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

        num_internal_layers = len(mcobj.layers) - 2
        rmax = np.zeros((num_internal_layers,))
        for index, layer in enumerate(mcobj.layers[1:-1]):
            g = layer.pf.pf().g(1)
            key = (rmax_estimator, float(layer.mua), float(layer.mus*(1.0 - g)))
            rmax_value = self._rmax_cache.get(key)
            if rmax_value is None:
                rmax_value = rmax_estimator(key[1], key[2])
                self._rmax_cache[key] = float(rmax_value)
            rmax[index] = rmax_value
        return rmax
