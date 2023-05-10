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

from .base import Sampler, PfSampler
from xopto.mcbase import mcmaterial


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
        sampler: MaterialSampler
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
