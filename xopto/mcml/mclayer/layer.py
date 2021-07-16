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

from typing import List

from xopto.mcml import mcobject
from xopto.mcml import mctypes
from xopto.mcml.mcutil import boundary
from xopto.mcml import cltypes
from xopto.mcml import mcpf


class Layer(mcobject.McObject):
    '''
    Class that represents a single sample layer.
    '''

    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Returns a structure data type that is used to represent one Layer
        instance in the OpenCL kernel of the Monte Carlo simulator.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.

        Returns
        -------
        opencl_t: ClLayer
            OpenCL Structure that represents a layer. 
        '''
        T = mc.types
        class ClLayer(cltypes.Structure):
            _fields_ = [
                ('thickness', T.mc_fp_t),
                ('top', T.mc_fp_t),
                ('bottom', T.mc_fp_t),
                ('n', T.mc_fp_t),
                ('cos_critical_top', T.mc_fp_t),
                ('cos_critical_bottom', T.mc_fp_t),
                ('mus', T.mc_fp_t),
                ('mua', T.mc_fp_t),
                ('inv_mut', T.mc_fp_t),
                ('mua_inv_mut', T.mc_fp_t),
                ('pf', self.pf.fetch_cl_type(mc))
            ]

        return ClLayer

    def __init__(self, d: float, n: float, mua: float, mus: float,
                 pf: mcpf.PfBase):
        '''
        Layer object constructor.

        Parameters
        ----------
        d: float
            Layer thickness (m).
        n: float
            Index of refraction.
        mua: float
            Absorption coefficient (1/m).
        mus: float
            Scattering (NOT reduced) coefficient (1/m).
        pf: mcpf.PfBase
            Scattering phase function object that is derived from PhBase
            class.


        The physical properties of the layer can be read or changed through
        member properties:

            - d: float - 
              Layer thickness (m).
            - n: float - 
              Index of refraction.
            - mua: float - 
              Absorption coefficient (1/m).
            - mus: float - 
              Scattering (NOT reduced) coefficient (1/m).
            - pf: mcpf.PfBase - 
              Scattering phase function object that is derived from PhBase
              class.

        Note
        ----
        The layer boundaries in the z direction are increasing in the
        direction of the layer stack.
        Z coordinate 0.0 bellongs to the top surface of the first sample layer!
        '''
        self._d = float(d)
        self._n = float(n)
        self._mua = float(mua)
        self._mus = float(mus)
        self._pf = pf

    def _set_d(self, d: float):
        self._d = float(d)
    def _get_d(self) -> float:
        return self._d
    d = property(_get_d, _set_d, None, 'Layer thickens (m).')

    def _set_n(self, n: float):
        self._n = float(n)
    def _get_n(self) -> float:
        return self._n
    n = property(_get_n, _set_n, None, 'Refractive index of the layer.')

    def _set_mua(self, mua: float):
        self._mua = float(mua)
    def _get_mua(self) -> float:
        return self._mua
    mua = property(_get_mua, _set_mua, None,
                   'Absorption coefficient of the layer (1/m).')

    def _set_mus(self, mus: float):
        self._mus = float(mus)
    def _get_mus(self) -> float:
        return self._mus
    mus = property(_get_mus, _set_mus, None,
                   'Scattering coefficient of the layer (1/m).')

    def _get_pf(self) -> mcpf.PfBase:
        return self._pf
    def _set_pf(self, pf: mcpf.PfBase):
        if type(self._pf) != type(pf):
            raise ValueError('The scattering phase function type '\
                             'of the layer must not change!')
        self._pf = pf
    pf = property(_get_pf, _set_pf, None, 'Phase function object.')

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'d':self._d, 'n':self._n, 'mua':self._mua, 'mus':self._mus,
                'pf':self._pf, 'type':'Layer'}

    @classmethod
    def fromdict(cls, data: dict) -> 'Layer':
        '''
        Create a new object from dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        t = data.pop('type')
        if t != 'Layer':
            raise ValueError('Cannot create a Layer instance from the data!')
        return cls(**data)

    def __str__(self):
        return 'Layer(d={}, n={}, mua={}, mus={}, pf={})'.format(
            self._d, self._n, self._mua, self._mus, self._pf)

    def __repr__(self):
        return  '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))


class Layers(mcobject.McObject):
    '''
    Class that represents a stack of layers forming the sample.

    Note
    ----
    The topmost and bottommost layers of the stack are used to describe the
    medium that surrounds the sample top and bottom surfaces, respectively.
    Therefore, at least three layers must be always specified,
    namely two layers of the surrounding medium and one sample layer!
    The thicknesses of the topmost and bottommost layers will be automatically
    set to infinity regardless of the layer thickness set by the user.
    '''
    def __init__(self, layers: List[Layer] or 'Layers'):
        '''
        Constructs a managed sample layer stack from a list of sample layers.

        Note
        ----
        The topmost and bottommost layers of the stack are used to describe the
        medium that surrounds the sample top and bottom surfaces, respectively.
        Therefore, at least three layers must be always specified,
        namely two layers of the surrounding medium and one sample layer!

        The bottom surface of the topmost layer (the surrounding medium) is
        located at z=0. The positive direction of the z axis points in the
        direction of the layer stack.

        The thicknesses of the topmost and bottommost layers will be
        automatically set to infinity when passed to the OpenCL kernel
        (regardless of the layer thickness set by the user).

        Note that all layers must use the same scattering phase function
        model.

        Parameters
        ----------
        layers: list or Layers
            A list of sample layers. Requires at least 3 items!
        '''
        self._pf_type = None

        if isinstance(layers, Layers):
            self._layers = layers.tolist()
            self._pf_type = type(self._layers[0].pf)
        else:
            self._layers = list(layers)
            self.check()

    def check(self):
        '''
        Check if the layers are consistent and using a single
        scattering phase function type.
        Raises exception on error.
        '''
        if len(self._layers) < 3:
                ValueError('At least three layers are required, '
                           'but got only {:d}!'.format(len(self._layers)))

        if self._pf_type is None:
            self._pf_type = type(self._layers[1].pf)

        for layer in self._layers:
            if not isinstance(layer, Layer):
                raise TypeError('All layers must be instances of Layer '
                                'but found {:s}!'.format(type(layer).__name__))
            if type(layer.pf) != self._pf_type:
                raise TypeError(
                    'All the sample layer must use the same scattering phase '
                    'function model! Found {} and {}!'.format(
                        self._pf_type.__name__, type(layer.pf).__name__))


    def layer(self, index: int) -> Layer:
        '''
        Returns layer at the specified index. Note that the first layer
        (index 0) and the last layer (index -1) represent the medium
        sourounding the top and bottom surfaces of the sample, respectively
        '''
        return self._layers[index]

    def layer_index(self, z: float) -> int:
        '''
        Returns the layer index that contains the given z coordinate. Note that
        the layer includes the top surface boundary but not the bottom surface
        boundary, i.e. the z extent of the layer is [z_top, z_bottom), where
        z_bottom > z_top.

        Parameters
        ----------
        z: float
            Z coordinate of a point.

        Returns
        -------
        layer_index: int
            Index of the layer that contains the give poit z.
        '''
        index = len(self._layers) - 1
        bottom = 0.0
        for pos, layer in enumerate(self._layers):
            if pos > 0:
                bottom += layer.d
            if z < bottom:
                index = pos
                break
        return index

    def thickness(self) -> float:
        '''
        Thickness of the layer stack excluding the topmost and bottommost
        layers that sourround the sample.

        Returns
        -------
        thickness: float
            The sample thickness excluding the topmost and bottommost layers of
            the surrounding medium.
        '''
        d = 0.0
        for layer in self._layers[1:-1]:
            d += layer.d
        return d

    def cl_options(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL options of the scattering phase function.
        '''
        return self._layers[1].pf.fetch_cl_options(mc)

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return self._layers[1].pf.fetch_cl_declaration(mc)

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return self._layers[1].pf.fetch_cl_implementation(mc)

    def cl_type(self, mc: mcobject.McObject) -> cltypes.Array:
        '''
        Returns an OpenCL array of ClLayer structures that is used to
        represent one instance of a layer stack.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.

        Returns
        ------- 
        clarray: cltypes.Structure*len(self)
            Array of ClLayers.
        '''
        return self._layers[0].fetch_cl_type(mc)*len(self._layers)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Array = None) \
            -> cltypes.Array:
        '''
        Pack the layers into an OpenCL data type. The OpenCL data
        type is returned by the :py:meth:`Layers.cl_type` method.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure*len(self)
            A structure representing an array of Layers.

        Returns
        -------
        target: cltypes.Structure
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''

        num_layers = len(self._layers)

        if target is None or len(target) != num_layers:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        for index, layer in enumerate(self._layers):
            cc_top = cc_bottom = 0.0
            if index > 0:
                cc_top = boundary.cos_critical(
                    layer.n, self._layers[index - 1].n)
            if index + 1 < num_layers:
                cc_bottom = boundary.cos_critical(
                    layer.n, self._layers[index + 1].n)

            mut = layer.mua + layer.mus
            if mut > 0.0:
                inv_mut = 1.0/mut
            else:
                inv_mut = float('inf')
            if layer.mus == 0.0:
                mua_inv_mut = 1.0
            else:
                mua_inv_mut = layer.mua*inv_mut

            target[index].thickness = layer.d
            if index == 0:
                target[index].top = -float('inf')
                target[index].bottom = 0.0
                target[index].thickness = float('inf')
            elif index == num_layers - 1:
                target[index].top = target[index - 1].bottom.value
                target[index].bottom = float('inf')
                target[index].thickness = float('inf')
            else:
                target[index].top = target[index - 1].bottom.value
                target[index].bottom = target[index - 1].bottom.value + layer.d
            target[index].n = layer.n
            target[index].cos_critical_top = cc_top
            target[index].cos_critical_bottom = cc_bottom
            target[index].mua = layer.mua
            target[index].mus = layer.mus
            target[index].inv_mut = inv_mut
            target[index].mua_inv_mut = mua_inv_mut

            layer.pf.cl_pack(mc, target[index].pf)

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'layers':self._layers, 'type': 'Layers'}

    def tolist(self) -> List[Layer]:
        '''
        Returns a weak copy of the list of managed layers.

        Returns
        -------
        layers: list[layers]
            List of managed layers
        '''
        return list(self._layers)

    @classmethod
    def fromdict(cls, data: dict) -> 'Layers':
        '''
        Create a new Layers object from a dict. The dict keys must match
        the parameter names defined by the constructor.
        '''
        T = data.pop('type')
        if T != 'Layers':
            raise ValueError(
            'Cannot create a Layers instance from the given data!')
        return cls(**data)

    def __getitem__(self, what):
        return self._layers[what]

    def __setitem__(self, what, value):
        self._layers[what] = value

    def __len__(self):
        return len(self._layers)

    def __str__(self):
        layers = ['    ' + str(layer) for layer in self._layers]
        layers[0] +=  ",  # medium above the sample"
        layers[-1] += "   # medium bellow the sample"
        layers_str = ',\n'.join(layers)
        return 'Layers([\n{}\n])'.format(layers_str)

    def __repr__(self):
        return  '{:s} # id 0x{:>08X}.'.format(self.__str__(), id(self))
