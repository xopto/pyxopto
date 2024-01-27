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

import numpy as np

from xopto.mcvox.mcsource.base import Source
from xopto.mcvox import cltypes, mcobject, mcoptions


class IsotropicVoxel(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClIsotropicVoxel(cltypes.Structure):
            '''
            Ctypes structure that is passed to the Monte Carlo kernel.

            Fields
            ------
            voxel: mc_point3f_t
                The source voxel (center) position as (x, y, z).
            voxel: mc_point3_t
                The source voxel indices as (ind_x, ind_y, ind_z).
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('voxel', T.mc_point3_t),
            ]
        return ClIsotropicVoxel

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource {',
            '	mc_point3f_t position;',
            '	mc_point3_t voxel;',
            '};',
            '',
            'void dbg_print_source(__mc_source_mem const McSource *src);',
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("IsotropicVoxel source:");',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
            '	dbg_print_point3(INDENT  "voxel   :", &src->voxel);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t sin_fi, cos_fi;',
            '	mc_fp_t sin_theta, cos_theta;',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		(source->voxel.x + mc_fmin(mcsim_random(mcsim), FP_1 - FP_EPS))*',
            '			mcsim_voxel_size_x(mcsim) + mcsim_top_left_x(mcsim),',
            '		(source->voxel.y + mc_fmin(mcsim_random(mcsim), FP_1 - FP_EPS))*',
            '			mcsim_voxel_size_y(mcsim) + mcsim_top_left_y(mcsim),',
            '		(source->voxel.z + mc_fmin(mcsim_random(mcsim), FP_1 - FP_EPS))*',
            '			mcsim_voxel_size_z(mcsim) + mcsim_top_left_z(mcsim)',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '',
            '	mcsim_set_direction_coordinates(',
            '		mcsim,',
            '		cos_fi*sin_theta,',
            '		sin_fi*sin_theta,',
            '		cos_theta',
            '	);',
            '',
            '	mcsim_set_weight(mcsim, FP_1);',
            '	dbg_print_status(mcsim, "Launch IsotropicVoxel");',
            '};',
        ))

    def __init__(self, voxel: Tuple[int, int, int]):
        '''
        Isotropic voxel source of photon packets.

        Parameters
        ----------
        voxel: (int, int, int)
            Source voxel position as an integer array-like object of size 3
            (ind_x, ind_y, ind_z).
        '''
        Source.__init__(self)

        self._voxel = np.zeros((3,), dtype=np.int32)
        self._set_voxel(voxel)

    def _get_voxel(self) -> Tuple[int, int, int]:
        return self._voxel
    def _set_voxel(self, voxel: Tuple[int, int, int]):
        self._voxel[:] = voxel
    voxel = property(_get_voxel, _set_voxel, None,
                     'Source voxel indices (ind_x, ind_y, ind_z).')

    def update(self, other: 'IsotropicVoxel' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: xopto.mcvox.mcsource.voxel.IsotropicVoxel or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, IsotropicVoxel):
            self.voxel = other.voxel
        elif isinstance(other, dict):
            self.voxel = other.get('voxel', self.voxel)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See :py:meth:`IsotropicVoxel.cl_type` for
        a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: mcypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: mctypes.Structures
            Filled ctypes structure received as an input argument or a new
            instance if the input argument target is None.
        topgeometry: None
            This source does not use advanced geometry at the top sample surface.
        bottomgeometry: None
            This source does not use advanced geometry at the bottom sample surface.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        if not mc.voxels.isvalid(self._voxel):
            raise ValueError('Voxel index ({}, {}, {}) is not '
                             'valid!'.format(*self._voxel))

        target.voxel.fromarray(self._voxel)
        target.position.fromarray(mc.voxels.center(self._voxel))

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'voxel': self._voxel.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'IsotropicVoxel(voxel=({}, {}, {}))'.format(*self._voxel)


class IsotropicVoxels(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClIsotropicVoxels(cltypes.Structure):
            '''
            Ctypes structure that is passed to the Monte Carlo kernel.

            Fields
            ------
            n: mc_size_t
                Number of voxels.
            offset: mc_size_t
                Offset of the first source voxel data in the read-only data
                buffer. Voxel data are organized as
                (weight_0, index_x_0, index_y_0, index_z_0),
                (weight_1, index_x_1, index_y_1, index_z_1), ...
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('n', T.mc_size_t),
                ('offset', T.mc_size_t),
            ]
        return ClIsotropicVoxels

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource {',
            '	mc_point3f_t position;',
            '	mc_size_t n;',
            '	mc_size_t offset;',
            '};',
            '',
            'void dbg_print_source(__mc_source_mem const McSource *src);',
        ))

    @staticmethod
    def cl_options(mc: mcobject.McObject) -> mcoptions.RawOptions:
        return [('MC_USE_FP_LUT', True)]

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("IsotropicVoxel source:");',
            '	dbg_print_point3(INDENT "n:", &src->position);',
            '	dbg_print_size_t(INDENT "n:", src->n);',
            '	dbg_print_size_t(INDENT  "offset   :", &src->offset);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t sin_fi, cos_fi;',
            '	mc_fp_t sin_theta, cos_theta;',
            '',
            '	mc_size_t src_voxel_index = mc_min(mcsim_random(mcsim)*source->n, source->n - 1);',
            '',
            '	mc_size_t offset = source->offset + src_voxel_index*4;',
            '	mc_fp_t weight = mcsim_fp_lut_array(mcsim)[offset++];',
            '	mc_fp_t voxel_index_x = mcsim_fp_lut_array(mcsim)[offset++];',
            '	mc_fp_t voxel_index_y = mcsim_fp_lut_array(mcsim)[offset++];',
            '	mc_fp_t voxel_index_z = mcsim_fp_lut_array(mcsim)[offset];',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		(voxel_index_x + mc_fmin(mcsim_random(mcsim), FP_1 - FP_EPS))*',
            '			mcsim_voxel_size_x(mcsim) + mcsim_top_left_x(mcsim),',
            '		(voxel_index_y + mc_fmin(mcsim_random(mcsim), FP_1 - FP_EPS))*',
            '			mcsim_voxel_size_y(mcsim) + mcsim_top_left_y(mcsim),',
            '		(voxel_index_z + mc_fmin(mcsim_random(mcsim), FP_1 - FP_EPS))*',
            '			mcsim_voxel_size_z(mcsim) + mcsim_top_left_z(mcsim)',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '',
            '	mcsim_set_direction_coordinates(',
            '		mcsim,',
            '		cos_fi*sin_theta,',
            '		sin_fi*sin_theta,',
            '		cos_theta',
            '	);',
            '',
            '	mcsim_set_weight(mcsim, weight);',
            '	dbg_print_status(mcsim, "Launch IsotropicVoxels");',
            '};',
        ))

    def __init__(self, voxels: np.ndarray, weights: np.ndarray = None):
        '''
        Isotropic voxel source of photon packets.

        Parameters
        ----------
        voxels: np.ndarray
            Array of voxel indices of shape (n, 3), where n is the number
            of voxels. Each row of the array must contain voxel
            indices as (x, y, z).
        weights: np.ndarray
            Initial weight of the photon packet for the corresponding
            source voxel (a value from 0.0 to 1.0 is required).
            Set to 1.0 if None.
        '''
        Source.__init__(self)

        voxels = np.asarray(voxels)
        if voxels.shape[-1] != 3:
            raise ValueError('Shape of the voxel array must be (n, 3)!')

        self._data = None
        self._set_voxels(voxels)
        if weights is not None:
            if weights.size != self._data.shape[0]:
                raise ValueError(
                    'The size of array with voxel weights/intensities '
                    'does not match the array of voxel indices!')
            self._data[:, 0] = np.clip(weights, 0.0, 1.0)

    def _get_voxels(self) -> np.ndarray:
        return self._data[:, 1:]
    def _set_voxels(self, voxels: np.ndarray):
        voxels = np.asarray(voxels)
        if voxels.ndim > 2 or voxels.shape[-1] != 3:
            raise ValueError('Shape of the voxel array must be (n, 3)!')
        if voxels.ndim == 1:
            voxels.shape = (1, 3)
        if self._data is None or self._data.shape[0] != voxels.shape[0]:
            dtype = self._data.dtype if self._data is not None else None
            self._data = np.ones((voxels.shape[0], 4), dtype=dtype)
            self._data[:, 1:] = voxels
        else:
            self._data[:, 1:] = voxels
    voxels = property(_get_voxels, _set_voxels, None,
                      'Indices of the source voxels as a numpy array of '
                      'shape (n, 3), ordered as (x, y, z).')

    def _get_weights(self) -> np.ndarray:
        return self._data[:, 0]
    def _set_weights(self, weights: np.ndarray):
        self._data[:, 0] = np.clip(weights, 0.0, 1.0)
    weights = property(_get_weights, _set_weights, None,
                       'Weight/intensity of the source voxels. '
                       'Must be from interval [0, 1].')

    def update(self, other: 'IsotropicVoxels' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: xopto.mcvox.mcsource.voxel.IsotropicVoxel or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, IsotropicVoxels):
            self.voxels = other.voxels
            self.weights = other.weights
        elif isinstance(other, dict):
            self.voxels = other.get('voxels', self.voxels)
            self.weights = other.get('weights', self.weights)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See :py:meth:`IsotropicVoxels.cl_type` for
        a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: mcypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: mctypes.Structures
            Filled ctypes structure received as an input argument or a new
            instance if the input argument target is None.
        topgeometry: None
            This source does not use advanced geometry at the top sample surface.
        bottomgeometry: None
            This source does not use advanced geometry at the bottom sample surface.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        for voxel in self.voxels:
            if not mc.voxels.isvalid(voxel):
                raise ValueError('Voxel index ({}, {}, {}) is not '
                                 'valid!'.format(*voxel))
        mc_np_type = mc.types.np_float
        if self._data.dtype != mc_np_type:
            self._data = self._data.astype(mc_np_type)

        lut_entry = mc.append_r_lut(np.reshape(self._data, (self._data.size, )))

        center_index = np.mean(self.voxels, 0)
        target.position.fromarray(mc.voxels.center(center_index))
        target.n = self._data.shape[0]

        target.offset = lut_entry.offset

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'voxels': self.voxels.tolist(),
                'weights': self.weights.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'IsotropicVoxels(voxels={}, weights={})'.format(
            self.voxels, self.weights)
