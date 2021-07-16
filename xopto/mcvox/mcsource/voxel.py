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
from xopto.mcvox import cltypes, mcobject


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
                The source voxel indices as (ind_z, ind_y, ind_x).
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
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	printf("IsotropicVoxel source:\\n");',
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
            (ind_z, ind_y, ind_x).
        '''
        Source.__init__(self)

        self._voxel = np.zeros((3,), dtype=np.int)
        self._set_voxel(voxel)

    def _get_voxel(self) -> Tuple[int, int, int]:
        return self._voxel
    def _set_voxel(self, voxel: Tuple[int, int, int]):
        self._voxel[:] = voxel
    voxel = property(_get_voxel, _set_voxel, None,
                        'Source voxel indices.')

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
            raise ValueError('The voxel index ({}, {}, {}) is not '
                             'valid!'.format(*self._voxel))

        target.voxel.fromarray(self._voxel[::-1])
        target.position.fromarray(mc.voxels.center(self._voxel))

        return target, None, None

    def __str__(self):
        return 'IsotropicVoxel(voxel=({}, {}, {}))'.format(*self._voxel)
