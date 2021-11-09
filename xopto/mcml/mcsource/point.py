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

from xopto.mcml.mcsource.base import Source
from xopto.mcml import cltypes, mcobject


class IsotropicPoint(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClIsotropicPoint(cltypes.Structure):
            '''
            Ctypes structure that is passed to the Monte Carlo kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position.
            layer_index: mc_size_t
                Index of the layer that contains the source.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('layer_index', T.mc_size_t),
            ]
        return ClIsotropicPoint

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource {',
            '	mc_point3f_t position;',
            '	mc_size_t layer_index;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	printf("IsotropicPoint source:\\n");',
            '	dbg_print_point3f(INDENT "position   :", &src->position);',
            '	dbg_print_size_t(INDENT  "layer_index:", src->layer_index);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_fp_t specular_r = FP_0;',
            '	mc_point3f_t position = source->position;',
            '	mc_point3f_t direction;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '',
            '	direction.x = cos_fi*sin_theta;',
            '	direction.y = sin_fi*sin_theta;',
            '	direction.z = cos_theta;',
            '',
            '	mc_point3f_t refracted_direction = direction;',
            '',
            '	if (position.z <= FP_0) {',
            '		mc_fp_t cc = mc_layer_cc_bottom(mcsim_top_layer(mcsim));',
            '',
            '		if (direction.z > cc) {',
            '			mc_point3f_t normal = {FP_0, FP_0, FP_1};',
            '			mc_fp_t n_sample = mc_layer_n(mcsim_top_sample_layer(mcsim));',
            '			mc_fp_t n_medium = mc_layer_n(mcsim_top_layer(mcsim));',
            '',
            '			refract(&direction, &normal, n_sample, n_medium,',
            '				&refracted_direction);',
            '',
            '			specular_r = reflectance(',
            '				n_medium, n_sample, direction.z, cc);',
            '',
            '			mc_fp_t t = mc_fdiv(-position.z, direction.z);',
            '			position.x += direction.x*t;',
            '			position.y += direction.y*t;',
            '		} else {',
            '			specular_r = FP_1;',
            '		};',
            '',
            '		position.z = FP_0;',
            '	}',
            '',
            '	mcsim_set_direction(mcsim, &refracted_direction);',
            '	mcsim_set_position(mcsim, &position);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &direction, specular_r);',
            '	#endif',
            '',
            '	mcsim_set_weight(mcsim, FP_1 - specular_r);',
            '',
            '	mcsim_set_current_layer_index(mcsim, source->layer_index);',
            '};',
        ))

    def __init__(self, position: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
        '''
        Isotropic point source of photon packets.

        Parameters
        ----------
        position: (float, float, float)
            Source position as an array-like object of size 3 (x, y, z).
            The position must be located above (z <= 0) or within
            the sample (z > 0) but not under the sample.

        Note
        ----
        If the source position lies outside of the sample (z <= 0), the entry
        point into the sample is determined by propagating the packet from the
        source position along the launch direction.
        The MC simulation will start after refracting the packet into the
        sample and subtracting the specular reflectance at the sample boundary
        from the initial weight of the packet. If the photon packet does not
        intersect the sample, the initial weight will be set to 0 (reflectance
        to 1) and the packet will be launched with the z coordinate set to 0.
        Such zero-weight packets are immediately
        terminated and have no contribution to the fluence and surface
        detectors, however will be included in the trace (have no effect on
        the sampling volume or other trace-based analysis due to zero-weight).
        Note that in case the position lies within the sample (z > 0),
        it will be used as the launch point and the packets will retain the
        full initial weight.
        '''
        Source.__init__(self)

        self._position = np.zeros((3,))
        self._set_position(position)

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, position: Tuple[float, float, float]):
        self._position[:] = position
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def update(self, other: 'IsotropicPoint' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: xopto.mcml.mcsource.point.IsotropicPoint or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, IsotropicPoint):
            self.position = other.position
        elif isinstance(other, dict):
            self.position = other.get('position', self.position)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See :py:meth:`IsotropicPoint.cl_type` for
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

        if self._position[2] >= mc.layers.thickness():
            raise ValueError('The source must be located above or within '
                             'the sample but not under the sample!')

        target.position.fromarray(self._position)

        if self._position[2] <= 0.0:
            target.layer_index = 1
        else:
            target.layer_index = layer_index = mc.layer_index(self._position[2])

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'position': self._position.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'IsotropicPoint(position=({}, {}, {}))'.format(*self._position)
