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

from xopto.mccyl.mcsource.base import Source
from xopto.mccyl import cltypes, mcobject


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
                Index of the layer that contains the source. Only valid if
                the source is located inside the sample or is outside and the
                photon packet is refracted into the sample.
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
            '	mc_fp_t initial_weight = FP_1;',
            '	mc_size_t layer_index = 1;',
            '	mc_point3f_t initial_pos = source->position;',
            '	mc_point3f_t source_dir;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '',
            '	source_dir.x = cos_fi*sin_theta;',
            '	source_dir.y = sin_fi*sin_theta;',
            '	source_dir.z = cos_theta;',
            '',
            '	mc_point3f_t initial_dir = source_dir;',
            '',
            '	mc_fp_t r2 = initial_pos.x*initial_pos.x + ',
            '		initial_pos.y*initial_pos.y;',
            '	mc_fp_t r_sample = mc_layer_r_inner(mcsim_layer(mcsim, 0));',
            '',
            '	if (r2 >= r_sample*r_sample) {',
            '		mc_fp_t d1, d2;',
            '',
            '		int has_intersection = ray_cylinder_intersections(',
            '			r_sample, &initial_pos, &source_dir, &d1, &d2);',
            '',
            '		if (has_intersection){',
            '			mc_fp_t k = mc_fmin(d1, d2);',
            '			initial_pos.x += k*source_dir.x;',
            '			initial_pos.y += k*source_dir.y;',
            '			initial_pos.z += k*source_dir.z;',
            '',
            '			mc_point3f_t normal;',
            '			radial_normal_inv(&initial_pos, &normal);',
            '			mc_fp_t cos1 = normal.x*source_dir.x + ',
            '				normal.y*source_dir.y;',
            '			mc_fp_t specular_r = reflectance(',
            '				mc_layer_n(mcsim_layer(mcsim, 0)),',
            '				mc_layer_n(mcsim_layer(mcsim, 1)), cos1,',
            '				mc_layer_cc_inner(mcsim_layer(mcsim, 0))',
            '			);'
            '			if (cos1 > mc_layer_cc_inner(mcsim_layer(mcsim, 0))) {',
            '				refract(',
            '					&source_dir, &normal,',
            '					mc_layer_n(mcsim_layer(mcsim, 0)),',
            '					mc_layer_n(mcsim_layer(mcsim, 1)),',
            '					&initial_dir',
            '				);',
            '			}',
            '			layer_index = 1;',
            '			initial_weight = FP_1 - specular_r;',
            '',
            '			#if MC_USE_SPECULAR_DETECTOR',
            '			if(specular_r > FP_0){',
            '				mc_point3f_t reflected_dir;',
            '				reflect(&source_dir, &normal, &reflected_dir);',
            '				mcsim_specular_detector_deposit(',
            '					mcsim, &initial_pos, &reflected_dir, specular_r);',
            '			}',
            '			#endif',
            '		} else {',
            '			initial_pos = (mc_point3f_t){.x=FP_0, .y=FP_0, .z=FP_0};',
            '			layer_index = mcsim_layer_count(mcsim) - 1;',
            '			initial_weight = FP_0;',
            '		}',
            '	} else {',
            '		initial_weight = FP_1;',
            '		layer_index = source->layer_index;',
            '	};',
            '',
            '	mcsim_set_current_layer_index(mcsim, layer_index);',
            '	mcsim_set_weight(mcsim, initial_weight);',
            '	mcsim_set_position(mcsim, &initial_pos);',
            '	mcsim_set_direction(mcsim, &initial_dir);',
            '',
            '	dbg_print_status(mcsim, "Launch IsotropicPoint");',
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
        other: xopto.mccyl.mcsource.point.IsotropicPoint or dict
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

        target.position.fromarray(self._position)

        r = (self._position[0]**2 + self._position[1]**2)**0.5
        if r >= mc.layers.diameter()*0.5:
            target.layer_index = 1
        else:
            target.layer_index = mc.layer_index(r)

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'position': self._position.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'IsotropicPoint(position=({}, {}, {}))'.format(*self._position)
