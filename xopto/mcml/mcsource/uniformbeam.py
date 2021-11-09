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
from xopto.mcml import cltypes, mcobject, mctypes
from xopto.mcml.mcutil import boundary, geometry


class UniformBeam(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformBeam(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.
            
            Parameters
            ----------
            mc: McObject
                A Monte Carlo simulator instance.

            Returns
            -------
            struct: cltypes.Structure
                A structure type that represents the uniform beam in
                the Monte Carlo kernel.

            Fields
            ------
            transformation: mc_matrix3f_t
                Transformation from the beam coordinate system to the Monte
                Carlo coordinate system.
            position: mc_point3f_t
                Source position (beam axis).
            direction: mc_point3f_t
                Source direction (beam axis) in the sample medium (after
                refraction).
            radius: mc_point2f_t
                Ellipse radii along the x and y axis. Use equal values for a
                circular beam.
            reflectance: mc_fp_t
                Precalculated reflectance at the source -> sample boundary.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point3f_t),
                ('direction', T.mc_point3f_t),
                ('radius', T.mc_point2f_t),
                ('reflectance', T.mc_fp_t),
            ]
        return ClUniformBeam

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource{',
            '	mc_matrix3f_t transformation;',
            '	mc_point3f_t position;',
            '	mc_point3f_t direction;',
            '	mc_point2f_t radius;',
            '	mc_fp_t reflectance;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return  '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	printf("UniformBeam source:\\n");',
            '	printf(INDENT "position: (%.3f, %.3f, %.3f) mm\\n",',
            '		src->position.x*1e3f, src->position.y*1e3f, src->position.z*1e3f);',
            '	printf(INDENT "direction: (%.3f, %.3f, %.3f)\\n",',
            '		src->direction.x, src->direction.y, src->direction.z);',
            '	printf(INDENT "radius: (%.3f, %.3f) mm\\n", src->radius.x*1e3f, src->radius.y*1e3f);',
            '	printf(INDENT "reflectance: %.3f\\n", src->reflectance);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '   __mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t rand_sqrt = mc_sqrt(mcsim_random(mcsim));',
            '	mc_fp_t fi = FP_2PI*mcsim_random(mcsim);',
            '	mc_point3f_t pt_src, pt_mc;',
            '',
            '	pt_src.x = rand_sqrt*mc_cos(fi)*source->radius.x;',
            '	pt_src.y = rand_sqrt*mc_sin(fi)*source->radius.y;',
            '	pt_src.z = FP_0;',
            '',
            '	transform_point3f(&source->transformation, &pt_src, &pt_mc);',
            '	mc_fp_t k = mc_fdiv(FP_0 - pt_mc.z, source->direction.z);',
            '	pt_mc.x += k*source->direction.x;',
            '	pt_mc.y += k*source->direction.y;',
            '	pt_mc.z = FP_0;',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + pt_mc.x,',
            '		source->position.y + pt_mc.y,',
            '		FP_0',
            '	);',
            '	mcsim_set_direction(mcsim, &source->direction);',
            '	mcsim_set_weight(mcsim, FP_1 - source->reflectance);',
            '	mcsim_set_current_layer_index(mcsim, 1);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = (mc_point3f_t){FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			mc_layer_n(mcsim_layer(mcsim, 0)), &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, source->reflectance);',
            '	#endif',
            '',
            '	dbg_print_status(mcsim, "Launch UniformBeam");',
            '};',
        ))

    def __init__(self, diameter: float or Tuple[float, float],
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Uniform intensity collimated beam photon packet source.

        Parameters
        ----------
        diameter: float or (float, float)
            Collimated beam diameter. Or diameters of the ellipse
            along the x and y axis
        position: (float, float, float)
            Center of the collimated beam as an array-like object of size 3
            (x, y, z). The beam will be always propagated to the top
            sample surface.
        direction: (float, float, float)
            Direction of the collimated beam as an array-like object of size 3
            (px, py, pz). The vector should be normalized to unit length and
            have a positive z coordinate (hitting the top sample surface).

        Note
        ----
        The beam will be first propagated from the given position to the
        entry point on the sample surface along the propagation
        direction (no interactions with the medium during this step).
        Note that in case the position lies within the sample, the
        beam will be propagated to the entry point using reversed direction.
        From there it will be refracted into the sample. The MC simulation
        will start after subtracting the specular reflectance at the
        sample boundary from the initial weight of the packet.
        '''
        Source.__init__(self)

        self._position = np.zeros((3,))
        self._direction = np.zeros((3,))
        self._direction[2] = 1.0
        self._diameter = np.zeros((2,))

        self._set_diameter(diameter)
        self._set_position(position)
        self._set_direction(direction)

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, position: Tuple[float, float, float]):
        self._position[:] = position
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('The norm/length of the propagation direction '
                             'vector must not be 0!')
        self._direction *= 1.0/norm
        if self._direction[-1] <= 0.0:
            raise ValueError('Z component of the propagation direction '
                             'must be positive!')
    direction = property(_get_direction, _set_direction, None,
                        'Source direction.')

    def _get_diameter(self) -> float:
        return self._diameter
    def _set_diameter(self, diameter: float or Tuple[float, float]):
        self._diameter[:] = diameter
        self._diameter = np.maximum(0.0, self._diameter)
        if np.any(self._diameter < 0.0):
            raise ValueError('Beam diameter must not be negative!')
    diameter = property(_get_diameter, _set_diameter, None,
                        'Beam diameter along the x and y axis (m).')

    def update(self, other: dict or 'UniformBeam'):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: UniformBeam or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, UniformBeam):
            self.diameter = other.diameter
            self.position = other.position
            self.direction = other.direction
        elif isinstance(other, dict):
            self.diameter = other.get('diameter', self.diameter)
            self.position = other.get('position', self.position)
            self.direction = other.get('direction', self.direction)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills a structure (target) with the data required by the
        Monte Carlo simulator kernel.
        See the :py:meth:`UniformBeam.cl_type` for a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: pyopyo.mcml.mcsource.UniformBeam.cl_type
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: pyopyo.mcml.mcsource.UniformBeam.cl_type
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

        # propagate the beam to the top sample surface
        k = (0.0 - self._position[2])/self._direction[2]
        position = self._position + k*self._direction
        position[2] = 0.0

        direction = boundary.refract(self._direction, (0.0, 0.0, 1.0),
                                     mc.layers[0].n, mc.layers[1].n)

        reflectance = boundary.reflectance(
            mc.layers[0].n, mc.layers[1].n, abs(self._direction[-1]))

        T = geometry.transform_base((0.0, 0.0, 1.0), self._direction)

        target.transformation.fromarray(T)

        target.position.fromarray(position)
        target.direction.fromarray(direction)

        target.radius.x = self._diameter[0]*0.5
        target.radius.y = self._diameter[1]*0.5

        target.reflectance = reflectance

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'diameter': self._diameter.tolist(), 
                'position': self._position.tolist(),
                'direction': self._direction.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'UniformBeam(diameter=({}, {}), position=({}, {}, {}). ' \
               'direction=({}, {}, {}))'.format(
            *self._diameter, *self._position, *self._direction)
