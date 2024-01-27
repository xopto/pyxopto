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
from xopto.mcvox import cltypes, mcobject, mctypes
from xopto.mcvox.mcutil import boundary

class Line(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClLine(cltypes.Structure):
            '''
            Ctypes structure that is passed to the Monte Carlo kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position.
            direction_medium: mc_point3f_t
                Source direction in the surrounding medium before refraction
                into the sample.
            direction_sample: mc_point3f_t
                Source direction in the sample after refraction from the
                surrounding medium.
            direction_reflected: mc_point3f_t
                Source direction in the surrounding medium after reflected from
                the sample surface.
            reflectance: mc_fp_t
                Precalculated specular reflectance at the source -> sample
                transition. Only relevant for sources located outside of
                the sample volume.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('direction_medium', T.mc_point3f_t),
                ('direction_sample', T.mc_point3f_t),
                ('direction_reflected', T.mc_point3f_t),
                ('reflectance', T.mc_fp_t)
            ]
        return ClLine

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource {',
            '	mc_point3f_t position;',
            '	mc_point3f_t direction_medium;',
            '	mc_point3f_t direction_sample;',
            '	mc_point3f_t direction_reflected;',
            '	mc_fp_t reflectance;',
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
            '	dbg_print("McSource Line source:");',
            '	dbg_print_point3f(INDENT "position:        ", &src->position);',
            '	dbg_print_point3f(INDENT "direction_medium:", &src->direction_medium);',
            '	dbg_print_point3f(INDENT "direction_sample:", &src->direction_sample);',
            '	dbg_print_point3f(INDENT "direction_reflected:", &src->direction_reflected);',
            '	dbg_print_float(INDENT   "reflectance:     ", src->reflectance);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mcsim_set_weight(mcsim, FP_1 - source->reflectance);',
            '',
            '	mcsim_set_position(mcsim, &source->position);',
            '	mcsim_set_direction(mcsim, &source->direction_sample);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t direction = source->direction_reflected;',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &direction, source->reflectance);',
            '	#endif',
            '};',
        ))

    def __init__(self, position: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Line photon packet source - equivalent of a spatial delta impulse.

        Parameters
        ----------
        position: (float, float, float)
            Line source position as an array-like object of size 3 (x, y, z).

        direction: (float, float, float)
            Line source direction vector as an array-like object of size 3
            (px, py, pz).

        Note
        ----
        If the position lies outside of the sample, the entry point on the
        sample surface will be used to launch the packets. The entry point
        will be determined by propagating the position along the given
        direction. The MC simulation will start after subtracting the specular
        reflectance at the sample boundary from the initial weight of the
        packet.
        Note that in case the position lies within the sample,
        it will be used as the launch point and the packets will retain the
        full initial weight.
        '''
        Source.__init__(self)

        self._position = np.zeros((3,))
        self._direction = np.zeros((3,))
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
        self._direction /= norm
    direction = property(_get_direction, _set_direction, None,
                         'Source direction.')

    def update(self, other: 'Line' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: Line or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, Line):
            self.position = other.position
            self.direction = other.direction
        elif isinstance(other, dict):
            self.position = other.get('position', self.position)
            self.direction = other.get('direction', self.direction)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See :py:meth:`Line.cl_type` for a detailed
        list of fields.

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

        if mc.voxels.contains(self._position):
            search_direction = -self._direction
        else:
            search_direction = self._direction

        position, normal = mc.voxels.intersect(
            self._position, search_direction)
        if position is None:
            raise ValueError(
                'The line source does not intersect the voxelized sample!')

        material_medium = mc.materials[0]
        material_sample = mc.material(position)
        reflectance = boundary.reflectance(
            material_medium.n, material_sample.n,
            np.abs(np.dot(normal, self._direction)))
        if reflectance >= 1.0:
            raise ValueError('The line source is fully reflected from the '
                             'surface of the voxelized sample!')

        sample_direction = boundary.refract(
            self._direction, normal, material_medium.n, material_sample.n)

        dir_reflected = boundary.reflect(self._direction, normal)

        target.position.fromarray(position)
        target.direction_medium.fromarray(self._direction)
        target.direction_sample.fromarray(sample_direction)
        target.direction_reflected.fromarray(dir_reflected)
        target.reflectance = reflectance

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'position': self._position.tolist(), 
                'direction': self._direction.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'Line (position=({}, {}, {}), direction=({}, {}, {}))'.\
            format(*self._position, *self._direction)

