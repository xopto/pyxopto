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

from xopto.mcml.mcdetector.base import Detector
from xopto.mcml import cltypes, mctypes, mcobject
from xopto.mcml.mcutil import axis


class Radial(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClRadial(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            direction: mc_point3f_t
                Reference direction/orientation of the detector.
            position: mc_point2f_t
                Position of the  center/origin of the radial detector.
            r_min: mc_fp_t
                The leftmost edge of the first concentric ring accumulator.
            inv_dr: mc_fp_t
                Inverse of the width of the concentric ring accumulators.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle (relative to the
                direction of the detector).
            n: mc_size_t
                The number of concentric ring accumulators.
            offset: mc_size_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            log_scale: mc_int_t
            A flag indicating logarithmic scale of the accumulator.

            Note
            ----
            Note that for logarithmic radial accumulators the values of r_min
            and inv_dr are passed in a logarithmic scale.
            '''
            _fields_ = [
                ('direction', T.mc_point3f_t),
                ('position', T.mc_point2f_t),
                ('r_min', T.mc_fp_t),
                ('inv_dr', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n', T.mc_size_t),
                ('offset', T.mc_size_t),
                ('log_scale', T.mc_int_t),
            ]
        return ClRadial

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_point3f_t direction;'
            '	mc_point2f_t position;'
            '	mc_fp_t r_min;',
            '	mc_fp_t inv_dr;',
            '	mc_fp_t cos_min;',
            '	mc_size_t n;',
            '	mc_size_t offset;',
            '	mc_int_t log_scale;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the detector accumulator in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_detector(__mc_detector_mem const Mc{}Detector *detector){{'.format(loc, Loc),
            '	dbg_print("Mc{}Detector - Radial detector:");'.format(Loc),
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_point2f(INDENT "position:", &detector->position);',
            '	dbg_print_float(INDENT "r_min (mm):", detector->r_min*1e3f);',
            '	dbg_print_float(INDENT "inv_dr (1/mm)", detector->inv_dr*1e-3f);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_size_t(INDENT "n:", detector->n);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '	dbg_print_int(INDENT "log_scale:", detector->log_scale);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} Radial detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const struct Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	mc_fp_t dx = pos->x - detector->position.x;',
            '	mc_fp_t dy = pos->y - detector->position.y;',
            '	mc_fp_t r = mc_sqrt(dx*dx + dy*dy);',
            '',
            '	if (detector->log_scale)',
            '		r = mc_log(mc_fmax(r, FP_RMIN));',
            '',
            '	mc_int_t r_index = mc_int((r - detector->r_min)*detector->inv_dr);',
            '	size_t index = mc_clip(r_index, 0, detector->n - 1);',
            '',
            '	address = mcsim_accumulator_buffer_ex(',
            '		mcsim, detector->offset + index);',
            '',
            '	mc_point3f_t detector_direction = detector->direction;',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(mc_dot_point3f(dir, &detector_direction)));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} Radial detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def __init__(self, raxis: axis.Axis or axis.RadialAxis,
                 position: Tuple[float, float] = (0.0, 0.0),
                 cosmin: float = 0.0,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Radial reflectance-transmittance accumulator.

        Parameters
        ----------
        raxis: axis.Axis or axis.RadialAxis
            Object that defines the accumulators along the radial axis
            (this axis supports log-scale).
        position: (float, float)
            Position of the center of the radial accumulator.
        cosmin: float
            Cosine of the maximum acceptance angle (relative to the direction
            of the detector) of the accumulator.
        direction: (float, float, float)
            Reverence direction / orientation of the detector.
        '''
        if isinstance(raxis, Radial):
            radial = raxis
            position = radial.position
            raxis = type(radial.raxis)(radial.raxis)
            cosmin = radial.cosmin
            direction = radial.direction
            raw_data = np.copy(radial.raw)
            nphotons = radial.nphotons
        else:
            raw_data = np.zeros((raxis.n,))
            nphotons = 0

        super().__init__(raw_data, nphotons)

        self._position = np.zeros((2,))
        self._cosmin = 0.0
        self._direction = np.zeros((3,))
        self._r_axis = raxis
        self._set_position(position)
        self._set_cosmin(cosmin)
        self._set_direction(direction)
        self._inv_accumulators_area = 1.0/(np.pi*(self._r_axis.edges[1:]**2 -
                                                  self._r_axis.edges[:-1]**2))

    def _get_raxis(self) -> axis.Axis or axis.RadialAxis:
        return self._r_axis
    raxis = property(_get_raxis, None, None, 'Radial axis object.')

    def _get_position(self) -> Tuple[float, float]:
        return self._position
    def _set_position(self, value: float or Tuple[float, float]):
        self._position[:] = value
    position = property(_get_position, _set_position, None,
                       'Position of the radial accumulator as a tuple (x, y).')

    def _get_cosmin(self) -> Tuple[float, float]:
        return self._cosmin
    def _set_cosmin(self, value: float or Tuple[float, float]):
        self._cosmin = min(max(float(value), 0.0), 1.0)
    cosmin = property(_get_cosmin, _set_cosmin, None,
                      'Cosine of the maximum acceptance angle.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('Direction vector norm/length must not be 0!')
        self._direction *= 1.0/norm
    direction = property(_get_direction, _set_direction, None,
                         'Detector reference direction.')

    def _get_r(self) -> np.ndarray:
        return self._r_axis.centers
    r = property(_get_r, None, None, 'Centers of the accumulators.')

    def _get_edges(self) -> np.ndarray:
        return self._r_axis.edges
    edges = property(_get_edges, None, None, 'Edges of the accumulators.')

    def _get_n(self) -> int:
        return self._r_axis.n
    n = property(_get_n, None, None, 'Number of accumulators.')

    def _get_logscale(self) -> bool:
        return self._r_axis.logscale
    logscale = property(_get_logscale, None, None, 'Axis scale.')

    def _get_normalized(self) -> np.ndarray:
        k = 1.0/max(self._nphotons, 1)
        return self.raw*self._inv_accumulators_area*k
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def cl_pack(self, mc: mcobject.McObject,
                target : cltypes.Structure = None) -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`Radial.cl_type` method for a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        allocation = mc.cl_allocate_rw_accumulator_buffer(self, self.shape)
        target.offset = allocation.offset

        target.position.fromarray(self._position)
    
        target.r_min = self._r_axis.scaled_start
        if self._r_axis.step != 0.0:
            target.inv_dr = 1.0/self._r_axis.step
        else:
            target.inv_dr = 0.0
        target.log_scale = self._r_axis.logscale
        target.n = self._r_axis.n

        target.cos_min = self._cosmin

        target.direction.fromarray(self._direction)

        return target

    def todict(self):
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`Radial.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'Radial',
            'position':self._position.tolist(),
            'r_axis':self._r_axis.todict(),
            'cosmin':self._cosmin,
            'direction':self._direction.tolist()
        }

    @staticmethod
    def fromdict(data):
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`Radial.todict` method.
        '''
        data = dict(data)
        rt_type = data.pop('type')
        if rt_type != 'Radial':
            raise TypeError('Expected "Radial" type bot got "{}"!'.format(rt_type))
        r_axis_data = data.pop('r_axis')
        r_axis_type = r_axis_data.pop('type')

        return Radial(getattr(axis, r_axis_type)(**r_axis_data), **data)

    def __str__(self):
        return 'Radial(raxis={}, position=({}, {}), cosmin={}, '\
               'direction=({}, {}, {}))'.format(
                   self._r_axis, *self._position, self._cosmin,
                   *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
