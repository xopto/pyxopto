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
from xopto.mcml import cltypes, mcobject, mcoptions
from xopto.mcml.mcutil import axis


class RadialPl(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClRadialPl(cltypes.Structure):
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
            pl_min: mc_fp_t
                The leftmost edge of the first optical path length accumulator.
            inv_dpl: mc_fp_t
                Inverse of the width of the optical path length accumulators.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle (relative to the
                direction of the detector).
            n_r: mc_size_t
                The number of concentric ring accumulators.
            n_pl: mc_size_t
                The number of path length accumulators.
            offset: mc_size_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            r_log_scale: mc_int_t
                A flag indicating logarithmic scale of the radial axis.
            pl_log_scale: mc_int_t
                A flag indicating logarithmic scale of the path length axis.

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
                ('l_min', T.mc_fp_t),
                ('inv_dpl', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n_r', T.mc_size_t),
                ('n_pl', T.mc_size_t),
                ('offset', T.mc_size_t),
                ('r_log_scale', T.mc_int_t),
                ('pl_log_scale', T.mc_int_t),
            ]
        return ClRadialPl

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
            '	mc_fp_t pl_min;',
            '	mc_fp_t inv_dpl;',
            '	mc_fp_t cos_min;',
            '	mc_size_t n_r;',
            '	mc_size_t n_pl;',
            '	mc_size_t offset;',
            '	mc_int_t r_log_scale;',
            '	mc_int_t pl_log_scale;',
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
            '	dbg_print("Mc{}Detector - RadialPl detector:");'.format(Loc),
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_point2f(INDENT "position:", &detector->position);',
            '	dbg_print_float(INDENT "r_min (mm):", detector->r_min*1e3f);',
            '	dbg_print_float(INDENT "inv_dr (1/mm)", detector->inv_dr*1e-3f);',
            '	dbg_print_float(INDENT "pl_min (um):", detector->r_min*1e6f);',
            '	dbg_print_float(INDENT "inv_dpl (1/um)", detector->inv_dr*1e-6f);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_size_t(INDENT "n_r:", detector->n_r);',
            '	dbg_print_size_t(INDENT "n_pl:", detector->n_pl);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '	dbg_print_int(INDENT "r_log_scale:", detector->r_log_scale);',
            '	dbg_print_int(INDENT "pl_log_scale:", detector->pl_log_scale);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} RadialPl detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const struct Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	mc_fp_t dx = pos->x - detector->position.x;',
            '	mc_fp_t dy = pos->y - detector->position.y;',
            '	mc_fp_t r = mc_sqrt(dx*dx + dy*dy);',
            '',
            '	if (detector->r_log_scale)',
            '		r = mc_log(mc_fmax(r, FP_RMIN));',
            '',
            '	mc_int_t r_index = mc_int((r - detector->r_min)*detector->inv_dr);',
            '	r_index = mc_clip(r_index, 0, detector->n_r - 1);',
            '',
            '	mc_fp_t pl = mcsim_optical_pathlength(mcsim);',
            '	if (detector->pl_log_scale)',
            '	    pl = mc_log(mc_fmax(pl, FP_PLMIN));',
            '	mc_int_t pl_index = mc_int((pl - detector->pl_min)*detector->inv_dpl);',
            '	pl_index = mc_clip(pl_index, 0, detector->n_pl - 1);',
            '',
            '	mc_size_t index = pl_index*detector->n_r + r_index;',
            '',
            '	address = mcsim_accumulator_buffer_ex(',
            '		mcsim, detector->offset + index);',
            '',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(dot3f(dir, &detector->direction)));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} RadialPl detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def cl_options(self, mc, target=None) -> mcoptions.RawOptions:
        '''
        OpenCL kernel options defined by this object.
        '''
        return [('MC_TRACK_OPTICAL_PATHLENGTH', True)]

    def __init__(self, raxis: axis.Axis or axis.RadialAxis,
                 plaxis: axis.Axis = None,
                 position: Tuple[float, float] = (0.0, 0.0),
                 cosmin: float = 0.0,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Radial path-length reflectance-transmittance accumulator.

        Parameters
        ----------
        raxis: axis.Axis or axis.RadialAxis
            Object that defines the accumulators along the radial axis
            (this axis supports log-scale).
        plaxis: axis.Axis
            Object that defines the accumulators along the optical path length
            axis (this axis supports log-scale).
        position: (float, float)
            Position of the center of the radial accumulator.
        cosmin: float
            Cosine of the maximum acceptance angle (relative to the direction
            of the detector) of the accumulator.
        direction: (float, float, float)
            Reverence direction / orientation of the detector.

        Note
        ----
        The first dimension of the accumulator represents the optical
        path length axis, the second dimension represents the radial axis. 
        '''
        if isinstance(raxis, RadialPl):
            radialpl = raxis
            position = radialpl.position
            raxis = type(radialpl.raxis)(radialpl.raxis)
            plaxis = type(radialpl.plaxis)(radialpl.plaxis)
            cosmin = radialpl.cosmin
            direction = radialpl.direction
            raw_data = np.copy(radialpl.raw)
            nphotons = radialpl.nphotons
        else:
            if plaxis is None:
                plaxis = axis.Axis(0.0, 1.0, 1)

            raw_data = np.zeros((plaxis.n, raxis.n))
            nphotons = 0

        super().__init__(raw_data, nphotons)

        self._position = np.zeros((2,))
        self._cosmin = 0.0
        self._direction = np.zeros((3,))
        self._r_axis = raxis
        self._pl_axis = plaxis
        self._set_position(position)
        self._set_cosmin(cosmin)
        self._set_direction(direction)
        self._inv_accumulators_area = 1.0/(np.pi*(self._r_axis.edges[1:]**2 -
                                                  self._r_axis.edges[:-1]**2))
        self._inv_accumulators_area.shape = (1, self._inv_accumulators_area.size)

    def _get_raxis(self) -> axis.Axis or axis.RadialAxis:
        return self._r_axis
    raxis = property(_get_raxis, None, None, 'Radial axis object.')

    def _get_plaxis(self) -> axis.Axis:
        return self._pl_axis
    plaxis = property(_get_plaxis, None, None, 'Path length axis object.')

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
    r = property(_get_r, None, None,
                 'Centers of the radial axis accumulators.')

    def _get_redges(self) -> np.ndarray:
        return self._r_axis.edges
    redges = property(_get_redges, None, None,
                     'Edges of the radial axis accumulators.')

    def _get_nr(self) -> int:
        return self._r_axis.n
    nr = property(_get_nr, None, None,
                  'Number of accumulators in the radial axis.')

    def _get_rlogscale(self) -> bool:
        return self._r_axis.logscale
    rlogscale = property(_get_rlogscale, None, None, 'Radial axis scale.')

    def _get_pl(self):
        return self._pl_axis.centers
    pl = property(_get_pl, None, None,
                  'Centers of the optical pathlength axis accumulators.')

    def _get_pledges(self):
        return self._pl_axis.edges
    pledges = property(_get_pledges, None, None,
                       'Edges of the optical pathlength axis accumulators.')

    def _get_npl(self):
        return self._pl_axis.n
    npl = property(_get_npl, None, None,
                   'Number of accumulators in the optical pathlength axis.')

    def _get_normalized(self) -> np.ndarray:
        return self.raw*self._inv_accumulators_area*(1.0/self.nphotons)
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def cl_pack(self, mc: mcobject.McObject, \
                target : cltypes.Structure = None) -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`RadialPl.cl_type` method for a detailed list of fields.

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
        target.r_log_scale = self._r_axis.logscale
        target.n_r = self._r_axis.n

        target.pl_min = self._pl_axis.scaled_start
        if self._pl_axis.step != 0.0:
            target.inv_dpl = 1.0/self._pl_axis.step
        else:
            target.inv_dpl = 0.0
        target.pl_log_scale = self._pl_axis.logscale
        target.n_pl = self._pl_axis.n

        target.cos_min = self._cosmin

        target.direction.fromarray(self._direction)

        return target

    def todict(self) -> dict:
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
            'type':'RadialPl',
            'position':self._position.tolist(),
            'r_axis':self._r_axis.todict(),
            'pl_axis':self._pl_axis.todict(),
            'cosmin':self._cosmin,
            'direction':self._direction.tolist()
        }

    @staticmethod
    def fromdict(data: dict) -> 'RadialPl':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`RadialPl.todict` method.
        '''
        rt_type = data.pop('type')
        if rt_type != 'RadialPl':
            raise TypeError('Expected "RadialPl" type bot got "{}"!'.format(rt_type))
        r_axis_data = data.pop('r_axis')
        r_axis_type = r_axis_data.pop('type')
        pl_axis_data = data.pop('pl_axis')
        pl_axis_type = pl_axis_data.pop('type')

        return RadialPl(
            getattr(axis, r_axis_type)(**r_axis_data),
            getattr(axis, pl_axis_type)(**pl_axis_data), 
            **data)

    def __str__(self):
        return 'RadialPl(raxis={}, plaxis={}, position=({}, {}), cosmin={}, '\
               'direction=({}, {}, {}))'.format(
                   self._r_axis, self._pl_axis, *self._position, self._cosmin,
                   *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
