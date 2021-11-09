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
from xopto.mcml import cltypes, mctypes, mcobject, mcoptions
from xopto.mcml.mcutil.fiber import MultimodeFiber
from xopto.mcml.mcutil import geometry
from xopto.mcml.mcutil import axis

class SixAroundOnePl(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClSixAroundOnePl(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            transformation: mc_point3f_t
                Transforms coordinates from Monte Carlo to the detector.
            position: mc_point2f_t
                Position of the  center/origin of the detector.
            core_r_squared: mc_fp_t
                Squared radius of the optical fibers.
            core_spacing: mc_fp_t
                Spacing between the optical fibers.
            pl_min: mc_fp_t
                The leftmost edge of the first optical path length accumulator.
            inv_dpl: mc_fp_t
                Inverse of the width of the optical path length accumulators.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle (in the detector
                coordinate space).
            n_pl: mc_size_t
                The number of path length accumulators.
            offset: mc_size_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            pl_log_scale: mc_int_t
                A flag indicating logarithmic scale of the path length axis.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point2f_t),
                ('core_r_squared', T.mc_fp_t),
                ('core_spacing', T.mc_fp_t),
                ('pl_min', T.mc_fp_t),
                ('inv_dpl', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n_pl', T.mc_size_t),
                ('offset', T.mc_size_t),
                ('pl_log_scale', T.mc_int_t),
            ]
        return ClSixAroundOnePl

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_matrix3f_t transformation;'
            '	mc_point2f_t position;'
            '	mc_fp_t core_r_squared;',
            '	mc_fp_t core_spacing;',
            '	mc_fp_t pl_min;',
            '	mc_fp_t inv_dpl;',
            '	mc_fp_t cos_min;',
            '	mc_size_t n_pl;'
            '	mc_size_t offset;',
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
            '	dbg_print("Mc{}Detector - SixAroundOnePl fiber array detector:");'.format(Loc),
            '	dbg_print_matrix3f(INDENT "transformation:", &detector->transformation);',
            '	dbg_print_point2f(INDENT "position:", &detector->position);',
            '	dbg_print_float(INDENT "core_r_squared (mm2):", detector->core_r_squared*1e6f);',
            '	dbg_print_float(INDENT "core_spacing (mm)", detector->core_spacing*1e3f);',
            '	dbg_print_float(INDENT "pl_min (um)", detector->pl_min*1e6f);',
            '	dbg_print_float(INDENT "inv_dpl (1/um)", detector->inv_dpl*1e-6f);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_size_t(INDENT "n_pl:", detector->n_pl);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
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
            '	dbg_print_status(mcsim, "{} SixAroundOnePl fiber array detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	mc_size_t fiber_index = 7; /* invalid index ... no fiber hit */',
            '',
            '	mc_point3f_t rel_pos = {',
            '		pos->x - detector->position.x,',
            '		pos->y - detector->position.y,',
            '		FP_0',
            '	};',
            '	mc_point3f_t mc_pos, detector_pos;',
            '	mc_fp_t dx, dy, r_squared;',
            '',
            '	/* The central detector fiber. */',
            '	mc_pos.x = rel_pos.x;',
            '	mc_pos.y = rel_pos.y;',
            '	mc_pos.z = FP_0;',
            '	transform_point3f(&detector->transformation, &mc_pos, &detector_pos);',
            '	dx = detector_pos.x;',
            '	dy = detector_pos.y;',
            '	r_squared = dx*dx + dy*dy;',
            '	if (r_squared <= detector->core_r_squared){',
            '		/* hit the central fiber */',
            '		fiber_index = 0;',
            '	};',
            '',
            '	/* The 1st or 4th fiber in the ring. */',
            '	mc_pos.x = mc_fabs(rel_pos.x) - detector->core_spacing;',
            '	mc_pos.y = rel_pos.y;',
            '	transform_point3f(&detector->transformation, &mc_pos, &detector_pos);',
            '	dx = detector_pos.x;',
            '	dy = detector_pos.y;',
            '	r_squared = dx*dx + dy*dy;',
            '	if (r_squared <= detector->core_r_squared){',
            '		/* hit the 1st or 4th fiber in the ring */',
            '		fiber_index = (rel_pos.x >= FP_0) ? 1 : 4;',
            '	};',
            '',
            '	/* The 2nd, 3rd, 5th or 6th fiber in the ring. */',
            '	mc_pos.x = mc_fabs(rel_pos.x) - detector->core_spacing*FP_0p5;',
            '	mc_pos.y = mc_fabs(rel_pos.y) - detector->core_spacing*FP_COS_30;',
            '	transform_point3f(&detector->transformation, &mc_pos, &detector_pos);',
            '	dx = detector_pos.x;',
            '	dy = detector_pos.y;',
            '	r_squared = dx*dx + dy*dy;',
            '	if (r_squared <= detector->core_r_squared){',
            '		/* hit the 2nd, 3rd, 5th or 6th fiber in the ring */',
            '		fiber_index = (rel_pos.x >= FP_0) ? ',
            '			((rel_pos.y >= FP_0) ? 2 : 6) : ',
            '			((rel_pos.y >= FP_0) ? 3 : 5);',
            '	}',
            '',
            '	if (fiber_index > 6)',
            '		return;',
            '',
            '	mc_fp_t pl = mcsim_optical_pathlength(mcsim);',
            '	if (detector->pl_log_scale)',
            '	    pl = mc_log(mc_fmax(pl, FP_PLMIN));',
            '	mc_int_t pl_index = mc_int((pl - detector->pl_min)*detector->inv_dpl);',
            '	pl_index = mc_clip(pl_index, 0, detector->n_pl - 1);',
            '',
            '	mc_size_t index = pl_index*7 + fiber_index;',
            '',
            '	address = mcsim_accumulator_buffer_ex(',
            '		mcsim, detector->offset + index);',
            '',
            '	/* Transfor direction vector component z into the detector space. */',
            '	mc_fp_t pz = transform_point3f_z(&detector->transformation, dir);',
            '	dbg_print_float("Packet direction z:", pz);',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(pz));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print("{} SixAroundOnePl fiber array detector depositing:");'.format(Loc),
            '		dbg_print_uint(INDENT "uint weight:", ui32w);',
            '		dbg_print_size_t(INDENT "to fiber index:", fiber_index);',
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};',
        ))

    def cl_options(self, mc, target=None) -> mcoptions.RawOptions:
        '''
        OpenCL kernel options defined by this object.
        '''
        return [('MC_TRACK_OPTICAL_PATHLENGTH', True)]

    def __init__(self, fiber: MultimodeFiber,
                 spacing: float = None,
                 plaxis: axis.Axis = None,
                 position: Tuple[float, float] = (0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        A six-around one optical fiber probe detector with optional
        tilted fibers (direction parameter). The optical fibers are always
        polished in a way to form a tight optical contact with the surface
        of the sample.

        Parameters
        ----------
        fiber: MultimodeFiber
            Properties of the optical fibers.
        spacing: float
            Spacing between the optical fibers. If spacing is None,
            a tight layout is used with spacing set to the outer diameter
            of the fiber cladding.
        plaxis: axis.Axis
            Object that defines the accumulators along the optical path length
            axis (this axis supports log-scale).
        position: (float, float)
            Position of the center of the probe / central finer.
        direction: (float, float, float)
            Reference direction / orientation of the detector. Fibers are
            oriented in this direction and polished to form a tight optical
            contact with the sample (the fiber cross sections are ellipsoids
            if the direction is not perpendicular, i.e different
            from (0, 0, 1).

        Note
        ----
        The first dimension of the accumulator represents the optical
        path length axis, the second dimension represents the optical fibers.
        '''
        if isinstance(fiber, SixAroundOnePl):
            sao = fiber
            fiber = sao.fiber
            spacing = sao.spacing
            position = sao.position
            direction = sao.direction
            nphotons = sao.nphotons
            raw_data = np.copy(sao.raw)
            plaxis = sao.plaxis
        else:
            if plaxis is None:
                plaxis = axis.Axis(0.0, 1.0, 1)
            nphotons = 0
            raw_data = np.zeros((plaxis.n, 7))
            if spacing is None:
                spacing = fiber.dcladding

        super().__init__(raw_data, nphotons)

        self._fiber = fiber
        self._spacing = 0.0
        self._set_spacing(spacing)
        self._direction = np.zeros((3,))
        self._position = np.zeros((2,))
        self._pl_axis = plaxis

        self._set_position(position)
        self._set_direction(direction)

    def fiber_position(self, index: int) -> Tuple[float, float]:
        '''
        Returns the position of the fiber center as a tuple (x, y).

        Parameters
        ----------
        index: int
            Fiber index from 0 to 6. Index zero corresponds to the central
            fiber. The remaining fibers are listed in a counter clockwise
            direction starting with the fiber located on the positive x axis.

        Returns
        -------
        position: (float, float)
            The position of the fiber center as a tuple (x, y).
        '''
        if index >= 7 or index < -7:
            raise IndexError('The fiber index is out of valid range!')

        x = self._spacing*0.5
        y = self._spacing*np.cos(np.pi/6.0)
        if index in (0, -7):
            pos = (0.0, 0.0)
        elif index in (1, -6):
            pos = (self._spacing, 0.0)
        elif index in (2, -5):
            pos = (x, y)
        elif index in (3, -4):
            pos = (-x, y)
        elif index in (4, -3):
            pos = (-self._spacing, 0.0)
        elif index in (5, -2):
            pos = (-x, y)
        else:
            pos = (x, -y)

        return self._position[0] + pos[0], self._position[1] + pos[1]

    def update(self, other: 'SixAroundOnePl' or dict):
        '''
        Update this detector configuration from the other detector. The
        other detector must be of the same type as this detector or a dict with
        appropriate fields.

        Parameters
        ----------
        other: SixAroundOnePl or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, SixAroundOnePl):
            self.fiber = other.fiber
            self.spacing = other.spacing
            self.position = other.position
            self.direction = other.direction
        elif isinstance(other, dict):
            self.fiber = other.get('fiber', self.fiber)
            self.spacing = other.get('spacing', self.spacing)
            self.position = other.get('position', self.position)
            self.direction = other.get('direction', self.direction)

    def _get_fiber(self) -> Tuple[float, float]:
        return self._fiber
    def _set_fiber(self, value: float or Tuple[float, float]):
        self._position[:] = value
    fiber = property(_get_fiber, _set_fiber, None,
                     'Properties of the optical fibers used by the detector.')

    def _get_position(self) -> Tuple[float, float]:
        return self._position
    def _set_position(self, value: float or Tuple[float, float]):
        self._position[:] = value
    position = property(_get_position, _set_position, None,
                       'Position of the detector as a tuple (x, y).')

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

    def _get_spacing(self) -> float:
        return self._spacing
    def _set_spacing(self, value:float):
        self._spacing = float(value)
    spacing = property(_get_spacing, _set_spacing, None,
                       'Spacing between the centers of the central and six '
                       'surrounding optical fiber.')

    def _get_plaxis(self) -> axis.Axis:
        return self._pl_axis
    plaxis = property(_get_plaxis, None, None, 'Path length axis object.')

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

    def check(self):
        '''
        Check if the configuration has errors and raise exceptions if so.
        '''
        if self._spacing < self.fiber.dcore:
            raise ValueError('Spacing between the optical fibers is smaller '
                             'than the diameter of the fiber core!')
        return True

    def _get_normalized(self) -> np.ndarray:
        return self.raw*(1.0/max(self.nphotons, 1.0))
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def cl_pack(self, mc: mcobject.McObject,
                target: cltypes.Structure = None) -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`SixAroundOnePl.cl_type` method for a detailed list
        of fields.

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

        adir = self._direction[0], self._direction[1], abs(self._direction[2])
        T = geometry.transform_base(adir, (0.0, 0.0, 1.0))
        target.transformation.fromarray(T)

        target.core_spacing = self.spacing
        target.core_r_squared = 0.25*self.fiber.dcore**2

        target.cos_min = (1.0 - (self._fiber.na/self._fiber.ncore)**2)**0.5

        target.position.fromarray(self._position)

        target.pl_min = self._pl_axis.scaled_start
        if self._pl_axis.step != 0.0:
            target.inv_dpl = 1.0/self._pl_axis.step
        else:
            target.inv_dpl = 0.0
        target.pl_log_scale = self._pl_axis.logscale
        target.n_pl = self._pl_axis.n

        return target

    def todict(self):
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`SixAroundOnePl.fromdict` method to
        create a new accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'SixAroundOnePl',
            'fiber': self._fiber.todict(),
            'spacing': self._spacing,
            'position':self._position.tolist(),
            'direction':self.direction.tolist(),
            'pl_axis': self.plaxis.todict(),
        }

    @staticmethod
    def fromdict(data: dict):
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`SixAroundOnePl.todict` method.
        '''
        detector_type = data.pop('type')
        if detector_type != 'SixAroundOnePl':
            raise TypeError(
                'Expected a "SixAroundOnePl" type bot got "{}"!'.format(
                    detector_type))
        fiber = data.pop('fiber')
        fiber = MultimodeFiber.fromdict(fiber)
        pl_axis_data = data.pop('pl_axis')
        pl_axis_type = pl_axis_data.pop('type')
        plaxis = getattr(axis, pl_axis_type)(**pl_axis_data)

        return SixAroundOnePl(fiber, plaxis=plaxis, **data)

    def __str__(self):
        return 'SixAroundOne(fiber={}, spacing={}, position=({}, {}), '\
               'direction=({}, {}, {}))'.format(
                   self._fiber, self._spacing, *self._position,
                   *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
