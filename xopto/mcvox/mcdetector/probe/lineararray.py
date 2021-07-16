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

from xopto.mcvox.mcdetector.base import Detector
from xopto.mcvox import cltypes, mctypes, mcobject
from xopto.mcvox.mcutil.fiber import MultimodeFiber
from xopto.mcvox.mcutil import geometry


class LinearArray(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClLinearArray(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            transformation: mc_point3f_t
                Transforms coordinates from Monte Carlo to the detector.
            first_position: mc_point2f_t
                The center of the first optical fiber in the array.
            delta_position: mc_point2f_t
                Distance vector between two neighboring optical fibers.
            core_r_squared: mc_fp_t
                Squared radius of the optical fibers.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle (in the detector
                coordinate space).
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('first_position', T.mc_point2f_t),
                ('delta_position', T.mc_point2f_t),
                ('core_r_squared', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('offset', T.mc_size_t),
            ]
        return ClLinearArray

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_matrix3f_t transformation;'
            '	mc_point2f_t first_position;'
            '	mc_point2f_t delta_position;'
            '	mc_fp_t core_r_squared;',
            '	mc_fp_t cos_min;',
            '	mc_size_t offset;',
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
            '	dbg_print("Mc{}Detector - LinearArray fiber array detector:");'.format(Loc),
            '	dbg_print_matrix3f(INDENT "transformation:", &detector->transformation);',
            '	dbg_print_point2f(INDENT "first_position:", &detector->first_position);',
            '	dbg_print_point2f(INDENT "delta_position:", &detector->delta_position);',
            '	dbg_print_float(INDENT "core_r_squared (mm2):", detector->core_r_squared*1e6f);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} LinearArray fiber array detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	mc_size_t fiber_index = {}; /* invalid index ... no fiber hit */'.format(self._n),
            '',
            '	mc_fp_t dx, dy, r_squared;',
            '	mc_point3f_t mc_pos, detector_pos;',
            '	mc_fp_t fiber_x = detector->first_position.x;',
            '	mc_fp_t fiber_y = detector->first_position.y;',
            '',
            '	pragma_unroll_hint({})'.format(self._n),
            '	for(mc_size_t index=0; index < {}; ++index){{'.format(self._n),
            '		mc_pos.x = pos->x - fiber_x;',
            '		mc_pos.y = pos->y - fiber_y;',
            '		mc_pos.z = FP_0;',
            '',
            '		transform_point3f(&detector->transformation, &mc_pos, &detector_pos);',
            '		dx = detector_pos.x;',
            '		dy = detector_pos.y;',
            '		r_squared = dx*dx + dy*dy;',
            '',
            '		if (r_squared <= detector->core_r_squared){',
            '			/* hit this fiber */',
            '			fiber_index = index;',
            '			break;',
            '		};',
            '',
            '		fiber_x += detector->delta_position.x;',
            '		fiber_y += detector->delta_position.y;',
            '	};',
            '',
            '	if (fiber_index >= {})'.format(self._n),
            '		return;',
            '',
            '	address = mcsim_accumulator_buffer_ex(',
            '		mcsim, detector->offset + fiber_index);',
            '',
            '	/* Transfor direction vector component z into the detector space. */',
            '	mc_fp_t pz = transform_point3f_z(&detector->transformation, dir);',
            '	dbg_print_float("Packet direction z:", pz);',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(pz));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print("{} LinearArray fiber array detector depositing:");'.format(Loc),
            '		dbg_print_uint(INDENT "uint weight:", ui32w);',
            '		dbg_print_size_t(INDENT "to fiber index:", fiber_index);',
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};',
        ))

    def __init__(self, fiber: MultimodeFiber, n=1,
                 spacing: float = None, orientation: Tuple[float, float] = (1.0, 0.0),
                 position: Tuple[float, float] = (0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Optical fiber probe detector with a linear array of optical fibers that
        are optionally tilted(direction parameter). The optical fibers are
        always polished in a way that forms a tight optical contact with
        the surface of the sample.

        Parameters
        ----------
        fiber: MultimodeFiber
            Optical properties of the fibers.
        n: int
            The number of optical fibers in the linear array. This option
            is a compile-time feature and caanot be changed once the detector
            object is created.
        spacing: float
            Spacing between the optical fibers. If spacing is None,
            a tight layout is used with spacing set to the outer diameter
            of the fiber cladding.
        orientation: (float, float)
            Vector that points in the direction of the linear fiber array.
            By default the fibers are place in the direction of x axis,
            i.e. vector (1.0, 0.0).
            The orientation must point in the direction from the first to
            the last optical fiber!
        position: (float, float)
            Position of the center of the linear fiber array.
        direction: (float, float, float)
            Reference direction / orientation of the detector fibers. Fibers
            are oriented in this direction and polished to form a tight optical
            contact with the sample (the fiber cross sections are ellipsoids
            if the direction is not perpendicular, i.e different
            from (0, 0, 1).
        '''
        if isinstance(fiber, LinearArray):
            la = fiber
            fiber = la.fiber
            n = la.n
            spacing = la.spacing
            orientation = la.orientation
            position = la.position
            direction = la.direction
            nphotons = la.nphotons
            raw_data = np.copy(la.raw)
        else:
            nphotons = 0
            n = max(int(n), 1)
            raw_data = np.zeros((n,))
            if spacing is None:
                spacing = fiber.dcladding

        super().__init__(raw_data, nphotons)

        self._n = 0
        self._fiber = None
        self._spacing = 0.0
        self._orientation = np.array((1.0, 0.0))
        self._direction = np.array((0.0, 0.0, 1.0))
        self._position = np.array((0.0, 0.0))

        self._set_fiber(fiber)
        self._n = max(int(n), 1)
        self._set_spacing(spacing)
        self._set_orientation(orientation)
        self._set_position(position)
        self._set_direction(direction)

    def update(self, other: 'LinearArray' or dict):
        '''
        Update this detector configuration from the other detector. The
        other detector must be of the same type as this detector or a dict with
        appropriate fields.

        Parameters
        ----------
        other: LinearArray or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, LinearArray):
            self.fiber = other.fiber
            self.spacing = other.spacing
            self.orientation = other.orientation
            self.position = other.position
            self.direction = other.direction
        elif isinstance(other, dict):
            self.fiber = other.get('fiber', self.fiber)
            self.spacing = other.get('spacing', self.spacing)
            self.orientation = other.get('orientation', self.orientation)
            self.position = other.get('position', self.position)
            self.direction = other.get('direction', self.direction)

    def _get_fiber(self) -> Tuple[float, float]:
        return self._fiber
    def _set_fiber(self, value: float or Tuple[float, float]):
        self._fiber = value
    fiber = property(_get_fiber, _set_fiber, None,
                     'Properties of the optical fibers used by the detector.')

    def _get_n_fiber(self) -> int:
        return self._n
    n = property(_get_n_fiber, None, None,
                 'Number of optical fiber in the linear array.')

    def _get_spacing(self) -> float:
        return self._spacing
    def _set_spacing(self, value:float):
        self._spacing = float(value)
    spacing = property(_get_spacing, _set_spacing, None,
                       'Spacing between the centers of the optical fibers')

    def _get_orientation(self) -> Tuple[float, float]:
        return self._orientation
    def _set_orientation(self, orientation: Tuple[float, float]):
        self._orientation[:] = orientation
        norm = np.linalg.norm(self._orientation)
        if norm == 0.0:
            raise ValueError('Orientation vector norm/length must not be 0!')
        self._orientation *= 1.0/norm
    orientation = property(_get_orientation, _set_orientation, None,
                         'Orientation / direction of the linear fiber array.')

    def _get_position(self) -> Tuple[float, float]:
        return self._position
    def _set_position(self, value: float or Tuple[float, float]):
        self._position[:] = value
    position = property(_get_position, _set_position, None,
                       'Position of the fiber array center as a tuple (x, y).')

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

    def check(self):
        '''
        Check if the configuration has errors and raise exceptions if so.
        '''
        if self._spacing < self.fiber.dcore:
            raise ValueError('Spacing between the optical fibers is smaller '
                             'than the diameter of the fiber core!')
        return True

    def fiber_position(self, index: int) -> Tuple[float, float]:
        '''
        Returns the position of the fiber center as a tuple (x, y).

        Parameters
        ----------
        index: int
            Fiber index from 0 to n - 1.

        Returns
        -------
        position: (float, float)
            The position of the fiber center as a tuple (x, y).
        '''
        if index >= self._n or index < -self._n:
            raise IndexError('The fiber index is out of valid range!')
        left = self._position - self._orientation*self._spacing*(self._n - 1)*0.5 
        return tuple(left + self._spacing*self._orientation*int(index))

    def _get_normalized(self) -> np.ndarray:
        return self.raw*(1.0/self.nphotons)
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def cl_pack(self, mc: mcobject.McObject,
                target : cltypes.Structure = None) -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`LinearArray.cl_type` method for a detailed list
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

        target.core_r_squared = 0.25*self.fiber.dcore**2

        target.cos_min = (1.0 - (self._fiber.na/self._fiber.ncore)**2)**0.5

        target.first_position.fromarray(self.fiber_position(0))
        target.delta_position.fromarray(self._orientation*self._spacing)

        return target

    def todict(self):
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`LinearArray.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'LinearArray',
            'fiber': self._fiber.todict(),
            'n': self._n,
            'orientation': self._orientation,
            'spacing': self._spacing,
            'position':self._position.tolist(),
            'direction':self.direction.todict(),
        }

    @staticmethod
    def fromdict(data):
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`LinearArray.todict` method.
        '''
        detector_type = data.pop('type')
        if detector_type != 'LinearArray':
            raise TypeError(
                'Expected a "LinearArray" type bot got "{}"!'.format(
                    detector_type))
        fiber = data.pop('fiber')
        fiber = MultimodeFiber.fromdict(fiber)

        return LinearArray(fiber, **data)

    def __str__(self):
        return 'LinearArray(fiber={}, n={}, spacing={}, orientation=({}, {}), '\
               'position=({}, {}), direction=({}, {}, {}))'.format(
                   self._fiber, self._n, self._spacing, *self._orientation,
                   *self._position, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
