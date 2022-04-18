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

from xopto.mcvox import mcobject
from xopto.mcvox import cltypes
from xopto.mcvox.mcutil.fiber import MultimodeFiber
from xopto.mcvox.mcutil import boundary, geometry

from ..base import SurfaceLayoutAny, TOP, BOTTOM


class LinearArray(SurfaceLayoutAny):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClLinearArray(cltypes.Structure):
            '''
            Structure that that represents a surface layout in the Monte Carlo
            simulator core.

            Fields
            ------
            transformation: mc_matrix3f_t
                Transforms coordinates from Monte Carlo to the surface layout.
            cutout_transformation: mc_matrix3f_t
                Transforms for the cutout to match the fiber array orientation.
            position: mc_point2f_t
                Position of the center of the fiber array.
            first_position: mc_point2f_t
                The center of the first optical fiber in the array.
            delta_position: mc_point2f_t
                Distance vector between two neighboring optical fibers.
            core_r_squared: mc_fp_t
                Squared radius of the optical fiber core.
            core_n: mc_fp_t
                Refractive index of the optical fiber core.
            cladding_r_squared: mc_fp_t
                Squared radius of the optical fiber cladding.
            cladding_n: mc_fp_t
                Refractive index of the optical fiber cladding.
            cutout_width_half: mc_fp_t
                half of the cutout width.
            cutout_height_half: mc_fp_t
                Half of the cutout height.
            cutout_n: mc_fp_t
                Refractive index of the probe cutout.
            probe_r_squared: mc_fp_t
                Squared radius of the optical fiber probe.
            probe_reflectivity: mc_fp_t
                Reflectivity of the probe stainless steel surface.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('cutout_transformation', T.mc_matrix2f_t),
                ('position', T.mc_point2f_t),
                ('first_position', T.mc_point2f_t),
                ('delta_position', T.mc_point2f_t),
                ('core_spacing', T.mc_fp_t),
                ('cladding_r_squared', T.mc_fp_t),
                ('cladding_n', T.mc_fp_t),
                ('core_r_squared', T.mc_fp_t),
                ('core_n', T.mc_fp_t),
                ('cutout_width_half', T.mc_fp_t),
                ('cutout_height_half', T.mc_fp_t),
                ('cutout_n', T.mc_fp_t),
                ('probe_r_squared', T.mc_fp_t),
                ('probe_reflectivity', T.mc_fp_t),
            ]
        return ClLinearArray

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the surface layout in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}SurfaceLayout{{'.format(Loc),
            '	mc_matrix3f_t transformation;'
            '	mc_matrix2f_t cutout_transformation;'
            '	mc_point2f_t position;'
            '	mc_point2f_t first_position;'
            '	mc_point2f_t delta_position;'
            '	mc_fp_t core_spacing;',
            '	mc_fp_t cladding_r_squared;',
            '	mc_fp_t cladding_n;',
            '	mc_fp_t core_r_squared;',
            '	mc_fp_t core_n;',
            '	mc_fp_t cutout_width_half;',
            '	mc_fp_t cutout_height_half;',
            '	mc_fp_t cutout_n;',
            '	mc_fp_t probe_r_squared;',
            '	mc_fp_t probe_reflectivity;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the surface layout in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_surface_layout('.format(loc),
            '		__mc_surface_mem const Mc{}SurfaceLayout *layout){{'.format(Loc),
            '	dbg_print("Mc{}SurfaceLayout - LinearArray surface layout:");'.format(Loc),
            '	dbg_print_matrix3f(INDENT "transformation:", &layout->transformation);',
            '	dbg_print_matrix2f(INDENT "cutout_transformation:", &layout->cutout_transformation);',
            '	dbg_print_point2f(INDENT "position:", &layout->position);',
            '	dbg_print_point2f(INDENT "first_position:", &layout->first_position);',
            '	dbg_print_point2f(INDENT "delta_position:", &layout->delta_position);',
            '',
            '	dbg_print_float(INDENT "core_r_squared (mm2):", layout->core_r_squared*1e6f);',
            '	dbg_print_float(INDENT "core_n:", layout->core_n);',
            '',
            '	dbg_print_float(INDENT "cladding_r_squared (mm2):", layout->cladding_r_squared*1e6f);',
            '	dbg_print_float(INDENT "cladding_n:", layout->cladding_n);',
            '',
            '	dbg_print_float(INDENT "cutout_width_half (mm):", layout->cutout_width_half*1e3f);',
            '	dbg_print_float(INDENT "cutout_height_half (mm):", layout->cutout_height_half*1e3f);',
            '	dbg_print_float(INDENT "cutout_n:", layout->cutout_n);',
            '',
            '	dbg_print_float(INDENT "probe_r_squared (mm2):", layout->probe_r_squared*1e6f);',
            '	dbg_print_float(INDENT "probe_reflectivity", layout->probe_reflectivity);',
            '};',
            '',
            'inline int mcsim_{}_surface_layout_handler('.format(loc),
            '		McSim *mcsim, mc_fp_t *n2){',
            '',
            '	dbg_print_status(mcsim, "{} LinearArray fiber array layout hit");'.format(Loc),
            '',
            '	__mc_surface_mem const Mc{}SurfaceLayout *layout = '.format(Loc),
            '		mcsim_{}_surface_layout(mcsim);'.format(loc),
            '',
            '	mc_fp_t dx, dy, r_squared;',
            '	mc_point3f_t mc_pos, layout_pos;',
            '',
            '	mc_fp_t fiber_x = layout->first_position.x;',
            '	mc_fp_t fiber_y = layout->first_position.y;',
            '',
            '	pragma_unroll_hint({})'.format(self._n),
            '	for(mc_size_t index=0; index < {}; ++index){{'.format(self._n),
            '		mc_pos.x = mcsim_position_x(mcsim) - fiber_x;',
            '		mc_pos.y = mcsim_position_y(mcsim) - fiber_y;',
            '		mc_pos.z = FP_0;',
            '',
            '		mc_matrix3f_t transformation = layout->transformation;',
            '		transform_point3f(&transformation, &mc_pos, &layout_pos);',
            '		dx = layout_pos.x;',
            '		dy = layout_pos.y;',
            '		r_squared = dx*dx + dy*dy;',
            '		if(r_squared <= layout->cladding_r_squared){',
            '			if(r_squared <= layout->core_r_squared){',
            '				/* hit the fiber core */',
            '				*n2 = layout->core_n;',
            '				dbg_print_size_t("{} LinearArray layout fiber core hit:", index);'.format(Loc),
            '				return MC_SURFACE_LAYOUT_CONTINUE;',
            '			};',
            '			*n2 = layout->cladding_n;',
            '			dbg_print_size_t("{} LinearArray layout fiber cladding hit:", index);'.format(Loc),
            '			return MC_SURFACE_LAYOUT_CONTINUE;',
            '		};',
            '',
            '		fiber_x += layout->delta_position.x;',
            '		fiber_y += layout->delta_position.y;',
            '	};',
            '',
            '	/* The cutout of the stainless steel probe. */',
            '	mc_pos.x = mcsim_position_x(mcsim) - layout->position.x;',
            '	mc_pos.y = mcsim_position_y(mcsim) - layout->position.y;',
            '	transform_point2f(&layout->cutout_transformation, &mc_pos, &layout_pos);',
            '	dx = mc_fabs(layout_pos.x);',
            '	dy = mc_fabs(layout_pos.y);',
            '	if(dx <= layout->cutout_width_half && dy < layout->cutout_height_half){',
            '		*n2 = layout->cutout_n;',
            '		dbg_print("{} LinearArray layout cutout hit");'.format(Loc),
            '		return MC_SURFACE_LAYOUT_CONTINUE;',
            '	};',
            '',
            '	/* The tip of the stainless steel probe. */',
            '	dx = mc_pos.x;',
            '	dy = mc_pos.y;',
            '	r_squared = dx*dx + dy*dy;'
            '	if(r_squared <= layout->probe_r_squared){',
            '		mcsim_reverse_direction_z(mcsim);',
            '		mcsim_set_weight(',
            '			mcsim, mcsim_weight(mcsim)*layout->probe_reflectivity);',
            '		dbg_print("{} LinearArray layout stainles steel hit");'.format(Loc),
            '		return MC_REFLECTED;',
            '	};',
            '',
            '	dbg_print("{} LinearArray layout missed");'.format(Loc),
            '	return MC_SURFACE_LAYOUT_CONTINUE;',
            '};',
        ))

    def __init__(self, fiber: MultimodeFiber, n=1, spacing: float = None,
                 orientation: Tuple[float, float] = (1.0, 0.0),
                 diameter: float = 0.0, reflectivity: float = 0.0,
                 cutout: Tuple[float, float] = (0.0, 0.0), cutoutn=1.0,
                 position: Tuple[float, float] = (0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Optical fiber probe layout for a linear array of optical fibers that
        are optionally tilted (direction parameter). The optical fibers are
        always polished in a way that forms a tight optical contact with
        the surface of the sample.

        Parameters
        ----------
        fiber: MultimodeFiber
            Optical properties of the fibers. All optical fibers of the array
            are considered the same.
        n: int
            The number of optical fibers in the linear array. This option
            is a compile-time feature and caanot be changed once the
            surface layout object is created.
        spacing: float
            Spacing between the optical fibers. If spacing is None,
            a tight layout is used with spacing set to the outer diameter
            of the fiber cladding.
        diameter: float
            Outer diameter of the optical fiber probe. Set to 0 if unused.
        reflectivity: float
            Reflectivity of the optical fiber probe stainless steeel surface.
        cutout: (float, float)
            Cutout size as a tuple (width, height). Setting any of the
            cutout width or height to 0 turns off the cutout.
            The cutout is rotated to match the orientation of the linear
            fiber array.
        cutoutn: float
            Refractive index of the cutout. Used only if the cutout surface is
            nonzero.
        orientation: (float, float)
            Vector that points in the direction of the linear fiber array.
            By default the fibers are place in the direction of x axis,
            i.e. vector (1.0, 0.0).
            The orientation must point in the direction from the first to
            the last optical fiber!
        position: (float, float)
            Position of the center of the linear fiber array and of the
            opticalfiber probe probe.
        direction: (float, float, float)
            Reference direction / orientation of the optical fibers. Fibers
            are oriented in this direction and polished to form a tight optical
            contact with the sample (the fiber cross sections are ellipsoids
            if the direction is not perpendicular, i.e different
            from (0, 0, 1).
        '''
        super().__init__()

        if isinstance(fiber, LinearArray):
            la = fiber
            fiber = la.fiber
            n = la.n
            spacing = la.spacing
            orientation = la.orientation
            diameter = la.diameter
            reflectivity = la.reflectivity
            cutout = la.cutout
            cutoutn = la.cutoutn
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

        self._n = 0
        self._fiber = None
        self._spacing = 0.0
        self._diameter = 0.0
        self._reflectivity = 1.0
        self._cutout = np.zeros((2,))
        self._cutout_n = 0.0
        self._orientation = np.array((1.0, 0.0))
        self._direction = np.array((0.0, 0.0, 1.0))
        self._position = np.array((0.0, 0.0))

        self._set_fiber(fiber)
        self._n = max(int(n), 1)
        self._set_spacing(spacing)
        self._set_diameter(diameter)
        self._set_reflectivity(reflectivity)
        self._set_cutout(cutout)
        self._set_cutout_n(cutoutn)
        self._set_orientation(orientation)
        self._set_position(position)
        self._set_direction(direction)

    def _get_fiber(self) -> Tuple[float, float]:
        return self._fiber
    def _set_fiber(self, value: float or Tuple[float, float]):
        self._fiber = value
    fiber = property(_get_fiber, _set_fiber, None,
                     'Properties of the optical fibers used by the layout.')

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

    def _get_diameter(self) -> float:
        return self._diameter
    def _set_diameter(self, diameter: float):
        self._diameter = max(float(diameter), 0.0)
    diameter = property(_get_diameter, _set_diameter, None,
                       'Outer diameter of the optical fiber probe tip.')

    def _get_reflectivity(self) -> float:
        return self._reflectivity
    def _set_reflectivity(self, reflectivity: float):
        self._reflectivity = min(max(float(reflectivity), 0.0), 1.0)
    reflectivity = property(_get_reflectivity, _set_reflectivity, None,
                           'Reflectivity of the stainless steel optical fiber '
                           'probe tip.')

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

    def _get_cutout(self) -> Tuple[float, float]:
        return self._scutout
    def _set_cutout(self, cutout: Tuple[float, float]):
        self._cutout[:] = np.maximum(0.0, cutout)
    cutout = property(_get_cutout, _set_cutout, None,
                       'Size of the cutout (width, height) that accommodates '
                       'the optical fibers. Set to (0, 0) if unused.')

    def _get_cutout_n(self) -> float:
        return self._cutout_n
    def _set_cutout_n(self, n: float):
        self._cutout_n = max(float(n), 1.0)
    cutoutn = property(_get_cutout_n, _set_cutout_n, None,
                       'Refractive index of the cutout fill.')

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
                         'Reference direction of the fibers in the layout.')

    def check(self) -> bool:
        '''
        Check if the configuration has errors and raise exceptions if so.
        '''
        if self._spacing < self.fiber.dcore:
            raise ValueError('Spacing between the optical fibers is smaller '
                             'than the diameter of the fiber core!')
        if self._diameter != 0.0 and \
                self._diameter < self._spacing*(self._n - 1) + \
                self._fiber.dcladding:
            raise ValueError('Diameter of the optical fiber probe too small '
                             'to accommodate the fibers!')
        if np.prod(self._cutout) != 0.0:
            if self.cutout[1] < self._fiber.dcladding or \
                self._cutout[0] < (self._n - 1)*self._spacing + \
                self._fiber.dcladding:
                raise ValueError('The cutout is too small to accommodate the '
                                 'optical fiber array!')
            if self._diameter != 0.0 and \
                    np.linalg.norm(self._cutout) >= self._diameter:
                raise ValueError('The optical fiber probe is too small to '
                                 'accomodate the cutout.')

        return True

    def fiber_position(self, index: int) -> Tuple[float, float]:
        '''
        Returns the position of the fiber center as a tuple (x, y).

        Parameters
        ----------
        index: int
            Fiber index from 0 to n -1.

        Returns
        -------
        position: (float, float)
            The position of the fiber center as a tuple (x, y).
        '''
        if index >= self._n or index < -self._n:
            raise IndexError('The fiber index is out of valid range!')
        left = self._position - self._orientation*self._spacing*(self._n - 1)*0.5 
        return tuple(left + self._spacing*self._orientation*int(index))

    def update(self, other: 'LinearArray' or dict):
        '''
        Update this surface layout configuration from the other surface layout.
        The other surface layout must be of the same type as this surface
        layout or a dict with appropriate fields.

        Parameters
        ----------
        other: LinearArray or dict
            This surface layout is updated with data from this parameter.
        '''
        if isinstance(other, LinearArray):
            self.fiber = other.fiber
            self.spacing = other.spacing
            self.diameter = other.diameter
            self.reflectivity = other.reflectivity
            self.cutout = other.cutout
            self.cutoutn = other.cutoutn
            self.position = other.position
            self.direction = other.direction
        elif isinstance(other, dict):
            self.fiber = other.get('fiber', self.fiber)
            self.spacing = other.get('spacing', self.spacing)
            self.diameter = other.get('diameter', self.diameter)
            self.reflectivity = other.get('reflectivity', self.reflectivity)
            self.cutout = other.get('cutout', self.cutout)
            self.cutoutn = other.get('cutoutn', self.cutoutn)
            self.position = other.get('position', self.position)
            self.direction = other.get('direction', self.direction)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
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

        adir = self._direction[0], self._direction[1], abs(self._direction[2])
        T = geometry.transform_base(adir, (0.0, 0.0, 1.0))

        target.transformation.fromarray(T)
        target.position.fromarray(self._position)

        target.first_position.fromarray(self.fiber_position(0))
        target.delta_position.fromarray(self._orientation*self._spacing)

        target.core_r_squared = 0.25*self._fiber.dcore**2
        target.core_n = self._fiber.ncore

        target.cladding_r_squared = 0.25*self._fiber.dcladding**2
        target.cladding_n = self._fiber.ncladding

        target.probe_r_squared = 0.25*self.diameter**2
        target.probe_reflectivity = self.reflectivity

        T = geometry.rotation_matrix_2d(self._orientation, [1.0, 0.0])
        target.cutout_transformation.fromarray(T)
        target.cutout_width_half = self._cutout[0]*0.5
        target.cutout_height_half = self._cutout[1]*0.5
        target.cutout_n = self._cutout_n

        return target

    def todict(self) -> dict:
        '''
        Save the surface layout configuration to a dictionary.
        Use the :meth:`LinearArray.fromdict` method to create a new
        surface layout instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'LinearArray',
            'fiber': self._fiber.todict(),
            'n': self._n,
            'spacing': self._spacing,
            'orientation': self._orientation.tolist(),
            'diameter': self._diameter,
            'reflectivity': self._reflectivity,
            'cutout': self._cutout.tolist(),
            'cutoutn': self._cutoutn,
            'position':self._position.tolist(),
            'direction':self.direction.tolist(),
        }

    @staticmethod
    def fromdict(data: dict) -> 'LinearArray':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`LinearArray.todict` method.
        '''
        layout_type = data.pop('type')
        if layout_type != 'LinearArray':
            raise TypeError(
                'Expected a "LinearArray" type bot got "{}"!'.format(
                    layout_type))
        fiber = data.pop('fiber')
        fiber = MultimodeFiber.fromdict(fiber)

        return LinearArray(fiber, **data)

    def __str__(self):
        return 'LinearArray(fiber={}, n={}, spacing={}, orientation=({}, {}) '\
               'diameter={}, reflectivity={}, cutout=({}, {}), cutoutn={}, '\
               'position=({}, {}), direction=({}, {}, {}))'.format(
                   self._fiber, self._n, self._spacing, *self._orientation,
                   self._diameter, self._reflectivity,
                   *self._cutout, self._cutoutn
                   *self._position, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
