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

from xopto.mcml import cltypes, mcobject
from xopto.mcml.mcutil import fiber as fiberutil
from xopto.mcml.mcutil import boundary, geometry

from ..base import SurfaceLayoutAny, TOP, BOTTOM


class SixAroundOne(SurfaceLayoutAny):
    @staticmethod
    def cl_type(mc: mcobject.McObject):
        '''
        Surface layout of a tilted six-around one probe with stainless steel
        housing.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents the six-around-one surface layout
            in the Monte Carlo kernel.

        Fields
        ------
        transformation: mc_matrix3f_t
            Transformation that transforms Monte Carlo coordinates relative to
            the probe center to probe coordinates
        position: mc_point3f_t
            Position of the probe center.
        core_spacing: float
            Spacing between the cores of the central and sourrounding optical
            fibers.
        core_r_squared: float
            Squared radius of the optical fiber core.
        core_n: float
            Refractive index of the fiber core.
        core_cos_critical: float
            Critical angle cosine for transition top layer => fiber core.
        cladding_r_squared: float
            Squared radius of the optical fiber cladding.
        cladding_n: float
            Refractive index of the fiber cladding.
        cladding_cos_critical: float
            Critical angle cosine for transition top layer => fiber cladding.
        cutout_r_squared: float
            Squared radius of the cutout in the probe tip that accommodates
            the fibers.
        cutout_n: float
            Refractive index of the probe cutout.
        cutout_cos_critical: float
            Critical angle cosine for transition top layer => cutout.
        probe_r_squared: float
            Squared radius of the probe.
        probe_reflectivity: float
            Reflectivity of the probe tip (a value from 0 to 1 is required).
        '''
        T = mc.types
        class ClSixAroundOne(cltypes.Structure):
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point2f_t),
                ('core_spacing', T.mc_fp_t),
                ('cladding_r_squared', T.mc_fp_t),
                ('cladding_n', T.mc_fp_t),
                ('cladding_cos_critical', T.mc_fp_t),
                ('core_r_squared', T.mc_fp_t),
                ('core_n', T.mc_fp_t),
                ('core_cos_critical', T.mc_fp_t),
                ('cutout_r_squared', T.mc_fp_t),
                ('cutout_n', T.mc_fp_t),
                ('cutout_cos_critical', T.mc_fp_t),
                ('probe_r_squared', T.mc_fp_t),
                ('probe_reflectivity', T.mc_fp_t),
            ]
        return ClSixAroundOne

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the layout in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}SurfaceLayout{{'.format(Loc),
            '	mc_matrix3f_t transformation;'
            '	mc_point2f_t position;'
            '	mc_fp_t core_spacing;',
            '	mc_fp_t cladding_r_squared;',
            '	mc_fp_t cladding_n;',
            '	mc_fp_t cladding_cos_critical;',
            '	mc_fp_t core_r_squared;',
            '	mc_fp_t core_n;',
            '	mc_fp_t core_cos_critical;',
            '	mc_fp_t cutout_r_squared;',
            '	mc_fp_t cutout_n;',
            '	mc_fp_t cutout_cos_critical;',
            '	mc_fp_t probe_r_squared;',
            '	mc_fp_t probe_reflectivity;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_surface_layout('.format(loc),
            '		__mc_surface_mem const Mc{}SurfaceLayout *layout){{'.format(Loc),
            '	dbg_print("Mc{}SurfaceLayout - SixAroundOne surface layout:");'.format(Loc),
            '	dbg_print_matrix3f(INDENT "transformation:", &layout->transformation);',
            '	dbg_print_point2f(INDENT "position:", &layout->position);',
            '	dbg_print_float(INDENT "core_spacing (mm):", layout->core_spacing*1e3f);',
            '',
            '	dbg_print_float(INDENT "core_r_squared (mm2):", layout->core_r_squared*1e6f);',
            '	dbg_print_float(INDENT "core_n:", layout->core_n);',
            '	dbg_print_float(INDENT "core_cos_critical:", layout->core_cos_critical);',
            '',
            '	dbg_print_float(INDENT "cladding_r_squared (mm2):", layout->cladding_r_squared*1e6f);',
            '	dbg_print_float(INDENT "cladding_n:", layout->cladding_n);',
            '	dbg_print_float(INDENT "cladding_cos_critical:", layout->cladding_cos_critical);',
            '',
            '	dbg_print_float(INDENT "cutout_r_squared (mm2):", layout->cutout_r_squared*1e6f);',
            '	dbg_print_float(INDENT "cutout_n:", layout->cutout_n);',
            '	dbg_print_float(INDENT "cutout_cos_critical:", layout->cutout_cos_critical);',
            '',
            '	dbg_print_float(INDENT "probe_r_squared (mm2):", layout->probe_r_squared*1e6f);',
            '	dbg_print_float(INDENT "probe_reflectivity", layout->probe_reflectivity);',
            '};',
            '',
            'inline int mcsim_{}_surface_layout_handler('.format(loc),
            '		McSim *mcsim, mc_fp_t *n2, mc_fp_t *cc){',
            '',
            '	__mc_surface_mem const struct Mc{}SurfaceLayout *layout = '.format(Loc),
            '		mcsim_{}_surface_layout(mcsim);'.format(loc),
            '',
            '	mc_fp_t dx, dy, r_squared;',
            '	mc_point3f_t mc_pos, layout_pos;',
            '',
            '	dbg_print_point3f("{} SixAroundOne layout hit: ", '.format(Loc),
            '		mcsim_direction(mcsim));',
            '',
            '	mc_point3f_t rel_pos = {',
            '		mcsim_position_x(mcsim) - layout->position.x,',
            '		mcsim_position_y(mcsim) - layout->position.y,',
            '		FP_0',
            '	};',
            '',
            '	/* The central fiber */',
            '	mc_pos.x = rel_pos.x;',
            '	mc_pos.y = rel_pos.y;',
            '	mc_pos.z = FP_0;',
            '	transform_point3f(&layout->transformation, &mc_pos, &layout_pos);',
            '	dx = layout_pos.x;',
            '	dy = layout_pos.y;',
            '	r_squared = dx*dx + dy*dy;',
            '	if(r_squared <= layout->cladding_r_squared){',
            '		if(r_squared <= layout->core_r_squared){',
            '			/* hit the fiber core */',
            '			*n2 = layout->core_n;',
            '			*cc = layout->core_cos_critical;',
            '			dbg_print("{} SixAroundOne layout fiber core 0 hit");'.format(Loc),
            '			return MC_SURFACE_LAYOUT_CONTINUE;',
            '		};',
            '		*n2 = layout->cladding_n;',
            '		*cc = layout->cladding_cos_critical;',
            '		dbg_print("{} SixAroundOne layout fiber cladding 0 hit");'.format(Loc),
            '		return MC_SURFACE_LAYOUT_CONTINUE;',
            '	};',
            '',
            '	/* The 1st and 4th fibers of the ring. */',
            '	mc_pos.x = mc_fabs(rel_pos.x) - layout->core_spacing;',
            '	mc_pos.y = rel_pos.y;',
            '	transform_point3f(&layout->transformation, &mc_pos, &layout_pos);',
            '	dx = layout_pos.x;',
            '	dy = layout_pos.y;',
            '	r_squared = dx*dx + dy*dy;',
            '	if(r_squared <= layout->cladding_r_squared){',
            '		if(r_squared <= layout->core_r_squared){',
            '			/* hit the fiber core */',
            '			*n2 = layout->core_n;',
            '			*cc = layout->core_cos_critical;',
            '			dbg_print("{} SixAroundOne layout fiber core 1 or 4 hit");'.format(Loc),
            '			return MC_SURFACE_LAYOUT_CONTINUE;',
            '		};',
            '		*n2 = layout->cladding_n;',
            '		*cc = layout->cladding_cos_critical;',
            '		dbg_print("{} SixAroundOne layout fiber cladding 1 or 4 hit");'.format(Loc),
            '		return MC_SURFACE_LAYOUT_CONTINUE;',
            '	};'
            '',
            '	/* The 2nd, 3rd, 5th and 6th fibers of the ring. */',
            '	mc_pos.x = mc_fabs(rel_pos.x) - layout->core_spacing*FP_0p5;',
            '	mc_pos.y = mc_fabs(rel_pos.y) - layout->core_spacing*FP_COS_30;',
            '	transform_point3f(&layout->transformation, &mc_pos, &layout_pos);',
            '	dx = layout_pos.x;',
            '	dy = layout_pos.y;',
            '	r_squared = dx*dx + dy*dy;',
            '	if(r_squared <= layout->cladding_r_squared){',
            '		if(r_squared <= layout->core_r_squared){',
            '			/* hit the fiber core */',
            '			*n2 = layout->core_n;',
            '			*cc = layout->core_cos_critical;',
            '			dbg_print("{} SixAroundOne layout fiber core 2, 3, 5 or 6 hit");'.format(Loc),
            '			return MC_SURFACE_LAYOUT_CONTINUE;',
            '		};',
            '		*n2 = layout->cladding_n;',
            '		*cc = layout->cladding_cos_critical;',
            '		dbg_print("{} SixAroundOne layout fiber cladding 2, 3, 5 or 6 hit");'.format(Loc),
            '		return MC_SURFACE_LAYOUT_CONTINUE;',
            '	};',
            '',
            '',
            '	/* The cutout of the stainless steel probe. */',
            '	dx = rel_pos.x;',
            '	dy = rel_pos.y;',
            '	r_squared = dx*dx + dy*dy;',
            '	if(r_squared <= layout->cutout_r_squared){',
            '		*n2 = layout->cutout_n;',
            '		*cc = layout->cutout_cos_critical;',
            '		dbg_print("{} SixAroundOne layout cutout hit");'.format(Loc),
            '		return MC_SURFACE_LAYOUT_CONTINUE;',
            '	};',
            '',
            '	/* The tip of the stainless steel probe. */',
            '	if(r_squared <= layout->probe_r_squared){',
            '		mcsim_reverse_direction_z(mcsim);',
            '		mcsim_set_weight(',
            '			mcsim, mcsim_weight(mcsim)*layout->probe_reflectivity);',
            '		dbg_print("{} SixAroundOne layout stainles steel hit");'.format(Loc),
            '		return MC_REFLECTED;',
            '	};',
            '',
            '	dbg_print("{} SixAroundOne layout missed");'.format(Loc),
            '	return MC_SURFACE_LAYOUT_CONTINUE;',
            '};',
        ))
    def __init__(self, fiber: fiberutil.MultimodeFiber or 'SixAroundOne',
                 spacing: float = None,
                 diameter: float = 0.0, reflectivity: float = 1.0,
                 cutout: float = 0.0, cutoutn: float = 1.0,
                 position: Tuple[float, float] = (0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        A six-around one optical fiber probe surface layout with optional
        tilted fibers (direction parameter) and a round cutout region that
        accommodates the optical fibers. The optical fibers are always
        polished in way that forms a tight optical contact with the surface
        of the sample.

        Parameters
        ----------
        fiber: MultimodeFiber or SixAroundOne
            Properties of the optical fibers or an instance of SixAroundOne.
            A copy is made if an instance of SixAroundOne.
        spacing: float
            Spacing between the optical fibers. If spacing is None,
            a tight layout is used with spacing set to the outer diameter
            of the fiber cladding.
        diameter: float
            Outer diameter of the probe tip.
        reflectivity: float
            Reflectivity of the probe tip.
        cutout: float
            Diameter of the cutout accommodating the optical fibers. Set to 0
            if not used.
        cutoutn: float
            Refractive index of the cutout. Only used if cutout is not 0.
        position: (float, float)
            Position of the center of the probe / central finer.
        direction: (float, float, float)
            Reference direction / orientation of the fibers. Fibers are
            oriented in this direction and polished to form a tight optical
            contact with the sample (the fiber cross sections are ellipsoids
            if the direction is not perpendicular, i.e different
            from (0, 0, 1).
        '''
        super().__init__()
        if isinstance(fiber, SixAroundOne):
            sao = fiber
            fiber = sao.fiber
            spacing = sao.spacing
            cutout = sao.cutout
            reflectivity = sao.reflectivity
            diameter = sao.diameter
            position = sao.position
            direction = sao.direction
        else:
            if spacing is None:
                spacing = fiber.dcladding

        self._fiber = fiber
        self._spacing = 0.0
        self._cutout = 0.0
        self._cutout_n = 1.0
        self._reflectivity = 1.0
        self._diameter = 1.0
        self._position = np.zeros((2,))
        self._direction = np.array((0.0, 0.0, 1.0))

        self._set_spacing(float(spacing))
        self._set_diameter(float(diameter))
        self._set_reflectivity(float(reflectivity))
        self._set_cutout(float(cutout))
        self._set_cutout_n(float(cutoutn))
        self._set_position(position) 
        self._set_direction(direction) 

    def _get_fiber(self) -> fiberutil.MultimodeFiber:
        return self._fiber
    def _set_fiber(self, fiber: fiberutil.MultimodeFiber):
        self._fiber = fiber
    fiber = property(_get_fiber, _set_fiber, None,
                     'Multimode optical fiber.')

    def _get_spacing(self) -> float:
        return self._spacing
    def _set_spacing(self, spacing: float):
        self._spacing = max(float(spacing), 0.0)
    spacing = property(_get_spacing, _set_spacing, None,
                       'Spacing of the optical fibers (m).')

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

    def _get_cutout(self) -> float:
        return self._scutout
    def _set_cutout(self, cutout: float):
        self._cutout = max(float(cutout), 0.0)
    cutout = property(_get_cutout, _set_cutout, None,
                       'Cutout diameter (m) accommodating the optical fibers.'
                       'Set to 0 if unused.')

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
                       'Position of the layout as a tuple (x, y).')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('Direction vector norm/length must not be 0!')
        self._direction *= 1.0/norm
    direction = property(_get_direction, _set_direction, None,
                         'Direction of the optical fibers.')

    def check(self) -> bool:
        '''
        Check if the properties of the optical fiber probe are correct and
        raise exception if not.
        '''
        if self._fiber.dcladding > self._spacing:
            raise ValueError('The optical fibers are overlapping!')
        if self._diameter < self._fiber.dcladding:
           raise ValueError('The probe diameter is too small to accommodate '
                            'the fibers')
        if self._cutout != 0.0 and \
                self._cutout < self._spacing + self._fiber.dcladding:
            raise ValueError(
                'The cutout is too small to accommodate the optical fibers!')

        return True

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

    def update(self, other: 'SixAroundOne' or dict):
        '''
        Update this surface layout configuration from the other surface layout.
        The other surface layout must be of the same type as this surface
        layout or a dict with appropriate fields.

        Parameters
        ----------
        other: SixAroundOne or dict
            This surface layout is updated with data from this parameter.
        '''
        if isinstance(other, SixAroundOne):
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
        See the :py:meth:`SixAroundOne.cl_type` method for a detailed list
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

        if self.location == TOP:
            n_sample = mc.layers[1].n
        else:
            n_sample = mc.layers[-2].n

        target.transformation.fromarray(T)
        target.position.fromarray(self._position)
        target.core_spacing = self._spacing

        target.core_r_squared = 0.25*self._fiber.dcore**2
        target.core_n = self._fiber.ncore
        target.core_cos_critical = boundary.cos_critical(
            n_sample, self._fiber.ncore)

        target.cladding_r_squared = 0.25*self._fiber.dcladding**2
        target.cladding_n = self._fiber.ncladding
        target.cladding_cos_critical = boundary.cos_critical(
            n_sample, self._fiber.ncladding)

        target.cutout_r_squared = 0.25*self._cutout**2
        target.cutout_n = self._cutout_n
        target.cutout_cos_critical = boundary.cos_critical(
            n_sample, self._cutout_n)

        target.probe_r_squared = 0.25*self.diameter**2
        target.probe_reflectivity = self.reflectivity

        return target

    def __str__(self):
        return 'SixAroundOne(fiber={}, spacing={}, '\
                'diameter={}, reflectivity={}, '\
                'cutout={}, cutoutn={}, '\
                'position=({}, {}), direction=({}, {}, {}))'.format(
                    self._fiber, self._spacing,
                    self._diameter, self._reflectivity,
                    self._cutout, self._cutout_n,
                    *self._position, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
