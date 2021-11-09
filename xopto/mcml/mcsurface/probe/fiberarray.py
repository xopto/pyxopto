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

from typing import Tuple, List

import numpy as np

from xopto.mcml import cltypes, mcobject
from xopto.mcml.mcutil.fiber import MultimodeFiber, FiberLayout
from xopto.mcml.mcutil import boundary, geometry

from ..base import SurfaceLayoutAny, TOP, BOTTOM


class FiberArray(SurfaceLayoutAny):

    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        n = self.n
        class ClFiberArray(cltypes.Structure):
            '''
            Structure that that represents a surface layout in the Monte Carlo
            simulator core.

            Fields
            ------
            transformation: mc_point3f_t*n
                Transforms coordinates from Monte Carlo to the layout / fiber
                coordinates for each fiber.
            fiber_position: mc_point2f_t*n
                Positions of the individual optical fibers.
            cladding_r_squared: mc_fp_t*n
                Squared radius of the optical fiber claddings.
            cladding_n: mc_fp_t*n
                Refractive index of the optical fiber claddings.
            cladding_cos_critical: mc_fp_t*n
                Cosine of the critical angle for refraction sample => fiber cladding.
            core_r_squared: mc_fp_t*n
                Squared radius of the optical fiber cores.
            core_n: mc_fp_t*n
                Refractive index of the optical fiber cores.
            core_cos_critical: mc_fp_t*n
                Cosine of the critical angle for refraction sample => fiber core.
            probe_position: mc_point2f_t
                Position of the probe.
            probe_r_squared: mc_fp_t
                Squared radius of the optical fiber probe tip.
            probe_reflectivity: mc_fp_t
                Reflectivity of the probe stainless steel tip.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t*n),
                ('fiber_position', T.mc_point2f_t*n),
                ('cladding_r_squared', T.mc_fp_t*n),
                ('cladding_n', T.mc_fp_t*n),
                ('cladding_cos_critical', T.mc_fp_t*n),
                ('core_r_squared', T.mc_fp_t*n),
                ('core_n', T.mc_fp_t*n),
                ('core_cos_critical', T.mc_fp_t*n),
                ('probe_position', T.mc_point2f_t),
                ('probe_r_squared', T.mc_fp_t),
                ('reflectivity', T.mc_fp_t),
            ]
        return ClFiberArray

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the surface layout in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        n = self.n
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}SurfaceLayout{{'.format(Loc),
            '	mc_matrix3f_t transformation[{}];'.format(n),
            '	mc_point2f_t fiber_position[{}];'.format(n),
            '	mc_fp_t cladding_r_squared[{}];'.format(n),
            '	mc_fp_t cladding_n[{}];'.format(n),
            '	mc_fp_t cladding_cos_critical[{}];'.format(n),
            '	mc_fp_t core_r_squared[{}];'.format(n),
            '	mc_fp_t core_n[{}];'.format(n),
            '	mc_fp_t core_cos_critical[{}];'.format(n),
            '	mc_point2f_t probe_position;',
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
        n = self.n
        return '\n'.join((
            'void dbg_print_{}_surface_layout('.format(loc),
            '		__mc_surface_mem const Mc{}SurfaceLayout *layout){{'.format(Loc),
            '	dbg_print("Mc{}SurfaceLayout - FiberArray surface layout:");'.format(Loc),
            '	for(mc_size_t index=0; index < {}; ++index){{'.format(n),
            '		dbg_print_size_t(INDENT "Fiber index:", index);',
            '		dbg_print_matrix3f(INDENT INDENT "transformation:", &layout->transformation[index]);',
            '		dbg_print_point2f(INDENT INDENT "fiber_position:", &layout->fiber_position[index]);',
            '		dbg_print_float(INDENT INDENT "cladding_r_squared (mm2):", layout->cladding_r_squared[index]*1e6f);',
            '		dbg_print_float(INDENT INDENT "cladding_n:", layout->cladding_n[index]);',
            '		dbg_print_float(INDENT INDENT "cladding_cos_critical:", layout->cladding_cos_critical[index]);',
            '		dbg_print_float(INDENT INDENT "core_r_squared (mm2):", layout->core_r_squared[index]*1e6f);',
            '		dbg_print_float(INDENT INDENT "core_n:", layout->core_n[index]);',
            '		dbg_print_float(INDENT INDENT "core_cos_critical:", layout->core_cos_critical[index]);',
            '	};',
            '',
            '	dbg_print_point2f(INDENT "probe_position:", &layout->probe_position);',
            '	dbg_print_float(INDENT "probe_r_squared:", layout->probe_r_squared);',
            '	dbg_print_float(INDENT "probe_reflectivity:", layout->probe_reflectivity);',
            '};',
            '',
            'inline int mcsim_{}_surface_layout_handler('.format(loc),
            '		McSim *mcsim, mc_fp_t *n2, mc_fp_t *cc){',
            '',
            '	dbg_print_status(mcsim, "{} FiberArray fiber array layout hit");'.format(Loc),
            '',
            '	__mc_surface_mem const Mc{}SurfaceLayout *layout = '.format(Loc),
            '		mcsim_{}_surface_layout(mcsim);'.format(loc),
            '',
            '	mc_fp_t dx, dy, r_squared;',
            '	mc_point3f_t mc_pos, layout_pos;',
            '',
            '	pragma_unroll_hint({})'.format(n),
            '	for(mc_size_t index=0; index < {}; ++index){{'.format(n),
            '		mc_pos.x = mcsim_position_x(mcsim) - layout->fiber_position[index].x;',
            '		mc_pos.y = mcsim_position_y(mcsim) - layout->fiber_position[index].y;',
            '		mc_pos.z = FP_0;',
            '',
            '		transform_point3f(&layout->transformation[index], &mc_pos, &layout_pos);',
            '		dx = layout_pos.x;',
            '		dy = layout_pos.y;',
            '		r_squared = dx*dx + dy*dy;',
            '		if(r_squared <= layout->cladding_r_squared[index]){',
            '			if(r_squared <= layout->core_r_squared[index]){',
            '				/* hit the fiber core */',
            '				*n2 = layout->core_n[index];',
            '				*cc = layout->core_cos_critical[index];',
            '				dbg_print_size_t("{} FiberArray layout fiber core hit:", index);'.format(Loc),
            '				return MC_SURFACE_LAYOUT_CONTINUE;',
            '			};',
            '			*n2 = layout->cladding_n[index];',
            '			*cc = layout->cladding_cos_critical[index];',
            '			dbg_print_size_t("{} FiberArray layout fiber cladding hit:", index);'.format(Loc),
            '			return MC_SURFACE_LAYOUT_CONTINUE;',
            '		};',
            '	};',
            '',
            '	/* The tip of the stainless steel probe. */',
            '	dx = mcsim_position_x(mcsim) - layout->probe_position.x;',
            '	dy = mcsim_position_y(mcsim) - layout->probe_position.y;',
            '	r_squared = dx*dx + dy*dy;'
            '	if(r_squared <= layout->probe_r_squared){',
            '		mcsim_reverse_direction_z(mcsim);',
            '		mcsim_set_weight(',
            '			mcsim, mcsim_weight(mcsim)*layout->probe_reflectivity);',
            '		dbg_print("{} FiberArray layout stainles steel hit");'.format(Loc),
            '		return MC_REFLECTED;',
            '	};',
            '',
            '	dbg_print("{} FiberArray layout missed");'.format(Loc),
            '	return MC_SURFACE_LAYOUT_CONTINUE;',
            '};',
        ))

    def __init__(self, fibers: List[FiberLayout],
                 diameter: float = 0.0, reflectivity: float = 0.0,
                 position: Tuple[float, float] = (0.0, 0.0)):
        '''
        Optical fiber probe layout for an array of optical fibers that
        are optionally tilted (direction parameter). The optical fibers are
        always polished in a way that forms a tight optical contact with
        the surface of the sample.

        Parameters
        ----------
        fibers: List[FiberLayout]
            A list of optical fiber layouts.
        diameter: float
            Outer diameter of the optical fiber probe. Set to 0 if unused.
        reflectivity: float
            Reflectivity of the optical fiber probe stainless steeel surface.
        position: (float, float)
            Position of the optical probe.
        '''
        super().__init__()

        if isinstance(fibers, FiberArray):
            la = fibers
            fibers = la.fibers
            diameter = la.diameter
            reflectivity = la.reflectivity
            position = la.position

        self._fibers = fibers
        self._diameter = 0.0
        self._reflectivity = 1.0
        self._position = np.array((0.0, 0.0))

        self._set_diameter(diameter)
        self._set_reflectivity(reflectivity)
        self._set_position(position)

    def _get_fibers(self) -> List[FiberLayout]:
        return self._fibers
    def _set_fibers(self, fibers: List[FiberLayout]):
        if len(self._fibers) != len(fibers):
            raise ValueError('The number of optical fibers must not change!')
        if type(self._fibers[0]) != type(fibers[0]):
            raise TypeError('The type of optical fibers must not change!')
        self._fibers[:] = fibers
    fibers = property(_get_fibers, _set_fibers, None,
                     'List of optical fibers.')

    def _get_n_fiber(self) -> int:
        return len(self._fibers)
    n = property(_get_n_fiber, None, None,
                 'Number of optical fiber in the array.')

    def __len__(self):
        return len(self._fibers)

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

    def _get_position(self) -> Tuple[float, float]:
        return self._position
    def _set_position(self, value: float or Tuple[float, float]):
        self._position[:] = value
    position = property(_get_position, _set_position, None,
                       'Position of the fiber array center as a tuple (x, y).')

    def check(self) -> bool:
        '''
        Check if the configuration has errors and raise exceptions if so.
        '''
        for fiber_cfg in self._fibers:
            for other_cfg in self._fibers:
                if other_cfg != fiber_cfg:
                    d = np.linalg.norm(fiber_cfg.position - other_cfg.position)
                    if d < max(fiber_cfg.fiber.dcladding, other_cfg.fiber.dcladding):
                        raise ValueError('Some of the fibers in the '
                                         'layout overlap!')

        if self._diameter != 0.0:
            for fiber_cfg in self._fibers:
                d = np.linalg.norm(fiber_cfg.position - self._position)
                if d + fiber_cfg.fire.dcladding*0.5 > self._diameter:
                    raise ValueError('The probe diameter is too small '
                                     'to accommodate the fibers!')

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
         
        return tuple(self._fibers[index].position[:2])

    def update(self, other: 'FiberArray' or dict):
        '''
        Update this surface layout configuration from the other surface layout.
        The other surface layout must be of the same type as this surface
        layout or a dict with appropriate fields.

        Parameters
        ----------
        other: FiberArray or dict
            This surface layout is updated with data from this parameter.
        '''
        if isinstance(other, FiberArray):
            self.fibers = other.fibers
            self.diameter = other.diameter
            self.reflectivity = other.reflectivity
            self.position = other.position
        elif isinstance(other, dict):
            self.fiber = other.get('fiber', self.fiber)
            self.diameter = other.get('diameter', self.diameter)
            self.reflectivity = other.get('reflectivity', self.reflectivity)
            self.position = other.get('position', self.position)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`FiberArray.cl_type` method for a detailed list
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

        if self.location == TOP:
            n_sample = mc.layers[1].n
        else:
            n_sample = mc.layers[-2].n

        for index, fiber_cfg in enumerate(self._fibers):
            adir = fiber_cfg.direction[0], fiber_cfg.direction[1], \
                   abs(fiber_cfg.direction[2])
            T = geometry.transform_base(adir, (0.0, 0.0, 1.0))

            target.transformation[index].fromarray(T)
            target.fiber_position[index].fromarray(fiber_cfg.position)

            target.cladding_r_squared[index] = 0.25*fiber_cfg.fiber.dcladding**2
            target.cladding_n[index] = fiber_cfg.fiber.ncladding
            target.cladding_cos_critical[index] = boundary.cos_critical(
                n_sample, fiber_cfg.fiber.ncladding)

            target.core_r_squared[index] = 0.25*fiber_cfg.fiber.dcore**2
            target.core_n[index] = fiber_cfg.fiber.ncore
            target.core_cos_critical[index] = boundary.cos_critical(
                n_sample, fiber_cfg.fiber.ncore)

        target.probe_position.fromarray(self._position)
        target.probe_r_squared = 0.25*self._diameter**2
        target.probe_reflectivity = self._reflectivity

        return target

    def todict(self) -> dict:
        '''
        Save the surface layout configuration to a dictionary.
        Use the :meth:`FiberArray.fromdict` method to create a new
        surface layout instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'FiberArray',
            'fibers': [item.todict() for item in self._fibers],
            'diameter': self._diameter,
            'reflectivity': self._reflectivity,
            'position':self._position.tolist(),
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'FiberArray':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`FiberArray.todict` method.
        '''
        layout_type = data.pop('type')
        if layout_type != cls.__name__:
            raise TypeError(
                'Expected a "{}" type bot got "{}"!'.format(
                    cls.__name__, layout_type))
        fibers = data.pop('fibers')
        fibers = [FiberLayout.fromdict(item) for item in fibers]

        return cls(fibers, **data)

    def __str__(self):
        return 'FiberArray(fibers={}, diameter={}, reflectivity={}, '\
               'position=({}, {}))'.format(
                   self._fibers, self._diameter, self._reflectivity,
                   *self._position)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
