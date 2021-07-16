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

from xopto.mcvox.mcsource.base import Source
from xopto.mcvox import mcobject
from xopto.mcvox import mctypes
from xopto.mcvox import mcoptions
from xopto.mcvox.mcutil import boundary
from xopto.mcvox.mcutil import geometry
from xopto.mcvox.mcutil import lut as lututil
from xopto.mcvox.mcutil import fiber as fiberutil
from xopto.mcvox import cltypes

class TopGeometryFiber(mcobject.McObject):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Top sample surface geometry representing an optical fiber.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents a top surface geometry in the
            Monte Carlo kernel.

            The returned structure type implements the following fields:

            - transformation mc_matrix3f_t
                Transformation from the Monte Carlo coordinate system to the
                fiber coordinate system (position must be expressed relative to
                the fiber center).
            - position: mc_point3f_t
                Position of the optical fiber in the Monte Carlo coordinate system.
            - core_radius: mc_fp_t
                Outer radius of the fiber core.
            - cladding_radius: mc_fp_t
                Outer radius of the fiber cladding.
            - core_n: mc_fp_t
                Refractive index of the fiber core.
            - cladding_n: mc_fp_t
                Refractive index of the fiber cladding.
        '''
        T = mc.types
        class ClTopGeometryFiber(cltypes.Structure):
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point3f_t),
                ('core_radius', T.mc_fp_t),
                ('cladding_radius', T.mc_fp_t),
                ('core_n', T.mc_fp_t),
                ('cladding_n', T.mc_fp_t),
            ]
        return ClTopGeometryFiber

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McComplexSurfaceTop{',
            '	mc_matrix3f_t transformation;',
            '	mc_point3f_t position;',
            '	mc_fp_t core_radius;',
            '	mc_fp_t cladding_radius;',
            '	mc_fp_t core_n;',
            '	mc_fp_t cladding_n;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        return '\n'.join((
            'void print_top_geometry(__mc_geometry_mem const McComplexSurfaceTop *geo){',
            '	dbg_print("TopGeometryFiber:");',
            '	dbg_print_matrix3f(INDENT "transformation:", &geo->transformation);',
            '	dbg_print_point3f(INDENT "position:", &geo->position);',
            '	dbg_print_float(INDENT "core_radius (mm):", geo->core_radius*1e3f);',
            '	dbg_print_float(INDENT "cladding_radius (mm):", geo->cladding_radius*1e3f);',
            '	dbg_print_float(INDENT "core_n:", geo->core_n);',
            '	dbg_print_float(INDENT "cladding_n:", geo->cladding_n);',
            '};',
            '',
            'inline mc_int_t mcsim_top_geometry_handler(McSim *psim, mc_fp_t *n2, mc_fp_t *cc){',
            '	__mc_geometry_mem const struct McComplexSurfaceTop *geometry = mcsim_top_geometry(mcsim);',
            '',
            '	mc_point3f_t pt_mc, pt_src;',
            '	mc_fp_t n = mc_material_n(mcsim_current_voxel_material(psim));',
            '	pt_mc.x = mcsim_position_x(psim) - geometry->position.x;',
            '	pt_mc.y = mcsim_position_y(psim) - geometry->position.y;',
            '	pt_mc.z = FP_0;',
            '',
            '	transform_point3f(&geometry->transformation, &pt_mc, &pt_src);',
            '	mc_fp_t r2 = pt_src.x*pt_src.x + pt_src.y*pt_src.y'
            '	if (dr2 <= geometry->cladding_radius*geometry->cladding_radius){',
            '		*n2 = geometry->cladding_n;',
            '		*cc = cos_critical(n, *n2);',
            '		if (dr2 < geometry->core_radius*geometry->core_radius){',
            '			*n2 = geometry->core_n;',
            '			*cc = cos_critical(n, *n2);',
            '		};'
            '	}',
            '	return MC_SURFACE_GEOMETRY_CONTINUE;',
            '};',
        ))

    def __init__(self, fiber: 'UniformFiber' or 'LambertianFiber' or
                              'UniformFiberLut'):
        self._fiber_src = fiber

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`FiberGeometryTop.cl_type` for a detailed
        list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
            Filled ctypes structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        fiber_src = self._fiber_src
        fiber = fiber_src.fiber
        T = geometry.transform_base(fiber_src.direction, (0.0, 0.0, 1.0))

        target.transformation.fromarray(T)
        target.position.fromarray(fiber_src.position)
        target.core_radius = fiber.dcore*0.5
        target.cladding_radius = fiber.dcladding*0.5
        target.core_n = fiber.ncore
        target.cladding_n = fiber.ncladding

        return target

class UniformFiber(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformFiber(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Parameters
            ----------
            mc: McObject
                A Monte Carlo simulator instance.

            Returns
            -------
            struct: cltypes.Structure
                A structure type that represents a uniform optical fiber source
                in the Monte Carlo kernel.
                
                The returned structure type implements the following fields:

                - transformation: mc_matrix3f_t
                    Transformation from the beam coordinate system to the Monte
                    Carlo coordinate system.
                - position: mc_point3f_t
                    Source position (axis).
                - direction: mc_point3f_t
                    Source direction (axis).
                - radius: mc_fp_t
                    Radius of the fiber core.
                - cos_min: mc_fp_t
                    Cosine of the fiber acceptance angle in air.
                - n: mc_fp_t
                    Refractive index of the fiber core.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point3f_t),
                ('direction', T.mc_point3f_t),
                ('radius', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n', T.mc_fp_t),
            ]
        return ClUniformFiber

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
            '	mc_fp_t radius;',
            '	mc_fp_t cos_min;',
            '	mc_fp_t n;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("UniformFiber source:");',
            '	dbg_print_matrix3f(INDENT "transformation:", &src->transformation);',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
             '	dbg_print_point3f(INDENT "direction:", &src->direction);',
            '	dbg_print_float(INDENT "radius (mm):", src->radius*1e3f);',
            '	dbg_print_float(INDENT "cos_min:", src->cos_min);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_point3f_t pt_src, pt_mc;',
            '	mc_fp_t specular_r = FP_0;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t r = mc_sqrt(mcsim_random(mcsim))*source->radius;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	pt_src.x = r*cos_fi;',
            '	pt_src.y = r*sin_fi;',
            '	pt_src.z = FP_0;',
            '',
            '	transform_point3f(&source->transformation, &pt_src, &pt_mc);',
            '	pt_mc.x += source->position.x;',
            '	pt_mc.y += source->position.y;',
            '	pt_mc.z += source->position.z;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	cos_theta = FP_1 - mcsim_random(mcsim)*(FP_1 - source->cos_min);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '	/* adjust the emission angle for the refractive index of the fiber */ ',
            '	sin_theta = mc_fdiv(sin_theta, source->n);',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);',
            '',
            '	pt_src.x = cos_fi*sin_theta;',
            '	pt_src.y = sin_fi*sin_theta;',
            '	pt_src.z = cos_theta;',
            '	mc_point3f_t packet_direction;',
            '	transform_point3f(&source->transformation, &pt_src, &packet_direction);',
            '',
            '	mc_point3f_t source_direction = source->direction;',
            '',
            '	if (mcsim_sample_box_contains(mcsim, &pt_mc)) {',
            '		reverse3f(&source_direction, &source_direction);',
            '		dbg_print("Move to the sample box surface:");',
            '		dbg_print_point3f(INDENT "from        :", &pt_mc);',
            '		dbg_print_point3f(INDENT "in direction:", &source_direction);',
            '	};'
            '',
            '	mc_point3f_t refracted_direction = packet_direction;',
            '	mc_point3f_t position = source->position;',
            '	mc_point3f_t intersection, normal;',
            '	if (mcsim_sample_box_intersect(',
            '			mcsim, &pt_mc, &source_direction, ',
            '			&intersection, &normal)) {',
            '		mcsim_set_position(mcsim, &intersection);',
            '',
            '		reverse3f(&normal, &normal);',
            '',
            '		dbg_print("Surface intersection:");',
            '		dbg_print_point3f(INDENT "position:", &intersection);',
            '		dbg_print_point3f(INDENT "normal:  ", &normal);',
            '',
            '		position = intersection;',
            '',
            '		mc_fp_t n_outside = source->n;',
            '',
            '		mc_point3_t voxel_index;',
            '		mcsim_position_to_voxel_index_safe(',
            '			mcsim, &intersection, &voxel_index);',
            '		dbg_print_point3("Launch voxel:", &voxel_index);',
            '		mc_size_t material_index = mcsim_voxel_material_index(',
            '			mcsim, &voxel_index);',
            '		dbg_print_size_t("Launch material:", material_index);',
            '		mc_fp_t n_inside = mc_material_n(',
            '			mcsim_material(mcsim, material_index));',
            '',
            '		refract_safe(&packet_direction, &normal, ',
            '			n_outside, n_inside, &refracted_direction);',
            '',
            '		specular_r = reflectance(n_outside, n_inside, ',
            '				dot3f(&normal, &packet_direction),',
            '				cos_critical(n_outside, n_inside));',
            '',
            '		dbg_print_status(mcsim, "Moved to the sample surface");',
            '	} else {',
            '		position = *mcsim_top_left(mcsim);',
            '		specular_r = FP_1;',
            '		dbg_print_status(mcsim, "Launch missed the sample box!");',
            '	}',
            '',
            '	mcsim_set_position(mcsim, &position);',
            '	mcsim_set_direction(mcsim, &refracted_direction)',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, &position, &packet_direction, specular_r);',
            '	#endif',
            '',
            '	mcsim_set_weight(mcsim, FP_1 - specular_r);',
            '	dbg_print_status(mcsim, "Launch UniformFiber");',
            '};',
        ))

    def __init__(self, fiber: fiberutil.MultimodeFiber,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        An optical fiber photon packet source.

        Parameters
        ----------
        fiber: MultimodeFiber
            Optical fiber parameters.
        position: (float, float, float)
            Position of the center of the fiber source array-like object as
            tuple (x, y, z). The z coordinate is ignored and set to 0.
        direction: (float, float, float)
            Direction of the fiber axis as an array-like object of size 3
            (px, py, pz). If not perpendicular to the sample surface, the
            fiber tip is cut at an angle so that the fiber surface is
            parallel with the sample surface at the bem entry point. The z
            component of the direction vector must be positive (fiber pointing
            towards the top surface of the sample) in order to use the top
            surface layout :py:class:`TopGeometryFiber`.

        Note
        ----
        The fiber will be always terminated in a way that forms a tight coupling
        between the sample surface and the fiber tip. If the incidence is
        not normal, the fiber will have an elliptical cross-section (
        cut at an angle).
        The entry point on the sample surface will be determined by propagating
        the position along the given direction (no interactions with the medium
        during this step). If the photon packet does not intersect the sample,
        the initial weight will be set to 0 (reflectance to 1) and the packet
        will be launched from the top-left corner of the voxelized sample box.
        Such zero-weight packets are immediately terminated and have no
        contribution to the fluence and surface detectors, however will be
        included in the trace (have no effect on the sampling volume or other
        trace-based analysis due to zero-weight).
        Note that in case the position lies within the sample, the
        position will be propagated to the entry point using reversed direction.
        From there the packets will be launched according to the NA of the
        fiber and refractive index of the sample surface. The MC simulation
        will start after subtracting the specular reflectance at the
        boundary from the initial weight of the packet.
        '''
        Source.__init__(self)

        self._fiber = fiber
        self._position = np.zeros((3,))
        self._direction = np.zeros((3,))
        self._direction[2] = 1.0

        self._set_position(position)
        self._set_direction(direction)

    def update(self, other: 'UniformFiber' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: UniformFiber or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, 'UniformFiber'):
            self.fiber = other.fiber
            self.position = other.position
            self.direction = other.direction
        elif isinstance(other, dict):
            self.fiber = other.get('fiber', self.fiber)
            self.position = other.get('position', self.position)
            self.direction = other.get('direction', self.direction)

    def _get_fiber(self) -> fiberutil.MultimodeFiber:
        return self._fiber
    def _set_fiber(self, fib: fiberutil.MultimodeFiber):
        self._fiber = fiberutil.MultimodeFiber(fib)
    fiber = property(_get_fiber, _set_fiber, None,
                     'Multimode optical fiber.')

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, position: Tuple[float, float, float]):
        self._position[:] = position
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        dir_norm = np.linalg.norm(direction)
        if dir_norm == 0.0:
            raise ValueError('The direction vector is singular!')
        self._direction[:] = direction
        self._direction *= 1.0/dir_norm
    direction = property(_get_direction, _set_direction, None,
                        'Source direction.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`UniformFiber.cl_type` for a detailed
        list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
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

        T = geometry.transform_base((0.0, 0.0, 1.0), self._direction)

        target.transformation.fromarray(T)

        target.n = self._fiber.ncore
        target.cos_min = (1.0 - (self._fiber.na)**2)**0.5

        target.position.fromarray(self.position)
        target.direction.fromarray(self.direction)

        target.radius = self._fiber.dcore*0.5

        return target, None, None

    def __str__(self):
        return 'UniformFiber(fiber={}, '\
               'position=({}, {}, {}), direction=({}, {}, {}))'.format(
                   self._fiber, *self._position, *self._direction)


class LambertianFiber(UniformFiber):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClLambertianFiber(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Parameters
            ----------
            mc: McObject
                A Monte Carlo simulator instance.

            Returns
            -------
            struct: cltypes.Structure
                A structure type that represents a lambertian fiber source in
                the Monte Carlo kernel.

                The returned structure type implements the following fields:

                - transformation: mc_matrix3f_t
                    Transformation from the beam coordinate system to the Monte
                    Carlo coordinate system.
                - position: mc_point3f_t
                    Source position (axis).
                - direction: mc_point3f_t
                    Source direction (axis).
                - radius: mc_fp_t
                    Radius of the fiber core.
                - na: mc_fp_t
                    Numerical aperture of the fiber core - sine of the acceptance
                    angle in air.
                - n: mc_fp_t
                    Refractive index of the fiber core.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point3f_t),
                ('direction', T.mc_point3f_t),
                ('radius', T.mc_fp_t),
                ('na', T.mc_fp_t),
                ('n', T.mc_fp_t),
            ]
        return ClLambertianFiber

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
            '	mc_fp_t radius;',
            '	mc_fp_t na;',
            '	mc_fp_t n;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("LambertianFiber source:");',
            '	dbg_print_matrix3f(INDENT "transformation:", &src->transformation);',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
             '	dbg_print_point3f(INDENT "direction:", &src->direction);',
            '	dbg_print_float(INDENT "radius (mm):", src->radius*1e3f);',
            '	dbg_print_float(INDENT "na:", src->na);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_point3f_t pt_src, pt_mc;',
            '	mc_fp_t specular_r = FP_0;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t r = mc_sqrt(mcsim_random(mcsim))*source->radius;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	pt_src.x = r*cos_fi;',
            '	pt_src.y = r*sin_fi;',
            '	pt_src.z = FP_0;',
            '',
            '	transform_point3f(&source->transformation, &pt_src, &pt_mc);',
            '	pt_mc.x += source->position.x;',
            '	pt_mc.y += source->position.y;',
            '	pt_mc.z += source->position.z;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	sin_theta = mc_sqrt(mcsim_random(mcsim))*source->na;',
            '	/* Adjust the propagation direction to fiber refractive index */',
            '	sin_theta = mc_fdiv(sin_theta, source->n);',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);',
            '',
            '	pt_src.x = cos_fi*sin_theta;',
            '	pt_src.y = sin_fi*sin_theta;',
            '	pt_src.z = cos_theta;',
            '	mc_point3f_t packet_direction;',
            '	transform_point3f(&source->transformation, &pt_src, &packet_direction);',
            '',
            '	 mc_point3f_t source_direction = source->direction;',
            '',
            '	if (mcsim_sample_box_contains(mcsim, &pt_mc)) {',
            '		reverse3f(&source_direction, &source_direction);',
            '		dbg_print("Move to the sample box surface:");',
            '		dbg_print_point3f(INDENT "from        :", &pt_mc);',
            '		dbg_print_point3f(INDENT "in direction:", &source_direction);',
            '	};'
            '',
            '	mc_point3f_t refracted_direction = packet_direction;',
            '	mc_point3f_t position = source->position;',
            '	mc_point3f_t intersection, normal;',
            '	if (mcsim_sample_box_intersect(',
            '			mcsim, &pt_mc, &source_direction, ',
            '			&intersection, &normal)) {',
            '		mcsim_set_position(mcsim, &intersection);',
            '',
            '		reverse3f(&normal, &normal);',
            '',
            '		dbg_print("Surface intersection:");',
            '		dbg_print_point3f(INDENT "position:", &intersection);',
            '		dbg_print_point3f(INDENT "normal:  ", &normal);',
            '',
            '		position = intersection;',
            '',
            '		mc_fp_t n_outside = source->n;',
            '',
            '		mc_point3_t voxel_index;',
            '		mcsim_position_to_voxel_index_safe(',
            '			mcsim, &intersection, &voxel_index);',
            '		dbg_print_point3("Launch voxel:", &voxel_index);',
            '		mc_size_t material_index = mcsim_voxel_material_index(',
            '			mcsim, &voxel_index);',
            '		dbg_print_size_t("Launch material:", material_index);',
            '		mc_fp_t n_inside = mc_material_n(',
            '			mcsim_material(mcsim, material_index));',
            '',
            '		refract_safe(&packet_direction, &normal, ',
            '			n_outside, n_inside, &refracted_direction);',
            '',
            '		specular_r = reflectance(n_outside, n_inside, ',
            '				dot3f(&normal, &packet_direction),',
            '				cos_critical(n_outside, n_inside));',
            '',
            '		dbg_print_status(mcsim, "Moved to the sample surface");',
            '	} else {',
            '		position = *mcsim_top_left(mcsim);',
            '		specular_r = FP_1;',
            '		dbg_print_status(mcsim, "Launch missed the sample box!");',
            '	}',
            '',
            '	mcsim_set_position(mcsim, &position);',
            '	mcsim_set_direction(mcsim, &refracted_direction)',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, &position, &packet_direction, specular_r);',
            '	#endif',
            '',
            '	mcsim_set_weight(mcsim, FP_1 - specular_r);',
            '	dbg_print_status(mcsim, "Launch LambertianFiber");',
            '};',
        ))

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`LambertianFiber.cl_type` for a detailed
        list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
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

        T = geometry.transform_base((0.0, 0.0, 1.0), self._direction)

        target.transformation.fromarray(T)

        target.n = self._fiber.ncore
        target.na = self._fiber.na

        target.position.fromarray(self.position)
        target.direction.fromarray(self.direction)

        target.radius = self._fiber.dcore*0.5

        return target, None, None

    def __str__(self):
        return 'LambertianFiber(fiber={}'\
               'position=({}, {}, {}), direction=({}, {}, {}))'.format(
                   self._fiber, *self._position, *self._direction)


class UniformFiberLut(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformFiber(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Parameters
            ----------
            mc: McObject
                A Monte Carlo simulator instance.

            Returns
            -------
            struct: cltypes.Structure
                A structure type that represents a uniform lookup table-based
                optical fiber source in the Monte Carlo kernel.

                The returned structure type implements the following fields:

                - transformation: mc_matrix3f_t
                    Transformation from the beam coordinate system to the Monte
                    Carlo coordinate system.
                - position: mc_point3f_t
                    Source position (axis).
                - direction: mc_point3f_t
                    Source direction (axis).
                - radius: mc_fp_t
                    Radius of the fiber core.
                - n: mc_fp_t
                    Refractive index of the fiber core.
                - lut: mc_fp_lut_t
                    Linear lookup table configuration.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point3f_t),
                ('direction', T.mc_point3f_t),
                ('radius', T.mc_fp_t),
                ('n', T.mc_fp_t),
                ('lut', lututil.LinearLut.cl_type(mc)),
            ]
        return ClUniformFiber

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
            '	mc_fp_t radius;',
            '	mc_fp_t n;',
            '	mc_fp_lut_t lut;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("UniformFiberLut source:");',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
             '	dbg_print_point3f(INDENT "direction:", &src->direction);',
            '	dbg_print_float(INDENT "radius (mm):", src->radius*1e3f);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '	dbg_print_fp_lut(INDENT "lut: ", &src->lut);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_point3f_t pt_src, pt_mc;',
            '	mc_fp_t specular_r = FP_0;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t r = mc_sqrt(mcsim_random(mcsim))*source->radius;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	pt_src.x = r*cos_fi;',
            '	pt_src.y = r*sin_fi;',
            '	pt_src.z = FP_0;',
            '',
            '	transform_point3f(&source->transformation, &pt_src, &pt_mc);',
            '	pt_mc.x += source->position.x;',
            '	pt_mc.y += source->position.y;',
            '	pt_mc.z += source->position.z;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	fp_linear_lut_rel_sample(mcsim_fp_lut_array(mcsim), ',
            '		&source->lut, mcsim_random(mcsim), &cos_theta);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '	/* Adjust the propagation direction to fiber refractive index */',
            '	sin_theta = mc_fdiv(sin_theta, source->n);',
            '',
            '	pt_src.x = cos_fi*sin_theta;',
            '	pt_src.y = sin_fi*sin_theta;',
            '	pt_src.z = cos_theta;',
            '',
            '	mc_point3f_t packet_direction;',
            '	transform_point3f(&source->transformation, &pt_src, &packet_direction);',
            '',
            '	mc_point3f_t source_direction = source->direction;',
            '',
            '	if (mcsim_sample_box_contains(mcsim, &pt_mc)) {',
            '		reverse3f(&source_direction, &source_direction);',
            '		dbg_print("Move to the sample box surface:");',
            '		dbg_print_point3f(INDENT "from        :", &pt_mc);',
            '		dbg_print_point3f(INDENT "in direction:", &source_direction);',
            '	};'
            '',
            '	mc_point3f_t refracted_direction = packet_direction;',
            '	mc_point3f_t position = source->position;',
            '	mc_point3f_t intersection, normal;',
            '	if (mcsim_sample_box_intersect(',
            '			mcsim, &pt_mc, &source_direction, ',
            '			&intersection, &normal)) {',
            '		mcsim_set_position(mcsim, &intersection);',
            '',
            '		reverse3f(&normal, &normal);',
            '',
            '		dbg_print("Surface intersection:");',
            '		dbg_print_point3f(INDENT "position:", &intersection);',
            '		dbg_print_point3f(INDENT "normal:  ", &normal);',
            '',
            '		position = intersection;',
            '',
            '		mc_fp_t n_outside = source->n;',
            '',
            '		mc_point3_t voxel_index;',
            '		mcsim_position_to_voxel_index_safe(',
            '			mcsim, &intersection, &voxel_index);',
            '		dbg_print_point3("Launch voxel:", &voxel_index);',
            '		mc_size_t material_index = mcsim_voxel_material_index(',
            '			mcsim, &voxel_index);',
            '		dbg_print_size_t("Launch material:", material_index);',
            '		mc_fp_t n_inside = mc_material_n(',
            '			mcsim_material(mcsim, material_index));',
            '',
            '		refract_safe(&packet_direction, &normal, ',
            '			n_outside, n_inside, &refracted_direction);',
            '',
            '		specular_r = reflectance(n_outside, n_inside, ',
            '				dot3f(&normal, &packet_direction),',
            '				cos_critical(n_outside, n_inside));',
            '',
            '		dbg_print_status(mcsim, "Moved to the sample surface");',
            '	} else {',
            '		position = *mcsim_top_left(mcsim);',
            '		specular_r = FP_1;',
            '		dbg_print_status(mcsim, "Launch missed the sample box!");',
            '	}',
            '',
            '	mcsim_set_position(mcsim, &position);',
            '	mcsim_set_direction(mcsim, &refracted_direction)',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, &position, &packet_direction, specular_r);',
            '	#endif',
            '',
            '	mcsim_set_weight(mcsim, FP_1 - specular_r);',
            '	dbg_print_status(mcsim, "Launch UniformFiberLut");',
            '};',
        ))

    @staticmethod
    def cl_option(mc: mcobject.McObject) -> mcoptions.RawOptions:
        '''
        This source uses lookup table of floating-point data.
        '''
        return [('MC_USE_FP_LUT', True)]

    def __init__(self, fiber: fiberutil.MultimodeFiberLut,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        An optical fiber photon packet source with emission characteristics
        defined by a lookup table. The lookup table is sampled using
        a uniform random variable and linear interpolation. The obtained
        value represents cosine of the propagation direction with respect to
        the fiber normal. The propagation direction cosines defined in the
        lookup table should be valid for a surrounding medium with refractive
        index 1 (air). The propagation direction is internally
        adjusted according to the refractive index of the medium
        surrounding the optical fiber:
        sin(theta_lut) = n_medium*sin(theta_medium).
        The reflectance at the optical fiber-medium boundary is subtracted
        from the initial weight of the photon packet.

        Parameters
        ----------
        lut: np.ndarray
            Lookup table data (cos(theta_air)).
        diameter: float
            Diameter of the fiber core.
        n: float
            Refractive index of the fiber core.
        position: (float, float, float)
            Center position the fiber source array-like object of size 3
            (x, y, z). The z coordinate is ignored and set to 0.
        direction: (float, float, float)
            Direction of the fiber axis as an array-like object of size 3
            (px, py, pz). If not perpendicular to the sample surface, the
            fiber tip is cut at an angle so that the fiber surface is
            parallel with the sample boundary at the entry point. The z
            component of the direction vector must be positive (fiber pointing
            towards the top surface of the sample) in order to use the top
            surface layout :py:class:`TopGeometryFiber`.

        Note
        ----
        The fiber will be always terminated in a way that forms a tight coupling
        between the sample surface and the fiber tip. If the incidence is
        not normal, the fiber will have an elliptical cross-section (
        cut at an angle).
        The entry point on the sample surface will be determined by propagating
        the position along the given direction (no interactions with the medium
        during this step). If the photon packet does not intersect the sample,
        the initial weight will be set to 0 (reflectance to 1) and the packet
        will be launched from the top-left corner of the voxelized sample box.
        Such zero-weight packets are immediately terminated and have no
        contribution to the fluence and surface detectors, however will be
        included in the trace (have no effect on the sampling volume or other
        trace-based analysis due to zero-weight).
        Note that in case the position lies within the sample, the
        position will be propagated to the entry point using reversed direction.
        From there the packets will be launched according to the NA of the
        fiber and refractive index of the sample surface. The MC simulation
        will start after subtracting the specular reflectance at the
        boundary from the initial weight of the packet.
        '''
        Source.__init__(self)

        self._fiber = fiber
        self._position = np.zeros((3,))
        self._direction = np.zeros((3,))
        self._direction[2] = 1.0

        self._set_position(position)
        self._set_direction(direction)

    def update(self, other: 'UniformFiberLut' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: UniformFiberLut or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, UniformFiberLut):
            self._fiber = other.fiber
            self.position = other.position
            self.direction = other.direction

        elif isinstance(other, dict):
            self._fiber = other.get('fiber', self.fiber)
            self.position = other.get('position', self._position)
            self.direction = other.get('direction', self._direction)

    def _get_fiber(self) -> fiberutil.MultimodeFiberLut:
        return self._fiber
    def _set_fiber(self, fib: fiberutil.MultimodeFiberLut):
        self._fiber = fiberutil.MultimodeFiber(fib)
    fiber = property(_get_fiber, _set_fiber, None,
                     'Multimode lookup table optical fiber.')

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, position: Tuple[float, float, float]):
        self._position[:] = position
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        dir_norm = np.linalg.norm(direction)
        if dir_norm == 0.0:
            raise ValueError('The direction vector is singular!')
        self._direction[:] = direction
        self._direction *= 1.0/dir_norm
    direction = property(_get_direction, _set_direction, None,
                        'Source direction.')

    def cl_pack(self, mc: mcobject.McObject,
                target: cltypes.Structure = None) \
                    -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See the UniformFiberLut.cl_type class for
        a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
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

        T = geometry.transform_base((0.0, 0.0, 1.0), self._direction)

        target.transformation.fromarray(T)

        target.position.fromarray(self._position)
        target.direction.fromarray(self._direction)

        target.radius = self._fiber.dcladding*0.5
        target.n = self._fiber.ncore

        self._fiber.emission.cl_pack(mc, target.lut)

        return target, None, None

    def __str__(self):
        return 'UniformFiberLut(fiber={}, ' \
               'position=({}, {}, {}), direction=({}, {}, {}))'.format(
                   self._lut, *self._position, *self._direction)
