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

from xopto.mcml.mcsource.base import Source
from xopto.mcml import mcobject
from xopto.mcml import mctypes
from xopto.mcml import mcoptions
from xopto.mcml.mcutil import boundary
from xopto.mcml.mcutil import geometry
from xopto.mcml.mcutil import lut as lututil
from xopto.mcml.mcutil import fiber as fiberutil
from xopto.mcml import cltypes

class TopGeometryFiberNI(mcobject.McObject):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Top sample surface geometry representing an optical fiber with
        normal incidence.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents the top geometry in
            the Monte Carlo kernel.

        Fields
        ------
        position: mc_point3f_t
            Position of the optical fiber in the Monte Carlo coordinate system.
        core_radius: mc_fp_t
            Outer radius of the fiber core.
        cladding_radius: mc_fp_t
            Outer radius of the fiber cladding.
        core_n: mc_fp_t
            Refractive index of the fiber core.
        cladding_n: mc_fp_t
            Refractive index of the fiber cladding.
        core_cos_critical: mc_fp_t
            Cosine of the total internal reflection angle for the transition
            sample -> fiber core.
        cladding_cos_critical: mc_fp_t
            Cosine of the total internal reflection angle for the transition
            sample -> fiber cladding.
        '''
        T = mc.types
        class ClTopGeometryFiberNI(cltypes.Structure):
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('core_radius', T.mc_fp_t),
                ('cladding_radius', T.mc_fp_t),
                ('core_n', T.mc_fp_t),
                ('cladding_n', T.mc_fp_t),
                ('core_cos_critical', T.mc_fp_t),
                ('cladding_cos_critical', T.mc_fp_t),
            ]
        return ClTopGeometryFiberNI

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McComplexSurfaceTop{',
            '	mc_point3f_t position;',
            '	mc_fp_t core_radius;',
            '	mc_fp_t cladding_radius;',
            '	mc_fp_t core_n;',
            '	mc_fp_t cladding_n;',
            '	mc_fp_t core_cos_critical;',
            '	mc_fp_t cladding_cos_critical;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        return '\n'.join((
            'void dbg_print_top_geometry(__mc_geometry_mem const McComplexSurfaceTop *geo){',
            '	dbg_print("FiberTopGeometryNI:");',
            '	dbg_print_point3f(INDENT "position:", &geo->position);',
            '	dbg_print_float(INDENT "core_radius (mm):", src->core_radius*1e3f);',
            '	dbg_print_float(INDENT "cladding_radius (mm):", src->cladding_radius*1e3f);',
            '	dbg_print_float(INDENT "core_n:", src->core_n);',
            '	dbg_print_float(INDENT "cladding_n:", src->cladding_n);',
            '	dbg_print_float(INDENT "core_cos_critical:", src->core_cos_critical);',
            '	dbg_print_float(INDENT "cladding_cos_critical:", src->cladding_cos_critical);',
            '};',
            '',
            'inline mc_int_t mcsim_top_geometry_handler(McSim *psim, mc_fp_t *n2, mc_fp_t *cc){',
            '	__mc_geometry_mem const struct McComplexSurfaceTop *geometry = mcsim_top_geometry(mcsim);',
            '',
            '	mc_point3f_t pt_mc;',
            '	pt_mc.x = mcsim_position_x(psim) - geometry->position.x;',
            '	pt_mc.y = mcsim_position_y(psim) - geometry->position.y;',
            '	pt_mc.z = FP_0;',
            '',
            '	mc_fp_t r2 = pt_mc.x*pt_mc.x + pt_mc.y*pt_mc.y;',
            '	if (dr2 <= geometry->cladding_radius*geometry->cladding_radius){',
            '		*n2 = geometry->cladding_n;',
            '		*cc = geometry->cladding_cos_critical;',
            '		if (dr2 < geometry->core_radius*geometry->core_radius){',
            '			*n2 = geometry->core_n;',
            '			*cc = geometry->core_cos_critical;',
            '		};'
            '	}',
            '	return MC_SURFACE_GEOMETRY_CONTINUE;',
            '};',
        ))

    def __init__(self, fiber: 'UniformFiberNI' or 'LambertianFiberNI' or
                              'UniformFiberLutNI'):
        self._fiber_src = fiber

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`FiberGeometryTopNI.cl_type` for a detailed
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

        target.position.fromarray(fiber_src.position)
        target.core_radius = fiber.dcore*0.5
        target.cladding_radius = fiber.dcladding*0.5
        target.core_n = fiber.ncore
        target.cladding_n = fiber.ncladding
        target.core_cos_critical = boundary.cos_critical(
            mc.layers[1].n, fiber.ncore)
        target.cladding_cos_critical = boundary.cos_critical(
            mc.layers[1].n, fiber.ncladding)

        return target

class UniformFiberNI(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformFiberNI(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Parameters
            ----------
            mc: McObject
                A Monte Carlo simulator instance.

            Returns
            -------
            struct: cltypes.Structure
                A structure type that represents the scattering phase function in
                the Monte Carlo kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position (axis).
            radius: mc_fp_t
                Radius of the fiber core.
            cos_min: mc_fp_t
                Cosine of the fiber acceptance angle in air.
            n: mc_fp_t
                Refractive index of the fiber core.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('radius', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n', T.mc_fp_t),
            ]
        return ClUniformFiberNI

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource{',
            '	mc_point3f_t position;',
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
            '	dbg_print("UniformFiberNI source:");',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
            '	dbg_print_float(INDENT "radius (mm):", src->radius*1e3f);',
            '	dbg_print_float(INDENT "cos_min:", src->cos_min);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '',
            '	mc_fp_t r = mc_sqrt(mcsim_random(mcsim))*source->radius;',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + r*cos_fi,',
            '		source->position.y + r*sin_fi,',
            '		FP_0',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	cos_theta = FP_1 - mcsim_random(mcsim)*(FP_1 - source->cos_min);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '	/* adjust the emission angle for the refractive index of the sample */ ',
            '	sin_theta = mc_fdiv(sin_theta, mc_layer_n(mcsim_layer(mcsim, 1)));',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);',
            '',
            '	mc_point3f_t pt_src = {',
            '		cos_fi*sin_theta,',
            '		sin_fi*sin_theta,',
            '		cos_theta',
            '	};'
            '	mcsim_set_direction(mcsim, &pt_src);',
            '',
            '	r = reflectance_cos2(',
            '		source->n, ',
            '		mc_layer_n(mcsim_layer(mcsim, 1)),',
            '		pt_src.z',
            '	);',
            '	mcsim_set_weight(mcsim, FP_1 - r);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = (mc_point3f_t){FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			source->n, &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, r);',
            '	#endif',
            '',
            '	mcsim_set_current_layer_index(mcsim, 1);',
            '};',
        ))

    def __init__(self, fiber: fiberutil.MultimodeFiber,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
        '''
        An optical fiber photon packet source.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents the scattering phase function in
            the Monte Carlo kernel.

        ----------
        fiber: MultimodeFiber
            Optical fiber parameters.
        position: (float, float, float)
            Center position the fiber source array-like object of size 3
            (x, y, z). The z coordinate is ignored and set to 0.
        '''
        Source.__init__(self)

        self._fiber = fiber
        self._position = np.zeros((3,))

        self._set_position(position)

    def update(self, other: 'UniformFiberNI' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: xopto.mcml.mcsource.UniformFiber or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, UniformFiberNI):
            self.fiber = other.fiber
            self.position = other.position
        elif isinstance(other, dict):
            self.fiber = other.get('fiber', self.fiber)
            self.position = other.get('position', self.position)

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
        self._position[2] = 0.0
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`UniformFiberNI.cl_type` for a detailed
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

        target.n = self._fiber.ncore
        target.cos_min = (1.0 - (self._fiber.na)**2)**0.5

        target.position.fromarray(self.position)

        target.radius = self._fiber.dcore*0.5

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'fiber': self._fiber.todict(), 
                'position': self._position.tolist(),
                'type': self.__class__.__name__}

    @classmethod
    def fromdict(cls, data: dict) -> 'UniformFiberNI':
        '''
        Create a new instance of a photon packet source from a dict that was
        created by the :py:meth:`todict` method. 
        '''
        data_ = dict(data)
        fiber_data = data_.pop('fiber')
        data_['fiber'] = getattr(
            fiberutil, fiber_data['type']).fromdict(fiber_data)
        return super().fromdict(data_)

    def __str__(self):
        return 'UniformFiberNI(fiber={}, position=({}, {}, {}))'.format(
            self._fiber *self._position)


class LambertianFiberNI(UniformFiberNI):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformFiberNI(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Parameters
            ----------
            mc: McObject
                A Monte Carlo simulator instance.

            Returns
            -------
            struct: cltypes.Structure
                A structure type that represents Lambertian fiber with normal
                incidence in the Monte Carlo kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position (axis).
            radius: mc_fp_t
                Radius of the fiber core.
            na: mc_fp_t
                Numerical aperture of the fiber core - sine of the acceptance
                angle in air.
            n: mc_fp_t
                Refractive index of the fiber core.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('radius', T.mc_fp_t),
                ('na', T.mc_fp_t),
                ('n', T.mc_fp_t),
            ]
        return ClUniformFiberNI

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource{',
            '	mc_point3f_t position;',
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
            '	dbg_print("LambertianFiberNI source:");',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
            '	dbg_print_float(INDENT "radius (mm):", src->radius*1e3f);',
            '	dbg_print_float(INDENT "na:", src->na);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_point3f_t pt_mc;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t r = mc_sqrt(mcsim_random(mcsim))*source->radius;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + r*cos_fi,',
            '		source->position.y + r*sin_fi,',
            '		FP_0',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	sin_theta = mc_sqrt(mcsim_random(mcsim))*source->na;',
            '	/* Adjust the propagation direction to the sample refractive index */',
            '	sin_theta = mc_fdiv(sin_theta, mc_layer_n(mcsim_layer(mcsim, 1)));',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);',
            '',
            '	pt_mc.x = cos_fi*sin_theta;',
            '	pt_mc.y = sin_fi*sin_theta;',
            '	pt_mc.z = cos_theta;',
            '	mcsim_set_direction(mcsim, &pt_mc);',
            '',
            '	r = reflectance_cos2(',
            '		source->n, ',
            '		mc_layer_n(mcsim_layer(mcsim, 1)),',
            '		pt_mc.z',
            '	);',
            '	mcsim_set_weight(mcsim, FP_1 - r);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = (mc_point3f_t){FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			source->n, &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, r);',
            '	#endif',
            '',
            '	mcsim_set_current_layer_index(mcsim, 1);',
            '};',
        ))

    def cl_pack(self, mc: mcobject.McObject,
                target: cltypes.Structure = None) \
                    -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`LambertianFiberNI.cl_type` for a detailed
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

        target.n = self._fiber.ncore
        target.na = self._fiber.na

        target.position.fromarray(self.position)

        target.radius = self._fiber.dcore*0.5

        return target, None, None

    def __str__(self):
        return 'LambertianFiberNI(fiber={}, position=({}, {}, {}))'.format(
            self._fiber, *self._position)


class UniformFiberLutNI(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformFiberLutNI(cltypes.Structure):
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
                fiber with normal incidence in the Monte Carlo kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position (axis).
            radius: mc_fp_t
                Radius of the fiber core.
            n: mc_fp_t
                Refractive index of the fiber core.
            lut: mc_fp_lut_t
                Linear lookup table configuration.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('radius', T.mc_fp_t),
                ('n', T.mc_fp_t),
                ('lut', lututil.LinearLut.cl_type(mc)),
            ]
        return ClUniformFiberLutNI

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource{',
            '	mc_point3f_t position;',
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
            '	dbg_print("UniformFiberLutNI source:");',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
            '	dbg_print_float(INDENT "radius (mm):", src->radius*1e3f);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '	dbg_print_fp_lut(INDENT "lut: ", &src->lut);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_point3f_t pt_mc;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t r = mc_sqrt(mcsim_random(mcsim))*source->radius;',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + r*cos_fi,',
            '		source->position.y + r*sin_fi,',
            '		FP_0',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	fp_linear_lut_rel_sample(mcsim_fp_lut_array(mcsim), ',
            '		&source->lut, mcsim_random(mcsim), &cos_theta);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '	/* Adjust the propagation direction to the sample refractive index */',
            '	sin_theta = mc_fdiv(sin_theta, mc_layer_n(mcsim_layer(mcsim, 1)));',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);',
            '',
            '	pt_mc.x = cos_fi*sin_theta;',
            '	pt_mc.y = sin_fi*sin_theta;',
            '	pt_mc.z = cos_theta;',
            '',
            '	r = reflectance_cos2(',
            '		source->n, ',
            '		mc_layer_n(mcsim_layer(mcsim, 1)),',
            '		pt_mc.z',
            '	);',
            '	mcsim_set_weight(mcsim, FP_1 - r);',
            '	mcsim_set_direction(mcsim, &pt_mc);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = (mc_point3f_t){FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			source->n, &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, r);',
            '	#endif',
            '',
            '	mcsim_set_current_layer_index(mcsim, 1);',
            '',
            '	dbg_print_status(mcsim, "Packet launched");',
            '};',
        ))

    @staticmethod
    def cl_options(mc: mcobject.McObject) -> mcoptions.RawOptions:
        '''
        This source uses lookup table of floating-point data.
        '''
        return [('MC_USE_FP_LUT', True)]

    def __init__(self, fiber: fiberutil.MultimodeFiberLut,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
        '''
        An optical fiber photon packet source with emission characteristics
        defined by a lookup table. The lookup table is sampled using
        a uniform random variable and linear interpolation. The obtained
        value represents cosine of the propagation direction with respect to
        the fiber normal. The propagation direction cosines defined in the
        lookup table should be valid for a surrounding medium with a refractive
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
            Center of the fiber source give as a tuple (x, y, z).
            The z coordinate is ignored and set to 0.
        '''
        Source.__init__(self)

        self._fiber = fiber
        self._position = np.zeros((3,))

        self._set_position(position)

    def update(self, other: 'UniformFiberLutNI' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: UniformFiberLutNI or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, UniformFiberLutNI):
            self._fiber = other.fiber
            self.position = other.position

        elif isinstance(other, dict):
            self._fiber = other.get('fiber', self.fiber)
            self.position = other.get('position', self._position)

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
        self._position[2] = 0.0
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See the UniformFiberLutNI.cl_type class for
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

        target.position.fromarray(self.position)

        target.radius = self._fiber.dcore*0.5
        target.n = self._fiber.ncore
        self._fiber.emission.cl_pack(mc, target.lut)

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'fiber': self._fiber.todict(), 
                'position': self._position.tolist(),
                'type': self.__class__.__name__}

    @classmethod
    def fromdict(cls, data: dict) -> 'UniformFiberNI':
        '''
        Create a new instance of a photon packet source from a dict that was
        created by the :py:meth:`todict` method. 
        '''
        data_ = dict(data)
        fiber_data = data_.pop('fiber')
        data_['fiber'] = getattr(
            fiberutil, fiber_data['type']).fromdict(fiber_data)
        return super().fromdict(data_)

    def __str__(self):
        return 'UniformFiberLut(fiber={}, position=({}, {}, {}))'.format(
            self._lut, *self._position)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id())
