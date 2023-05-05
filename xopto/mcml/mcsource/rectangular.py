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

from xopto.mcml.mcsource.base import Source
from xopto.mcml import mcobject, mctypes, cltypes, mcoptions
from xopto.mcml.mcutil import boundary, geometry
from xopto.mcml.mcutil.lut import EmissionLut

class UniformRectangular(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformRectangular(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position (center of the source),
            size: mc_point2f_t
                Source size along the x and y axis.
            n: mc_fp_t
                Refractive index of the source.
            cos_critical: mc_fp_t
                Cosine of the total internal reflection angle for the
                source -> sample boundary transition.
            cos_min: mc_fp_t
                Cosine of the maximum emission angle in air. 
            layer_index: mc_size_t
                Layer index in which the source is located.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('size', T.mc_point2f_t),
                ('n', T.mc_fp_t),
                ('cos_critical', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('layer_index', T.mc_size_t),
            ]
        return  ClUniformRectangular

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource{',
            '	mc_point3f_t position;',
            '	mc_point2f_t size;',
            '	mc_fp_t n;',
            '	mc_fp_t cos_critical;',
            '	mc_fp_t cos_min;',
            '	mc_size_t layer_index;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("UniformRectangular source:");',
            '	dbg_print_point3(INDENT "position:", &src->position);',
            '	dbg_print_point2f(INDENT "size:", &src->size);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '	dbg_print_float(INDENT "cos_critical:", src->cos_critical);',
            '	dbg_print_float(INDENT "cos_min:", src->cos_min);',
            '	dbg_print_size_t(INDENT "layer_index:", src->layer_index);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_fp_t r;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + (mcsim_random(mcsim) - FP_0p5)*source->size.x,',
            '		source->position.y + (mcsim_random(mcsim) - FP_0p5)*source->size.y,',
            '		source->position.z',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	/* compute propagation direction in air, n_medium = 1 */',
            '	cos_theta = FP_1 - mcsim_random(mcsim)*(FP_1 - source->cos_min);',
            '	sin_theta = mc_sqrt(FP_1 - cos_theta*cos_theta);',
            '	/* adjust the emission angle for the refractive index of the medium */ ',
            '	sin_theta = mc_fdiv(sin_theta, ',
            '		mc_layer_n(mcsim_layer(mcsim, source->layer_index)));',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);'
            '	mcsim_set_direction_coordinates(',
            '		mcsim,',
            '		cos_fi*sin_theta,',
            '		sin_fi*sin_theta,',
            '		cos_theta',
            '	);',
            '',
            '	r = reflectance_cos2(',
            '		source->n, ',
            '		mc_layer_n(mcsim_layer(mcsim, source->layer_index)),',
            '		cos_theta',
            '	);',
            '	mcsim_set_weight(mcsim, FP_1 - r);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = {FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			source->n, &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, source->reflectance);',
            '	#endif',
            '',
            '	mcsim_set_current_layer_index(mcsim, source->layer_index);',
            '};',
        ))

    def __init__(self, width: float, height: float, n: float, na: float,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
        '''
        A rectangular optical photon packet source with emission characteristics
        defined by the numerical aperture 
        The propagation direction is internally adjusted according to the
        refractive index of the medium surrounding the source:
        sin(theta_lut) = n_medium*sin(theta_medium).
        The reflectance at the source-medium boundary is subtracted
        from the initial weight of the photon packet.

        Parameters
        ----------
        width: float
            Source width (x) in (m).

        height: float
            Source height (y) in (m).

        n: float
            Refractive index of the source.

        na: float
            Numerical aperture of the source in air (emission is cut off at
            the NA).

        position: (float, float, float)
            Center of the source beam/rectangle as an array-like object
            of size 3 (px, py, pz).
        '''
        Source.__init__(self)

        self._na = float(na)
        self._width = float(width)
        self._height = float(height)
        self._n = max(1.0, float(n))
        self._position = np.zeros((3,))
        self._set_position(position)

    def update(self, other: 'UniformRectangular', dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: UniformRectangular or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, UniformRectangular):
            self._na = other.na
            self.width = other.width
            self.height = other.height
            self.n = other.n
            self.position = other.position

        elif isinstance(other, dict):
            self.na = other.grt('na', self._na)
            self.width = other.get('width', self._width)
            self.height = other.get('height', self._height)
            self.n = other.get('n', self._n)
            self.position = other.get('position', self._position)

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, position: Tuple[float, float, float]):
        self._position[:] = position
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def _get_width(self) -> float:
        return self._width
    def _set_width(self, width: float):
        self._width = float(width)
    width = property(_get_width, _set_width, None,
                     'Source width (m).')

    def _get_height(self) -> float:
        return self._height
    def _set_height(self, height: float):
        self._height = float(height)
    height = property(_get_height, _set_height, None,
                      'Source height (m).')

    def _get_n(self) -> float:
        return self._n
    def _set_n(self, n: float):
        self._n = max(1.0, float(n))
    n = property(_get_n, _set_n, None,
                 'Source refractive index.')


    def _get_na(self) -> float:
        return self._na
    def _set_na(self, na: float):
        self._na = max(0.0, float(na))
    na = property(_get_na, _set_na, None,
                  'Source numerical aperture (NA).')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator kernel. 
        See the :py:meth:`RectangularLut.cl_type` method for a detailed
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

        if self._position[2] <= 0.0:
            position = (self._position[0], self._position[1], 0.0)
            layer_index = 1
        else:
            position = self._position
            layer_index = mc.layer_index(self._position[2])

        target.position.fromarray(position)

        target.size.fromarray([self._width, self._height])

        target.n = self._n

        target.cos_critical = boundary.cos_critical(
            self._n, mc.layers[layer_index].n)

        target.cos_min = (1 - (self._na)**2)**0.5

        target.layer_index = layer_index

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'width': self._width, 
                'height': self._height,
                'n': self._n,
                'na':  self._na,
                'position': self._position.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'UniformRectangular(width={}, height={}, n={}, na={},'\
               'position=({}, {}, {}))'.format(
                   self._lut, self._width, self._height,
                   self._n, self._na, *self._position)


class LambertianRectangular(UniformRectangular):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClLambertianRectangular(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position (center of the source),
            size: mc_point2f_t
                Source size along the x and y axis.
            n: mc_fp_t
                Refractive index of the source.
            cos_critical: mc_fp_t
                Cosine of the total internal reflection angle for the
                source -> sample boundary transition.
            na: mc_fp_t
                Numerical aperture of the source in air. 
            layer_index: mc_size_t
                Layer index in which the source is located.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('size', T.mc_point2f_t),
                ('n', T.mc_fp_t),
                ('cos_critical', T.mc_fp_t),
                ('na', T.mc_fp_t),
                ('layer_index', T.mc_size_t),
            ]
        return  ClLambertianRectangular

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource{',
            '	mc_point3f_t position;',
            '	mc_point2f_t size;',
            '	mc_fp_t n;',
            '	mc_fp_t cos_critical;',
            '	mc_fp_t na;',
            '	mc_size_t layer_index;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("LambertianRectangular source:");',
            '	dbg_print_point3(INDENT "position:", &src->position);',
            '	dbg_print_point2f(INDENT "size:", &src->size);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '	dbg_print_float(INDENT "cos_critical:", src->cos_critical);',
            '	dbg_print_float(INDENT "na:", src->na);',
            '	dbg_print_size_t(INDENT "layer_index:", src->layer_index);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_fp_t r;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + (mcsim_random(mcsim) - FP_0p5)*source->size.x,',
            '		source->position.y + (mcsim_random(mcsim) - FP_0p5)*source->size.y,',
            '		source->position.z',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	/* compute propagation direction in air, n_medium = 1 */',
            '	sin_theta = mc_sqrt(mcsim_random(mcsim))*source->na;',
            '	/* Adjust the propagation direction to refractive index of the medium */',
            '	sin_theta = mc_fdiv(sin_theta, ',
            '		mc_layer_n(mcsim_layer(mcsim, source->layer_index)));',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);',
            '',
            '	mcsim_set_direction_coordinates(',
            '		mcsim,',
            '		cos_fi*sin_theta,',
            '		sin_fi*sin_theta,',
            '		cos_theta',
            '	);',
            '',
            '	r = reflectance_cos2(',
            '		source->n, ',
            '		mc_layer_n(mcsim_layer(mcsim, source->layer_index)),',
            '		cos_theta',
            '	);',
            '	mcsim_set_weight(mcsim, FP_1 - r);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = {FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			source->n, &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, source->reflectance);',
            '	#endif',
            '',
            '	mcsim_set_current_layer_index(mcsim, source->layer_index);',
            '};',
        ))

    def __init__(self, width: float, height: float, n: float, na: float,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
        '''
        A rectangular optical photon packet source with a Lambertian
        emission characteristics defined by the numerical aperture.
        The propagation direction is internally adjusted according to the
        refractive index of the medium surrounding the source:
        sin(theta_lut) = n_medium*sin(theta_medium).
        The reflectance at the source-medium boundary is subtracted
        from the initial weight of the photon packet.

        Parameters
        ----------
        width: float
            Source width (x) in (m).

        height: float
            Source height (y) in (m).

        n: float
            Refractive index of the source.

        na: float
            Numerical aperture of the source in air (emission is cut off at
            the NA).

        position: (float, float, float)
            Center of the source beam/rectangle as an array-like object
            of size 3 (px, py, pz).
        '''
        Source.__init__(self)

        self._na = float(na)
        self._width = float(width)
        self._height = float(height)
        self._n = max(1.0, float(n))
        self._position = np.zeros((3,))
        self._set_position(position)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator kernel. 
        See the :py:meth:`RectangularLut.cl_type` method for a detailed
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

        if self._position[2] <= 0.0:
            position = (self._position[0], self._position[1], 0.0)
            layer_index = 1
        else:
            position = self._position
            layer_index = mc.layer_index(self._position[2])

        target.position.fromarray(position)

        target.size.fromarray([self._width, self._height])

        target.n = self._n

        target.cos_critical = boundary.cos_critical(
            self._n, mc.layers[layer_index].n)

        target.na = self._na

        target.layer_index = layer_index

        return target, None, None

    def __str__(self):
        return 'LambertianRectangular(width={}, height={}, n={}, na={},'\
               'position=({}, {}, {}))'.format(
                   self._width, self._height,
                   self._n, self._na, *self._position)


class UniformRectangularLut(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClUniformRectangularLut(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Fields
            ------
            position: mc_point3f_t
                Source position (center of the source),
            size: mc_point2f_t
                Source size along the x and y axis.
            n: mc_fp_t
                Refractive index of the source.
            cos_critical: mc_fp_t
                Cosine of the total internal reflection angle for the
                source -> sample boundary transition.
            lut: cltypes.Structure
                Linear lookup table structure. 
            layer_index: mc_size_t
                Layer index in which the source is located.
            '''
            _fields_ = [
                ('position', T.mc_point3f_t),
                ('size', T.mc_point2f_t),
                ('n', T.mc_fp_t),
                ('cos_critical', T.mc_fp_t),
                ('lut', EmissionLut.cl_type(mc)),
                ('layer_index', T.mc_size_t),
            ]
        return  ClUniformRectangularLut

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McSource{',
            '	mc_point3f_t position;',
            '	mc_point2f_t size;',
            '	mc_fp_t n;',
            '	mc_fp_t cos_critical;',
            '	mc_fp_lut_t lut;',
            '	mc_size_t layer_index;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("UniformRectangularLut source:");',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
            '	dbg_print_point2f(INDENT "size:", &src->size);',
            '	dbg_print_float(INDENT "n:", src->n);',
            '   dbg_print_fp_lut(INDENT "lut: ", &src->lut);',
            '	dbg_print_float(INDENT "cos_critical:", src->cos_critical);',
            '	dbg_print_size_t(INDENT "layer_index:", src->layer_index);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '	mc_fp_t r;',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + (mcsim_random(mcsim) - FP_0p5)*source->size.x,',
            '		source->position.y + (mcsim_random(mcsim) - FP_0p5)*source->size.y,',
            '		source->position.z',
            '	);',
            '',
            '	mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '	/* compute propagation direction in air, n_medium = 1 */',
            '	cos_theta = FP_0;',
            '	fp_linear_lut_rel_sample(',
            '		mcsim_fp_lut_array(mcsim), &source->lut, ',
            '		mcsim_random(mcsim), &cos_theta);',
            '',
            '	/* Compute sine of the propagation direction and adjust ',
            '	 * for the refractive index of the surrounding medium.  */ ',
            '	sin_theta = mc_fdiv(mc_sqrt(FP_1 - cos_theta*cos_theta), ',
            '		mc_layer_n(mcsim_layer(mcsim, source->layer_index)));',
            '	cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);'
            '	mcsim_set_direction_coordinates(',
            '		mcsim,',
            '		cos_fi*sin_theta,',
            '		sin_fi*sin_theta,',
            '		cos_theta',
            '	);',
            '',
            '	r = reflectance_cos2(',
            '		source->n, ',
            '		mcsim_layer(mcsim, source->layer_index)->n,',
            '		cos_theta',
            '	);',
            '	mcsim_set_weight(mcsim, FP_1 - r);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = {FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			source->n, &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, source->reflectance);',
            '	#endif',
            '',
            '	mcsim_set_current_layer_index(mcsim, source->layer_index);',
            '};',
        ))

    @staticmethod
    def cl_options(mc: mcobject.McObject) -> mcoptions.RawOptions:
        return [('MC_USE_FP_LUT', True)]

    def __init__(self, lut: EmissionLut, width: float, height: float, n: float,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
        '''
        A rectangular optical photon packet source with emission characteristics
        defined by a lookup table. The lookup table is sampled using
        a uniform random variable and linear interpolation. The obtained
        value represents cosine of the propagation direction with respect to
        the source normal. The propagation direction cosines defined in the
        lookup table should be valid for a surrounding medium with a
        refractive index 1 (air).
        The propagation direction is internally adjusted according to the
        refractive index of the medium surrounding the source:
        sin(theta_lut) = n_medium*sin(theta_medium).
        The reflectance at the source-medium boundary is subtracted
        from the initial weight of the photon packet.

        Parameters
        ----------
        lut: EmissionLut
            Lookup table of the angular emission.

        width: float
            Source width (x) in (m).

        height: float
            Source height (y) in (m).

        n: float
            Refractive index of the source.

        position: (float, float, float)
            Center of the source beam/rectangle as an array-like object
            of size 3 (px, py, pz).
        '''
        Source.__init__(self)

        self._lut = lut
        self._width = float(width)
        self._height = float(height)
        self._n = max(1.0, float(n))
        self._position = np.zeros((3,))
        self._set_position(position)

    def update(self, other: 'UniformRectangularLut' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: UniformRectangularLut or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, UniformRectangularLut):
            self._lut = EmissionLut(other.lut)
            self._offset = 0
            self.width = other.width
            self.height = other.height
            self.n = other.n
            self.position = other.position

        elif isinstance(other, dict):
            new_lut = other.get('lut')
            if new_lut is not None:
                self._lut = new_lut
                self._offset = 0
            self.width = other.get('width', self._width)
            self.height = other.get('height', self._height)
            self.n = other.get('n', self._n)
            self.position = other.get('position', self._position)

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, position: Tuple[float, float, float]):
        self._position[:] = position
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def _get_width(self) -> float:
        return self._width
    def _set_width(self, width: float):
        self._width = float(width)
    width = property(_get_width, _set_width, None,
                     'Source width (m).')

    def _get_height(self) -> float:
        return self._height
    def _set_height(self, height: float):
        self._height = float(height)
    height = property(_get_height, _set_height, None,
                      'Source height (m).')

    def _get_lut(self) -> EmissionLut:
        return self._lut
    lut = property(_get_lut, None, None,
                   'Lookup table of angular sensitivity.')

    def _get_n(self) -> float:
        return self._n
    def _set_n(self, n: float):
        self._n = max(1.0, float(n))
    n = property(_get_n, _set_n, None,
                 'Source refractive index.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> Tuple[cltypes.Structure, None, None]:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator kernel. 
        See the :py:meth:`UniformRectangularLut.cl_type` method for a detailed
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

        if self._position[2] <= 0.0:
            position = (self._position[0], self._position[1], 0.0)
            layer_index = 1
        else:
            position = self._position
            layer_index = mc.layer_index(self._position[2])

        self._lut.cl_pack(mc, target.lut)

        target.position.fromarray(position)

        target.size.fromarray([self._width, self._height])

        target.n = self._n

        target.cos_critical = boundary.cos_critical(
            self._n, mc.layers[layer_index].n)

        target.layer_index = layer_index

        return target, None, None

    def todict(self) -> dict:
        '''
        Save the source configuration data to a dictionary.
        Use the :meth:`UniformRectangularLut.fromdict` method to create a new
        source instance from the returned data.

        Returns
        -------
        data: dict
            Source configuration as a dictionary.
        '''
        return {
            'type':'UniformRectangularLut',
            'lut':self._lut.todict(),
            'width':self._width,
            'height':self._height,
            'n':self._n,
            'position':self._position.tolist()
        }

    @staticmethod
    def fromdict(data) -> 'UniformRectangularLut':
        '''
        Create a source instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`Radial.todict` method.
        '''
        source_type = data.pop('type')
        if source_type != 'UniformRectangularLut':
            raise TypeError(
                'Expected "UniformRectangularLut" type bot got "{}"!'.format(
                    source_type))
        lut = EmissionLut.fromdict(data.pop('lut'))

        return UniformRectangularLut(lut, **data)

    def __str__(self):
        return 'UniformRectangularLut(lut={}, width={}, height={}, n={}, '\
               'position=({}, {}, {}))'.format(
                   self._lut, self._width, self._height,
                   self._n, *self._position)

    def __repr__(self):
        return '{} #{}'.format(self, id(self))
