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
from xopto.mcml import mcobject
from xopto.mcml import cltypes
from xopto.mcml import mctypes
from xopto.mcml.mcutil import boundary, geometry

class GaussianBeam(Source):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Structure that is passed to the Monte carlo simulator kernel.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents the Gaussian beam in
            the Monte Carlo kernel.

            The returned structure type implements the following fields:

            - transformation: mc_matrix3f_t
                Transformation from the beam coordinate system to the Monte
                Carlo coordinate system.
            - position: (float, float, float)
                Center of the collimated beam as an array-like object of size 3
                (x, y, z). The beam will be always propagated to the top
                sample surface.
            - direction: (float, float, float)
                Direction of the collimated beam as an array-like object of size 3
                (px, py, pz). The vector should be normalized to unit length and
                have a positive z coordinate (hitting the top sample surface).
            - sigma: mc_point2f_t
                Width of the Gaussian beam in terms of standard deviation along
                the x and y axis,
            - clip: mc_fp_t
                Clips the gaussian beam at the specified number of sigmas.
            - reflectance: mc_fp_t
                Precalculated reflectance at the source-sample boundary (c_float).
            - layer_index: mc_int_t
                Index of the sample layer into which this source is launching
                photon packets.
        '''
        T = mc.types
        class ClGaussianBeam(cltypes.Structure):
            _pack_ = 1
            _fields_ = [
                ('transformation', T.mc_matrix3f_t),
                ('position', T.mc_point3f_t),
                ('direction', T.mc_point3f_t),
                ('sigma', T.mc_point2f_t),
                ('clip', T.mc_fp_t),
                ('reflectance', T.mc_fp_t),
            ]
        return ClGaussianBeam

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
            '	mc_point2f_t sigma;',
            '	mc_fp_t clip;',
            '	mc_fp_t reflectance;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the source in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_source(__mc_source_mem const McSource *src){',
            '	dbg_print("GaussianBeam source:");',
            '	dbg_print_matrix3f(INDENT "transformation:", &src->transformation);',
            '	dbg_print_point3f(INDENT "position:", &src->position);',
            '	dbg_print_point3f(INDENT "direction:", &src->direction);',
            '	dbg_print_point2f(INDENT "sigma:", &src->sigma);',
            '	dbg_print_float(INDENT "clip:", src->clip);',
            '	dbg_print_float(INDENT "reflectance:", src->reflectance);',
            '};',
            '',
            'inline void mcsim_launch(McSim *mcsim){',
            '	__mc_source_mem const struct McSource *source = mcsim_source(mcsim);',
            '	mc_fp_t cos_fi, sin_fi, r;',
            '	mc_point3f_t pt_src, pt_mc;',
            '',
            '	/* r = sigma*np.sqrt(-FP_2*np.log(FP_1 - uniform_random)) */',
            '	r = mc_sqrt(-FP_2*mc_log(FP_1 - mcsim_random(mcsim)));',
            '	r = mc_fmin(r, source->clip);',
            '',
            '	mc_sincos(FP_2PI*mcsim_random(mcsim), &sin_fi, &cos_fi);',
            '	pt_src.x = r*cos_fi*source->sigma.x;',
            '	pt_src.y = r*sin_fi*source->sigma.y;',
            '	pt_src.z = FP_0;',
            '',
            '	mc_matrix3f_t transformation = source->transformation;',
            '	transform_point3f(&transformation, &pt_src, &pt_mc);',
            '	mc_fp_t k = mc_fdiv(FP_0 - pt_mc.z, source->direction.z);',
            '	pt_mc.x += k*source->direction.x;',
            '	pt_mc.y += k*source->direction.y;',
            '	pt_mc.z = FP_0;',
            '',
            '	mcsim_set_position_coordinates(',
            '		mcsim,',
            '		source->position.x + pt_mc.x,',
            '		source->position.y + pt_mc.y,',
            '		FP_0',
            '	);',
            '	mcsim_set_direction(mcsim, &source->direction);',
            '	mcsim_set_current_layer_index(mcsim, 1);',
            '	mcsim_set_weight(mcsim, FP_1 - source->reflectance);',
            '',
            '	#if MC_USE_SPECULAR_DETECTOR',
            '		mc_point3f_t dir_in = {mcsim_direction_x(mcsim), ',
            '			mcsim_direction_y(mcsim), -mcsim_direction_z(mcsim)};',
            '		mc_point3f_t dir;',
            '		mc_point3f_t normal = (mc_point3f_t){FP_0, FP_0, -FP_1};',
            '		refract(&dir_in, &normal, mc_layer_n(mcsim_layer(mcsim, 1)),',
            '			mc_layer_n(mcsim_layer(mcsim, 0)), &dir);',
            '		mcsim_specular_detector_deposit(',
            '			mcsim, mcsim_position(mcsim), &dir, source->reflectance);',
            '	#endif',
            '',
            '	dbg_print_status(mcsim, "Launch GaussianBeam");',
            '};',
        ))

    @staticmethod
    def fwhm2sigma(fwhm: float) -> float:
        '''
        Converts sigma to Full width at half maximum (FWHM).

        Parameters
        ----------
        sigma: float
            Standard deviation of the Gaussian.

        Returns
        -------
        fwhm: float
            FWHM parameter of the Gaussian
        '''
        return fwhm/(8*np.log(2))**0.5

    @staticmethod
    def sigma2fwhm(sigma: float) -> float:
        '''
        Converts Full width at half maximum (FWHM) to sigma (standard deviation).

        Parameters
        ----------
        fwhm: float
            Full width at half maximum (FWHM).

        Returns
        -------
        sigma: float
            Standard deviation of the Gaussian.
        '''
        return sigma*(8*np.log(2))**0.5
    
    def __init__(self, sigma: float or Tuple[float, float], clip: float = 5.0,
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Collimated Gaussian beam photon packet source.

        Parameters
        ----------
        sigma: float or (float, float)
            Collimated beam width in terms of standard deviation given along
            the x and y axis. If a single value is provided, the same width
            is used along the x and y axis.
        clip: float
            Clip the beam at clip*sigma distance from the beam axis.
        position: (float, float, float)
            Center of the collimated beam as an array-like object of size 3
            (x, y, z).
        direction: (float, float, float)
            Direction of the collimated beam as an array-like object of size 3
            (px, py, pz). The vector should be normalized to unit length and
            have a positive z coordinate (hitting the top sample surface).

        Note
        ----
        The beam will be first propagated from the given position to the
        entry point on the sample surface along the propagation
        direction (no interactions with the medium during this step).
        Note that in case the position lies within the sample, the
        beam will be propagated to the entry point using reversed direction.
        From there it will be refracted into the sample. The MC simulation
        will start after subtracting the specular reflectance at the
        sample boundary from the initial weight of the packet.
        '''
        Source.__init__(self)
        self._clip = self._sigma = None
        self._position = np.zeros((3,))
        self._direction = np.zeros((3,))
        self._direction[2] = 1.0
        self._sigma = np.zeros((2,))

        self._set_sigma(sigma)
        self._set_clip(clip)
        self._set_position(position)
        self._set_direction(direction)

    def update(self, other: 'GaussianBeam' or dict):
        '''
        Update this source configuration from the other source. The
        other source must be of the same type as this source or a dict with
        appropriate fields.

        Parameters
        ----------
        other: GaussianBeam or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, GaussianBeam):
            self.sigma = other.sigma
            self.clip = other.clip
            self.position = other.position
            self.direction = other.direction
        elif isinstance(other, dict):
            self.sigma = other.get('sigma', self.sigma)
            self.clip = other.get('clip', self.clip)
            self.position = other.get('position', self.position)
            self.direction = other.get('direction', self.direction)

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, position: Tuple[float, float, float]):
        self._position[:] = position
    position = property(_get_position, _set_position, None,
                        'Source position.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('The norm/length of the propagation direction '
                             'vector must not be 0!')
        self._direction *= 1.0/norm
        if self._direction[-1] <= 0.0:
            raise ValueError('Z component of the propagation direction '
                             'must be positive!')
    direction = property(_get_direction, _set_direction, None,
                        'Source direction.')

    def _get_sigma(self) -> Tuple[float, float]:
        return self._sigma
    def _set_sigma(self, sigma: float or Tuple[float, float]):
        self._sigma[:] = sigma
        if np.any(self._sigma < 0.0):
            raise ValueError('Beam diameter/sigma must not be negative!')
    sigma = property(_get_sigma, _set_sigma, None,
                     'Beam standard deviation (m).')

    def _get_clip(self) -> float:
        return self._clip
    def _set_clip(self, clip: Tuple[float, float]):
        self._clip = float(clip)
        if self._clip < 0.0:
            raise ValueError('Clip diameter/sigma must be greater than zero!.')
    clip = property(_get_clip, _set_clip, None,
                    'Number of standard deviations at which '
                    'the beam is clipped.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`GaussianBeam.cl_type`
        for a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: pyopyo.mcml.mcsource.GaussianBeam.cl_type
            Structure that is filled with the source data.

        Returns
        -------
        target: pyopyo.mcml.mcsource.GaussianBeam.cl_type
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        topgeometry: None
            This source does not use advanced geometry at the top sample surface.
        bottomgeometry: None
            This source does not use advanced geometry at the bottom sample surface.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        # propagate the beam to the top sample surface
        k = (0.0 - self._position[2])/self._direction[2]
        position = self._position + k*self._direction
        position[2] = 0.0

        direction = boundary.refract(self._direction, (0.0, 0.0, 1.0),
                                     mc.layers[0].n, mc.layers[1].n)

        reflectance = boundary.reflectance(
            mc.layers[0].n, mc.layers[1].n, abs(self._direction[-1]))

        T = geometry.transform_base((0.0, 0.0, 1.0), self._direction)

        target.transformation.fromarray(T)

        target.position.fromarray(position)
        target.direction.fromarray(direction)

        target.sigma.fromarray(self._sigma)

        target.clip = self._clip

        target.reflectance = reflectance

        return target, None, None

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'sigma': self._sigma.tolist(), 
                'clip': self._clip,
                'position': self._position.tolist(),
                'direction': self._direction.tolist(),
                'type': self.__class__.__name__}

    def __str__(self):
        return 'GaussianBeam(sigma=({}, {}), clip={}, '\
               'position=({}, {}, {}), direction=({}, {}, {}))'.format(
                   *self._sigma, self._clip, *self._position, *self._direction)
