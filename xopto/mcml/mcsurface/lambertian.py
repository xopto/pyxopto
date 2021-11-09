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

import numpy as np

from xopto.mcml import mcobject
from xopto.mcml import cltypes
from .base import SurfaceLayoutAny


class LambertianReflector(SurfaceLayoutAny):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClLambertianReflector(cltypes.Structure):
            '''
            Structure that is passed to the Monte carlo simulator kernel.

            Fields
            ------
            reflectance: mc_fp_t
                Total reflectance of the Lambertian reflector. Use a value between
                0 and 1.
            specular: mc_fp_t
                Fraction of specular reflections. Use a value of 0 for a
                pure lambertian reflector.
            '''
            _fields_ = [
                ('reflectance', T.mc_fp_t),
                ('specular', T.mc_fp_t)
            ]
        return ClLambertianReflector

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the reflector in the Monte Carlo simulator.
        '''
        loc = self.location.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}SurfaceLayout{{'.format(loc),
            '	mc_fp_t reflectance;',
            '	mc_fp_t specular;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the reflector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_surface_layout('.format(loc),
            '		__mc_surface_mem const Mc{}SurfaceLayout *layout){{'.format(Loc),
            '	dbg_print("Mc{}SurfaceLayout - Lambertian reflector:");'.format(Loc),
            '	dbg_print_float(INDENT "reflectance:", layout->reflectance);',
            '	dbg_print_float(INDENT "specular:", layout->specular);',
            '};',
            '',
            'inline int mcsim_{}_surface_layout_handler('.format(loc),
            '		McSim *mcsim, mc_fp_t *n2, mc_fp_t *cc){',
            '	__mc_surface_mem const struct Mc{}SurfaceLayout *layout = '.format(Loc),
            '		mcsim_{}_surface_layout(mcsim);'.format(loc),
            '	mc_fp_t sin_fi, cos_fi, sin_theta, cos_theta;',
            '',
            '	dbg_print_point3f("{} LambertianReflector hit: ", '.format(Loc),
            '		mcsim_direction(mcsim));',
            '',
            '	if (mcsim_random(mcsim) > layout->specular){',
            '		dbg_print("{} LambertianReflector: Lambertian reflection");'.format(Loc),
            '		/* Lambertian reflection */',
            '		sin_theta = mc_sqrt(mcsim_random(mcsim));',
            '		cos_theta = mc_sqrt(FP_1 - sin_theta*sin_theta);',
            '',
            '		mc_sincos(mcsim_random(mcsim)*FP_2PI, &sin_fi, &cos_fi);',
            '',
            '		mcsim_set_direction_coordinates(',
            '			mcsim,',
            '			cos_fi*sin_theta,',
            '			sin_fi*sin_theta,',
            '			mc_fsign(-mcsim_direction_z(mcsim))*cos_theta',
            '		);',
            '	} else {',
            '		dbg_print("{} LambertianReflector: specular reflection");'.format(Loc),
            '		/* Specular reflection */',
            '		mcsim_reverse_direction_z(mcsim);',
            '	};',
            ''
            '	dbg_print_point3f("{} LambertianReflector reflected dir: ", '.format(Loc),
            '		mcsim_direction(mcsim));',
            '',
            '	mcsim_set_weight(mcsim, mcsim_weight(mcsim)*layout->reflectance);',
            '',
            '	return MC_REFLECTED;',
            '};',
        ))

    def __init__(self, reflectance: float = 1.0, specular: float = 0.0):
        '''
        Base class of Lambertian reflectors.

        Parameters
        ----------
        reflectance: float
            Total reflectance of the surface. Use a value between 0 and 1.
        specular: float
            Fraction of specular reflections at the surface. Use a value of 0
            for an ideal Lambertian reflector or a value of 1.0 for an ideal
            mirror (with the given reflectance). Any value between 0 and 1
            will result in surface that is a combination of a Lambertian
            and specular reflector.
        '''
        super().__init__()
        self._reflectance = 1.0
        self._specular = 0.0
        self._set_reflectance(reflectance)
        self._set_specular(specular)

    def _get_reflectance(self) -> float:
        return self._reflectance
    def _set_reflectance(self, reflectance: float):
        self._reflectance = min(max(float(reflectance), 0.0), 1.0)
    reflectance = property(_get_reflectance, _set_reflectance, None,
                           'Surface reflectance of the Lambertian reflector')

    def _get_specular(self) -> float:
        return self._specular_fraction
    def _set_specular(self, fraction: float):
        self._specular_fraction = min(max(float(fraction), 0.0), 1.0)
    specular = property(_get_specular, _set_specular, None,
                        'Specular fraction of surface reflectance. A value '
                        'of 0 indicate an ideal Lambertian reflector, and '
                        'a value of 1 a perfect mirror.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`LambertianReflector.cl_type` for a detailed
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

        target.reflectance = self._reflectance
        target.specular = self._specular

        return target

    def todict(self) -> dict:
        '''
        Export object to dict
        '''
        return {'type': self.__class__.__name__,
                'reflectance': self._reflectance, 'specular': self._specular}

    def __str__(self):
        return 'LambertianReflector(reflectance={}, specular={})'.format(
            self._reflectance, self._specular)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
