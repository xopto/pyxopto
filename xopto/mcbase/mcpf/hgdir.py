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
from .pfbase import PfBase, cltypes, McObject
from .. import mcoptions
import xopto.pf

import numpy as np


class HgDir(PfBase):
    @staticmethod
    def cl_type(mc: McObject) -> cltypes.Structure:
        '''
        Returns an OpenCL structure that can be passed to the Monte carlo
        simulator.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents the scattering phase function in
            the Monte Carlo kernel.

        Structure fields
        ----------------
        direction: mc_point3f_t
            Direction around which to scatter the packet.
        g: mc_fp_t
            Parameter of the Henyey-Greenstein scattering phase function,
        p: mc_fp_t
            Probability of directional scattering.
        '''
        T = mc.types
        class ClHgDir(cltypes.Structure):
            _fields_ = [
                ('direction', T.mc_point3f_t),
                ('g', T.mc_fp_t),
                ('p', T.mc_fp_t)
            ]

        return ClHgDir

    @staticmethod
    def cl_options(mc: McObject) -> mcoptions.RawOptions:
        return [('MC_PF_SAMPLE_DIRECTION', True)]

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return  '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	mc_point3f_t direction;',
            '	mc_fp_t g;',
            '	mc_fp_t p;',
            '};',
            '',
            'void dbg_print_pf(const McPf *pf);',
        ))

    @staticmethod
    def cl_implementation(mc: McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return '\n'.join((
            'void dbg_print_pf(const McPf *pf) {',
            '	dbg_print("HgDir scattering phase function:");',
            '	dbg_print_point3f(INDENT "direction:", pf->direction);',
            '	dbg_print_float(INDENT "g:", pf->g);',
            '	dbg_print_float(INDENT "p", pf->p);',
            '};',
            '',
            'inline void mcsim_pf_sample_dir(McSim *mcsim, mc_point3f_t *out_dir){',
            '	mc_fp_t k, cos_theta;',
            '	mc_fp_t g = mcsim_current_pf(mcsim)->g;',
            '',
            '	k = mc_fdiv(',
            '		(FP_1 - g*g),',
            '		(FP_1 + g*(FP_2*mcsim_random(mcsim) - FP_1))',
            '	);',
            '	cos_theta = mc_fdiv((FP_1 + g*g - k*k), FP_2*g);'
            '',
            '	if(g == FP_0)',
            '		cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '',
            '	cos_theta = mc_fmax(mc_fmin(cos_theta, FP_1), -FP_1);',
            '',
            '	*out_dir = (mcsim_random(mcsim) < mcsim_current_pf(mcsim)->p) ?',
            '		mcsim_current_pf(mcsim)->direction : *mcsim_direction(mcsim);',
            '	cos_theta = (mc_dot_point3f(mcsim_direction(mcsim), out_dir) < FP_0) ?',
            '		-cos_theta : cos_theta;',
            '	scatter_direction(out_dir, cos_theta, FP_2PI*mcsim_random(mcsim));',
            '};'
        ))

    def __init__(self, g: float,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0),
                 p: float = 1.0):
        '''
        Directional Henyey-Greenstein scattering phase function object constructor.

        Parameters
        ----------
        g: float
            Anisotropy factor. Assigned to the first user-defined
            scattering phase function parameter.
        direction: Tuple[float, float, float]
            Direction around which this scattering phase function scatters light.
        p: float
            Probability of directional scattering.
        '''
        super().__init__()
        self._g = 0.0
        self._direction = np.array((0.0, 0.0, 1.0))
        self._set_direction(direction)
        self._set_g(g)
        self._p = 0.0
        self._set_p(p)

    def _get_g(self) -> float:
        return self._g
    def _set_g(self, g: float):
        self._g = min(max(float(g), -1.0), 1.0)
    g = property(_get_g, _set_g, None, 'Anisotropy factor.')

    def _get_p(self) -> float:
        return self._p
    def _set_p(self, p: float):
        self._p = min(max(float(p), 0.0), 1.0)
    p = property(_get_p, _set_p, None, 'Probability of directional scattering.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('The norm/length of the scattering direction '
                             'vector must not be 0!')
        self._direction *= 1.0/norm
    direction = property(_get_direction, _set_direction, None,
                        'Scattering direction.')

    def pf(self) -> xopto.pf.Hg:
        '''
        Returns a new instance of the related utility scattering phase function
        class that can be used to compute Legendre moments and other
        scattering phase function quantifiers.

        Returns
        -------
        pf: xopto.pf.Hg
            Instance of the related utility scattering phase function.
        '''
        return None

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the an OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~Hg.cl_type` method
        for a detailed list of fields.

        Parameters
        ----------
        mc: McObject
            Simulator instance.
        target: cltypes.Structure
            Target OpenCL structure for packing.

        Returns
        -------
        target: cltypes.Structure
            Target structure received as an input argument or a new
            instance of ClHg if the input argument target is None.
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        target.g = self._g
        target.p = self._p
        target.direction.fromarray(self._direction)

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'g': self._g, 'direction': tuple(self._direction), 'p': self._p,
                'type': self.__class__.__name__}

    def __str__(self):
        return 'HgDir(g={}, direction=({}, {}, {}), p={})'.format(
            self._g, *self._direction, self._p)
