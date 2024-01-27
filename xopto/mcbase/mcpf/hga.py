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

from .pfbase import PfBase, cltypes, McObject
from xopto.pf import hga


class Hga(PfBase):
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
        g: mc_matrix3f_t
            Tensor (3 x 3) of the Henyey-Greenstein scattering phase function,
        '''
        T = mc.types
        class ClHga(cltypes.Structure):
            _fields_ = [('g', T.mc_matrix3f_t)]

        return ClHga

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return  '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	mc_matrix3f_t g;',
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
            '	dbg_print_matrix3f("Hga scattering phase function:", &pf->g);',
            '};',
            '',
            'inline mc_fp_t mcsim_pf_sample_angles(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t k, cos_theta;',
            '	mc_fp_t g = tensor3f_project(',
            '		&mcsim_current_pf(mcsim)->g, mcsim_direction(mcsim));',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
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
            '	return mc_fmax(mc_fmin(cos_theta, FP_1), -FP_1);',
            '};'
        ))

    def __init__(self, g: float or Tuple[float, float, float] or np.ndarray):
        '''
        Henyey-Greenstein scattering phase function object constructor.

        Parameters
        ----------
        g: float or Tuple[float, float, float] or np.ndarray
            Anisotropy tensor.
        '''
        super().__init__()
        self._g = np.zeros((3, 3))
        self._set_g(g)

    def _get_g(self) -> np.ndarray:
        return self._g
    def _set_g(self, g: float or np.ndarray):
        if isinstance(g, (float, int)):
            g = min(max(g, -1.0), 1.0)
            self._g[0, 0] = g
            self._g[1, 1] = g
            self._g[2, 2] = g
        else:
            g = np.asarray(g, dtype=float)
            if g.size == 3:
                g = np.clip(g, -1.0, 1.0)
                self._g[0, 0] = g[0] 
                self._g[1, 1] = g[1] 
                self._g[2, 2] = g[2] 
            else:
                self._g[:] = min(max(g, -1.0), 1.0)

    g = property(_get_g, _set_g, None, 'Anisotropy tensor.')

    def pf(self) -> hga.Hga:
        '''
        Returns a new instance of the related utility scattering phase function
        class that can be used to compute Legendre moments and other
        scattering phase function quantifiers.

        Returns
        -------
        pf: xopto.pf.hga.Hga
            Instance of the related utility scattering phase function.
        '''
        return hga.Hga(self._g)

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the OpenCL Structure (target) with the data required by the
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

        target.g.fromarray(self._g)

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'g': self._g, 'type': self.__class__.__name__}

    def __str__(self):
        return 'Hga(g={})'.format(self._g)
