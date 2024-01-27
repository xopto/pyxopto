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

from .pfbase import PfBase, cltypes, McObject
from . import Hg
import xopto.pf


class Hg2(PfBase):
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
        hg1: mc_pf_hg_t
            The first Henyey-Greenstein scattering phase function,
        hg2: mc_pf_hg_t
            The second Henyey-Greenstein scattering phase function,
        b: mc_fp_t
            Relative contribution of the second Henyey-Greenstein
            scattering phase function. Must be a value from [0, 1].
        '''
        T = mc.types
        T_HG = Hg.cl_type(mc)
        class ClHg2(cltypes.Structure):
            _fields_ = [
                ('hg_1', T_HG),
                ('hg_2', T_HG),
                ('b',    T.mc_fp_t),
            ]


        return ClHg2

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return  '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES mc_pf_hg_t {'
            '	mc_fp_t g;',
            '};',
            '',
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	struct mc_pf_hg_t hg_1;',
            '	struct mc_pf_hg_t hg_2;',
            '	mc_fp_t b;',
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
            '	dbg_print("Hg2 scattering phase function:");',
            '	dbg_print(INDENT "The first HG term:");',
            '	dbg_print_float(INDENT "g1:", pf->hg_1.g);',
            '',
            '	dbg_print(INDENT "The second HG term:");',
            '	dbg_print_float(INDENT "g2:", pf->hg_2.g);',
            '};',
            '',
            'inline mc_fp_t mcsim_pf_sample_angles(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t k, cos_theta;',
            '	mc_fp_t b = mcsim_current_pf(mcsim)->b;',
            '',
            '	mc_fp_t g;',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
            '',
            '	if (mcsim_random(mcsim) >= b) {',
            '		g = mcsim_current_pf(mcsim)->hg_1.g;',
            '		k = mc_fdiv(',
            '			(FP_1 - g*g),',
            '			(FP_1 + g*(FP_2*mcsim_random(mcsim) - FP_1))',
            '		);',
            '		cos_theta = mc_fdiv((FP_1 + g*g - k*k), FP_2*g);'
            '	} else {',
            '		g = mcsim_current_pf(mcsim)->hg_2.g;',
            '		k = mc_fdiv(',
            '			(FP_1 - g*g),',
            '			(FP_1 + g*(FP_2*mcsim_random(mcsim) - FP_1))',
            '		);',
            '		cos_theta = mc_fdiv((FP_1 + g*g - k*k), FP_2*g);'
            '	};',
            '',
            '	if(g == FP_0)',
            '		cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '',
            '	return mc_fmax(mc_fmin(cos_theta, FP_1), -FP_1);',
            '};'
        ))

    def __init__(self, g1: float, g2: float, b: float):
        '''
        Henyey-Greenstein (HG) scattering phase function object constructor.

        Parameters
        ----------
        g1: float
            Anisotropy factor of the first HG scattering phase function.
        g2: float
            Anisotropy factor of the second HG scattering phase function.
        b: float
            A value between 0 and 1 that determines the relative contribution of
            the second HG scattering scattering phase function.
        '''
        super().__init__()
        self._hg1 = Hg(g1)
        self._hg2 = Hg(g2)
        self._b = 0.0

        self._set_g1(g1)
        self._set_g2(g2)
        self._set_b(b)

    def _get_g1(self) -> float:
        return self._hg1.g
    def _set_g1(self, g: float):
        self._hg1.g = min(max(float(g), 0.0), 1.0)
    g1 = property(_get_g1, _set_g1, None, 'Anisotropy factor of the first HG.')

    def _get_g2(self) -> float:
        return self._hg2.g
    def _set_g2(self, g: float):
        self._hg2.g = min(max(float(g), -1.0), 0.0)
    g2 = property(_get_g2, _set_g2, None, 'Anisotropy factor of the second HG.')

    def _get_b(self) -> float:
        return self._b
    def _set_b(self, b: float):
        self._b = min(max(float(b), 0.0), 1.0)
    b = property(_get_b, _set_b, None,
                 'Relative contribution of the second HG term.')

    def pf(self) -> xopto.pf.Hg:
        '''
        Returns a new instance of the related utility scattering phase function
        class that can be used to compute Legendre moments and other
        scattering phase function quantifiers.

        Returns
        -------
        pf: xopto.pf.Hg2
            Instance of the related utility scattering phase function.
        '''
        return xopto.pf.Hg2(self.g1, self.g2, self.b)

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the an OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~Hg2.cl_type` method
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
            instance of ClHg2 if the input argument target is None.
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        self._hg1.cl_pack(mc, target.hg_1)
        self._hg2.cl_pack(mc, target.hg_2)
        target.b = self._b


        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'g1': self._hg1.g,
                'g2': self._hg2.g,
                'b': self._b,
                'type': self.__class__.__name__}

    def __str__(self):
        return 'Hg2(g1={}, g2={}, b={})'.format(
            self._hg1.g, self._hg2.g, self._b)
