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
from . import Gk
import xopto.pf


class Gk2(PfBase):
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
        gk_1: mc_pf_gk_t
            Parameters of the first Gegenbauer kernel scattering phase function.
        gk_2: mc_pf_gk_t
            Parameters of the second Gegenbauer kernel scattering phase function.
        b: mc_fp_t
            Relative contribution of the second Gegenbauer kernel scattering
            phase function.
        '''
        T = mc.types
        T_GK = Gk.cl_type(mc)
        class ClGk2(cltypes.Structure):
            _fields_ = [
                ('gk_1', T_GK),
                ('gk_2', T_GK),
                ('b',    T.mc_fp_t),
            ]

        return ClGk2

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return  '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES mc_pf_gk_t {'
            '	mc_fp_t g;',
            '	mc_fp_t a;',
            '	mc_fp_t inv_a;',
            '	mc_fp_t a1;',
            '	mc_fp_t a2;',
            '};',
            '',
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	struct mc_pf_gk_t gk_1;',
            '	struct mc_pf_gk_t gk_2;',
            '	mc_fp_t b;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return '\n'.join((
            'void dbg_print_pf(const McPf *pf) {',
            '	dbg_print("Gk2 scattering phase function:");',
            '	dbg_print(INDENT "The first Gk term:");',
            '	dbg_print_float(INDENT INDENT "gk_1.g:",     pf->gk_1.g_1);',
            '	dbg_print_float(INDENT INDENT "gk_1.a:",     pf->gk_1.a_1);',
            '	dbg_print_float(INDENT INDENT "gk_1.inv_a:", pf->gk_1.inv_a_1);',
            '	dbg_print_float(INDENT INDENT "gk_1.a1:",    pf->gk_1.a1_1);',
            '	dbg_print_float(INDENT INDENT "gk_1.a2:",    pf->gk_1.a2_1);',
            '',
            '	dbg_print(INDENT "The second Gk term:");',
            '	dbg_print_float(INDENT INDENT "gk_2.g:",     pf->gk_2.g_1);',
            '	dbg_print_float(INDENT INDENT "gk_2.a:",     pf->gk_2.a_1);',
            '	dbg_print_float(INDENT INDENT "gk_2.inv_a:", pf->gk_2.inv_a_1);',
            '	dbg_print_float(INDENT INDENT "gk_2.a1:",    pf->gk_2.a1_1);',
            '	dbg_print_float(INDENT INDENT "gk_2.a2:",    pf->gk_2.a2_1);',
            '',
            '	dbg_print(INDENT "Relative contribution of the second Gk term:");',
            '	dbg_print_float(INDENT INDENT "b:", pf->b);',
            '};',
            '',
            'inline mc_fp_t mcsim_sample_pf(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t tmp, cos_theta;',
            '	mc_fp_t b = mcsim_current_pf(mcsim)->b;',
            '',
            '	mc_fp_t g;',
            '	mc_fp_t a;',
            '	mc_fp_t inv_a;',
            '	mc_fp_t a1;',
            '	mc_fp_t a2;',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
            '',
            '	if (mcsim_random(mcsim) >= b) {',
            '		g = mcsim_current_pf(mcsim)->gk_1.g;',
            '		a = mcsim_current_pf(mcsim)->gk_1.a;',
            '		inv_a = mcsim_current_pf(mcsim)->gk_1.inv_a;',
            '		a1 = mcsim_current_pf(mcsim)->gk_1.a1;',
            '		a2 = mcsim_current_pf(mcsim)->gk_1.a2;',
            '	} else {',
            '		g = mcsim_current_pf(mcsim)->gk_2.g;',
            '		a = mcsim_current_pf(mcsim)->gk_2.a;',
            '		inv_a = mcsim_current_pf(mcsim)->gk_2.inv_a;',
            '		a1 = mcsim_current_pf(mcsim)->gk_2.a1;',
            '		a2 = mcsim_current_pf(mcsim)->gk_2.a2;',
            '	};',
            '',
            '	if (g == FP_0) {',
            '		cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '	} else if(a == FP_0) {',
            '		cos_theta = a1 + mc_pow(',
            '			mc_fdiv(FP_1 - g, FP_1 + g),',
            '			FP_2*mcsim_random(mcsim)',
            '		)*a2;',
            '	} else {',
            '		tmp = a1*mcsim_random(mcsim) + a2;',
            '		tmp = FP_1 + g*g - mc_pow(tmp, -inv_a);',
            '		cos_theta = mc_fdiv(tmp, FP_2*g);',
            '	};',
            '',
            '	return mc_fclip(cos_theta, -FP_1, FP_1);',
            '};'
        ))

    def __init__(self, g1: float, a1: float, g2: float, a2: float, b: float):
        '''
        Gegenbauer kernel scattering phase function constructor.

        Parameters
        ----------
        g1: float
            Parameter of the first Gegenbauer kernel scattering scattering phase
            function.
            :math:`0 <= g <= 1`
        a1: float
            Parameter of the first Gegenbauer kernel scattering phase function.
            :math:`a > - 1/2`
            A value of 0.5 produces the Henyey-Greenstein scattering
            phase function.
        g2: float
            Parameter of the second Gegenbauer kernel scattering scattering phase
            function.
            :math:`-1 <= g < 0`
        a2: float
            Parameter of the second Gegenbauer kernel scattering phase function.
            :math:`a > - 1/2`
            A value of 0.5 produces the Henyey-Greenstein scattering
            phase function.
        b: float
            A value between 0 and 1 that determines the relative contribution of
            the second Gegenbauer kernel scattering scattering phase function.
        '''
        super().__init__()
        self._gk1 = Gk(g1, a1)
        self._gk2 = Gk(g2, a2)
        self._b = 0.0

        self._set_b(b)

    def _get_g1(self) -> float:
        return self._gk1.g
    def _set_g1(self, g: float):
        self._gk1.g = min(max(float(g), 0.0), 1.0)
    g1 = property(_get_g1, _set_g1, None,
                  'Parameter g (anisotropy factor when a=0.5) '
                  'of the first GK term.')

    def _get_a1(self) -> float:
        return self._gk1.a
    def _set_a1(self, a: float):
        self._gk1.a = max(float(a), -0.5)
    a1 = property(_get_a1, _set_a1, None,
                  'Parameter alpha of the first GK term.')

    def _get_g2(self) -> float:
        return self._gk2.g
    def _set_g2(self, g: float):
        self._gk2.g = min(max(float(g), -1.0), 0.0)
    g2 = property(_get_g2, _set_g2, None,
                  'Parameter g (anisotropy factor when a=0.5) '
                  'of the second GK term.')

    def _get_a2(self) -> float:
        return self._gk2.a
    def _set_a2(self, a: float):
        self._gk2.a = max(float(a), -0.5)
    a2 = property(_get_a2, _set_a2, None,
                  'Parameter alpha of the second GK term.')

    def _get_b(self) -> float:
        return self._b
    def _set_b(self, b: float):
        self._b = min(max(float(b), 0.0), 1.0)
    b = property(_get_b, _set_b, None,
                 'Relative contribution of the second GK term.')

    def pf(self) -> xopto.pf.Gk2:
        '''
        Returns a new instance of the related utility scattering phase function
        class that can be used to compute Legendre moments and other
        scattering phase function quantifiers.

        Returns
        -------
        pf: xopto.pf.Gk2
            Instance of the related utility scattering phase function.
        '''
        return xopto.pf.Gk2(self._gk1.g, self._gk1.a,
                            self._gk2.g, self._gk2.a, self.b)

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the an OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~Gk2.cl_type` method
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
            instance of ClGk if the input argument target is None.
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        self._gk1.cl_pack(mc, target.gk_1)
        self._gk2.cl_pack(mc, target.gk_2)
        target.b = self._b

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'g1': self._gk1.g, 'a1': self._gk1.a,
                'g2': self._gk2.g, 'a2': self._gk2.a,
                'b': self._b, 'type': self.__class__.__name__}

    def __str__(self):
        return 'Gk2(g1={}, a1={}, g2={}, a2={}, b={})'.format(
            self._gk1.g, self._gk1.a, self._gk2.g, self._gk2.a, self._b)
