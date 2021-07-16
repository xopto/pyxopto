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

class MHg(PfBase):
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
        g: mc_fp_t
            Parameter of the Henyey-Greenstein scattering phase function.
        beta: mc_fp_t
            Contribution of the Henyey-Greenstein scattering phase function.
        '''
        T = mc.types
        class ClMHg(cltypes.Structure):
            _fields_ = [
                ('g', T.mc_fp_t),
                ('beta', T.mc_fp_t),
            ]

        return ClMHg

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return  '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	mc_fp_t g;',
            '	mc_fp_t beta;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return  '\n'.join((
            'void dbg_print_pf(const McPf *pf) {',
            '	dbg_print("MHg scattering phase function:");',
            '	dbg_print_float(INDENT "g:", pf->g);',
            '	dbg_print_float(INDENT "beta:", pf->beta);',
            '};',
            '',
            'inline mc_fp_t mcsim_sample_pf(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t k, cos_theta;',
            '	mc_fp_t g = mcsim_current_pf(mcsim)->g;',
            '	mc_fp_t beta = mcsim_current_pf(mcsim)->beta;',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
            '',
            '	if (mcsim_random(mcsim) <= beta) {',
            '		k = mc_fdiv(',
            '			(FP_1 - g*g),',
            '			(FP_1 + g*(FP_2*mcsim_random(mcsim) - FP_1))',
            '		);',
            '		cos_theta = mc_fdiv((FP_1 + g*g - k*k), FP_2*g);'
            '',
            '		if(g == FP_0)',
            '			cos_theta = FP_1 - FP_2*mcsim_random(mcsim);',
            '	} else {',
            '		cos_theta = mc_cbrt(FP_2*mcsim_random(mcsim) - FP_1);',
            '	};',
            '',
            '	/* printf("cos_theta: %.3f\\n", cos_theta); */'
            '	return mc_fclip(cos_theta, -FP_1, FP_1);',
            '};'
        ))

    def __init__(self, g: float, b: float):
        '''
        Modified Henyey-Greenstein scattering phase function constructor.

        b*HG(g) + \\frac{2}{3} (1 - b) \\cos^{2}(\\theta)

        Parameters
        ----------
        g: float
            Anisotropy factor. Assigned to the first user-defined
            scattering phase function parameter.
        b: float
            Contribution of the Henyey-Greenstein scattering phase function.
        '''
        super().__init__()
        self._g = self._b = 0.0
        self._set_g(g)
        self._set_b(b)

    def _get_g(self) -> float:
        return self._g
    def _set_g(self, g: float):
        self._g = min(max(float(g), -1.0), 1.0)
    g = property(_get_g, _set_g, None, 'Anisotropy factor.')

    def _get_b(self) -> float:
        return self._b
    def _set_b(self, b: float):
        self._b = min(max(float(b), 0.0), 1.0)
    b = property(_get_b, _set_b, None, 'Contribution of the Henyey-Greenstein '\
                                       'scattering phase function.')

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills an OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~MHg.cl_type` method
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
            instance of ClMHg if the input argument target is None.
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        target.g = self._g
        target.beta = self._b

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'g': self._g, 'b': self._b, 'type': self.__class__.__name__}

    def __str__(self):
        return 'MHg(g={}, b={})'.format(self._g, self._b)
