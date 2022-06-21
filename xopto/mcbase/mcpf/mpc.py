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
import xopto.pf


class MPc(PfBase):
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
        n: mc_fp_t
            Power of the cosine.

        beta: mc_fp_t
            Contribution of the power of cosine scattering phase function.
        '''
        T = mc.types
        class ClMPc(cltypes.Structure):
            _fields_ = [
                ('n', T.mc_fp_t),
                ('beta', T.mc_fp_t)
            ]

        return ClMPc

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	mc_fp_t n;',
            '	mc_fp_t beta;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return '\n'.join((
            'void dbg_print_pf(const McPf *pf) {',
            '	dbg_print("MPc scattering phase function:");',
            '	dbg_print_float(INDENT "n:", pf->n);',
            '	dbg_print_float(INDENT "beta:", pf->beta);',
            '};',
            '',
            'inline mc_fp_t mcsim_sample_pf(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t tmp, cos_theta;',
            '	mc_fp_t n = mcsim_current_pf(mcsim)->n;',
            '	mc_fp_t beta = mcsim_current_pf(mcsim)->beta;',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
            '',
            '	if (mcsim_random(mcsim) <= beta){',
            '		cos_theta = FP_2*mc_pow(',
            '			mcsim_random(mcsim),',
            '			mc_fdiv(FP_1, n + FP_1)',
            '		) - FP_1;',
            '	} else {',
            '		/* Isotropic rayleigh scattering phase fun.: 3.0/(2.0)*costheta**2 */',
            '		cos_theta = mc_cbrt(FP_2*mcsim_random(mcsim) - FP_1);',
            '	};',
            '',
            '	return mc_fclip(cos_theta, -FP_1, FP_1);',
            '};'
        ))

    def __init__(self, n: float, b: float):
        '''
        Modified power of cosines scattering phase function.

        .. math::

            Pc(\\cos(\\theta)) &= \\frac{n + 1}{2^{n + 1}}(1 + \\cos(\\theta))^{n}

            MPc(\\cos(\\theta)) &= Pc(\\cos(\\theta))b + \\frac{3}{2}\\cos^{2}(\\theta)(1 - b)

        Parameters
        ----------
        n: float
            Parameter of the power of cosine scattering phase function.
        b: float
            Contribution of the Rayleigh scattering component defined as
            (1 - b).
        '''
        super().__init__()
        self._n = self._b = 0.0
        self._set_n(n)
        self._set_b(b)

    def _get_n(self) -> float:
        return self._n
    def _set_n(self, n: float):
        self._n = float(n)
    n = property(_get_n, _set_n, None, 'Power of cosine.')

    def _get_b(self) -> float:
        return self._b
    def _set_b(self, b: float):
        self._b = min(max(float(b), 0.0), 1.0)
    b = property(_get_b, _set_b, None, 'Contribution of the power of cosine '\
                                       'scattering phase function.')

    def pf(self) -> xopto.pf.MPc:
        '''
        Returns a new instance of the related utility scattering phase function
        class that can be used to compute Legendre moments and other
        scattering phase function quantifiers.

        Returns:
        -------
        pf: xopto.pf.MPc
            Instance of the related utility scattering phase function.
        '''
        return xopto.pf.MPc(self.n, self.b)

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the an OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~MPc.cl_type` method
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
            instance of ClMPc if the input argument target is None.
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()


        target.n = self._n
        target.beta = self._b

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'n': self._n, 'b':self._b, 'type': self.__class__.__name__}

    def __str__(self):
        return 'MPc(n={}, b={})'.format(self._n, self._b)
