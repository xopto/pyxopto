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


class Pc(PfBase):
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
            Power of the cosine,
        '''
        T = mc.types
        class ClPc(cltypes.Structure):
            _fields_ = [('n', T.mc_fp_t)]

        return ClPc


    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	mc_fp_t n;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return '\n'.join((
            'void dbg_print_pf(const McPf *pf) {',
            '	dbg_print("Pc scattering phase function:");',
            '	dbg_print_float(INDENT "n:", pf->n);',
            '};',
            '',
            'inline mc_fp_t mcsim_sample_pf(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t tmp, cos_theta;',
            '	mc_fp_t n = mcsim_current_pf(mcsim)->n;',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
            '',
            '	cos_theta = FP_2*mc_pow(',
            '		mcsim_random(mcsim),',
            '		mc_fdiv(FP_1, n + FP_1)',
            '	) - FP_1;',
            '',
            '	return mc_fclip(cos_theta, -FP_1, FP_1);',
            '};'
        ))

    def __init__(self, n: float):
        '''
        Power of cosines (PC) scattering phase function.

        .. math::

            Pc(\\cos(\\theta)) &= \\frac{n + 1}{2^{n + 1}}(1 + \\cos(\\theta))^{n}

        Parameters
        ----------
        n: float >= 0
            Parameter of the power of cosine scattering phase function.

        Note
        ----
        The cumulative density function of PC follows:

        .. math::

            CDF(\\cos(\\theta)) &= \\frac{(1 + \\cos(\\theta))^{n + 1}}{2^{n + 1}}

        For a given uniformly distributed random number :math`\\xi \\in [0, 1]`,
        the scattering angle cosine can be computed as:

        .. math::

            \\cos(\\theta) &= 2 e^{\\frac{1}{n + 1}} - 1

        Since this expression becomes singular for n=-1, only positive values
        of :math:`n` are allowed.
        '''
        super().__init__()
        self._n = 0.0
        self._set_n(n)

    def _get_n(self) -> float:
        return self._n
    def _set_n(self, n: float):
        self._n = float(n)
    n = property(_get_n, _set_n, None, 'Power of cosine.')

    def pf(self) -> xopto.pf.Pc:
        '''
        Returns a new instance of the related utility scattering phase function
        class that can be used to compute Legendre moments and other
        scattering phase function quantifiers.

        Returns:
        -------
        pf: xopto.pf.Pc
            Instance of the related utility scattering phase function.
        '''
        return xopto.pf.PPc(self.n)

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~Pc.cl_type` method
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
            instance of ClPc if the input argument target is None.
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        target.n = self._n

        return target

    def todict(self):
        '''
        Export object to a dict.
        '''
        return {'n': self._n, 'type': self.__class__.__name__}

    def __str__(self):
        return 'Pc(n={})'.format(self._n)
