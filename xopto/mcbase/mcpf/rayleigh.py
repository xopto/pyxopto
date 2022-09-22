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

from .pfbase import PfBase, McObject
from xopto.mcbase import cltypes
import xopto.pf


class Rayleigh(PfBase):
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

        Fields
        ------
        gamma: mc_fp_t
            Parameter of the Rayleigh scattering phase function,
        a, b: mc_fp_t
            precalculated constants used to speed up the computation.
        '''
        T = mc.types
        class ClRayleigh(cltypes.Structure):
            _fields_ = [
                ('gamma', T.mc_fp_t),
                ('a', T.mc_fp_t),
                ('b', T.mc_fp_t)
            ]
        return ClRayleigh

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McPf {',
            '	mc_fp_t gamma;',
            '	mc_fp_t a;',
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
            '	dbg_print("Rayleigh scattering phase function:");',
            '	dbg_print_float(INDENT "gamma:", pf->gamma);',
            '	dbg_print_float(INDENT "a:", pf->a);',
            '	dbg_print_float(INDENT "b:", pf->b);',
            '};',
            '',
            'inline mc_fp_t mcsim_sample_pf(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t tmp, cos_theta;',
            '	mc_fp_t gamma = mcsim_current_pf(mcsim)->gamma;',
            '	mc_fp_t a = mcsim_current_pf(mcsim)->a;',
            '	mc_fp_t b = mcsim_current_pf(mcsim)->b;',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
            '',
            '	if (gamma == FP_1) {',
            '		/* isotropic scattering */',
            '		cos_theta = FP_2*mcsim_random(mcsim) - FP_1;',
            '	} else {',
            '		/* Cardano formula */',
            '		b = b*(FP_1 - FP_2*mcsim_random(mcsim))',
            '		tmp = mc_sqrt(b*b*FP_0p25 + a*a*a*FP_1d27);',
            '		cos_theta = mc_cbrt(-FP_0p5*b + tmp) + mc_cbrt(-FP_0p5*b - tmp);',
            '	);',
            '',
            '	return mc_fclip(cos_theta, -FP_1, FP_1);',
            '};'
        ))

    def __init__(self, gamma):
        '''
        The Rayleigh scattering phase function is defined as:

        .. math::

            \\frac{3}{8}((1 + 3\\gamma) + (1 - \\gamma)*\\cos^{2}(\\theta))/(1 + 2\\gamma)

        [Anthony Bucholtz, Applied Optics, Vol. 34, No. 15
        Rayleigh-scattering calculations for the terrestrial atmosphere].
        The scattering angle cosine is sampled/computed by using the Cardans
        formula to solve :math:`x^3 + ax + b = 0`, where:

        .. math::

            a &= 3 (1 + 3\\gamma)/(1 - \\gamma)

            b &= 4 (1 + 2\\gamma)/(1 - \\gamma)(1 - 2\\xi)

        Jeppe Revall Frisvad, Journal of the Optical Society of America A,
        Vol. 28, Issue 12, pp. 2436-2441 (2011),
        Importance sampling the Rayleigh phase function.

        For a special case when :math:`\\gamma = 1`, i.e. isotropic scattering,
        the scattering angle cosine becomes :math:`2\\xi - 1`, where
        :math:`\\xi` is a random number from :math:`[0, 1]`.
        For gamma != 1 the solution (:math:`\\cos(\\theta)`) is defined in
        terms of Cardano formula:

        .. math::

            x = \\cos(\\theta) &= (-b/2 + (b^{2/4} + a^{3/27})^{1/2})^{1/3} +
                                  (-b/2 - (b^{2/4} + a^{3/27})^{1/2})^{1/3}

        Parameters
        ----------
        gamma: float
            Molecular anisotropy. Scattering is isotropic for gamma = 0.
        '''
        super().__init__()
        self._gamma = 0.0
        self._precalculated = None

        self._set_gamma(gamma)

        self._recalculate()

    def _get_gamma(self) -> float:
        return self._gamma
    def _set_gamma(self, gamma: float):
        self._gamma = min(max(float(gamma), 0.0), 1.0)
        self._recalculate()
    gamma = property(_get_gamma, _set_gamma, None,
                     'Molecular anisotropy factor.')

    def _recalculate(self):
        gamma = self._gamma
        if gamma == 1.0:
            a = b = 0.0
        else:
            a = 3.0*(1.0 + 3.0*gamma)/(1.0 - gamma)
            b = 4.0*(1.0 + 2.0*gamma)/(1.0 - gamma)

        self._precalculated = (a, b)

    def pf(self) -> xopto.pf.Rayleigh:
        '''
        Returns a new instance of the related utility scattering phase function
        class that can be used to compute Legendre moments and other
        scattering phase function quantifiers.

        Returns
        -------
        pf: xopto.pf.Rayleigh
            Instance of the related utility scattering phase function.
        '''
        return xopto.pf.Rayleigh(self.gamma)

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~Rayleigh.cl_type` method
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
            instance of ClRayleigh if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        target.gamma = self._gamma
        target.a = self._precalculated[0]
        target.b = self._precalculated[1]

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'gamma': self._gamma, 'type':self.__class__.__name__}

    def __str__(self):
        return 'Rayleigh(gamma={})'.format(self._gamma)
