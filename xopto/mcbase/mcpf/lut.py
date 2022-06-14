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

from .. import mcoptions
from .pfbase import PfBase, cltypes, McObject

from xopto import pf


class Lut(PfBase):
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

            The returned structure type implements the following fields:

            - a: mc_fp_t
                Parameter of the nonlinear lookup table index transformation.
            - b: mc_fp_t
                Parameter of the nonlinear lookup table index transformation.
            - c: mc_fp_t
                Parameter of the nonlinear lookup table index transformation.
            - offset: mc_int_t
                Location of the lookup table as an offset from the start of
                the lookup table buffer.
            - size: mc_int_t
                Number of entries in the lookup table.

        Note
        ----
        Note that this implementation of the lookup
        table uses the following nonlinear transformation to compute
        the lookup table index from a uniform random number
        :math:`\\xi \\in [0, 1]`.

        .. math::

            index &= \\frac{1}{2}(\\frac{a}/{\\xi - c} - b + 1)N,

        where :math:`a, b, c` are the parameters of the nonlinear index
        transformation and :math:`N` is the size of the lookup table
        '''
        T = mc.types
        class ClLut(cltypes.Structure):
            _fields_ = [
                ('a', T.mc_fp_t),
                ('b', T.mc_fp_t),
                ('c', T.mc_fp_t),
                ('offset', T.mc_size_t),
                ('size', T.mc_size_t)
            ]

        return ClLut

    @staticmethod
    def cl_options(mc: McObject) -> mcoptions.RawOptions:
        return [('MC_USE_FP_LUT', True)]

    @staticmethod
    def cl_declaration(mc: McObject) -> str:
        '''
        OpenCL declarations of the scattering phase function.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McPf{',
            '	mc_fp_t a;',
            '	mc_fp_t b;',
            '	mc_fp_t c;',
            '	mc_size_t offset;',
            '	mc_size_t size;',
            '};'
        ))

    @staticmethod
    def cl_implementation(mc: McObject) -> str:
        '''
        OpenCL implementation of the scattering phase function.
        '''
        return '\n'.join((
            'void dbg_print_pf(const McPf *pf) {',
            '	dbg_print("Lut scattering phase function:");',
            '	dbg_print_float(INDENT "a:",  pf->a);',
            '	dbg_print_float(INDENT "b:",  pf->b);',
            '	dbg_print_float(INDENT "c:",  pf->c);',
            '	dbg_print_size_t(INDENT "offset:",  pf->offset);',
            '	dbg_print_size_t(INDENT "size:",  pf->size);',
            '};',
            '',
            'inline mc_fp_t mcsim_sample_pf(McSim *mcsim, mc_fp_t *azimuth){',
            '	mc_fp_t fp_index, fp_index_floor, d;',
            '   size_t index;',
            '	mc_fp_t a = mcsim_current_pf(mcsim)->a;',
            '	mc_fp_t b = mcsim_current_pf(mcsim)->b;',
            '	mc_fp_t c = mcsim_current_pf(mcsim)->c;',
            '	size_t offset = mcsim_current_pf(mcsim)->offset;',
            '	size_t lut_size = mcsim_current_pf(mcsim)->size - 1;',
            '',
            '	*azimuth = FP_2PI*mcsim_random(mcsim);',
            '',
            '	/* calculate the floating point lookup table index */',
            '	fp_index = (mc_fdiv(a, mcsim_random(mcsim) - c) - b + FP_1)*',
            '		lut_size*FP_0p5;',
            '',
            '	/* calculate the weight of second point */',
            '	fp_index_floor = mc_floor(fp_index);',
            '	d = fp_index - fp_index_floor;',
            '',
            '	/* calculate the integer index of the first point in the lookup table */',
            '	index = mc_int(fp_index_floor);',
            '',
            '	/* now do the linear interpolation - dont forget the lookup table offset */ ',
            '	mc_fp_t cos_theta = mcsim_pf_lut_array(mcsim)[offset + index]*(FP_1 - d) +',
            '		mcsim_pf_lut_array(mcsim)[offset + mc_min(index + 1, lut_size)]*d;',
            '	//dbg_print_float("cos_theta:", cos_theta);',
            '	//dbg_print_float("lut[0]:", mcsim_pf_lut_array(mcsim)[offset]);',
            '	//dbg_print_float("lut[1]:", mcsim_pf_lut_array(mcsim)[offset + 1]);',
            '	//dbg_print_float("lut[2]:", mcsim_pf_lut_array(mcsim)[offset + 2]);',
            '	//dbg_print_float("lut[3]:", mcsim_pf_lut_array(mcsim)[offset + 3]);',
            '	//dbg_print_float("lut[-2]:", mcsim_pf_lut_array(mcsim)[offset + lut_size - 1]);',
            '	//dbg_print_float("lut[-1]:", mcsim_pf_lut_array(mcsim)[offset + lut_size]);',
            '	return cos_theta;',
            '};'
        ))

    def __init__(self, params: list, lut: np.ndarray):
        '''
        Lookup table-based scattering phase function constructor.
        Lookup table is used for efficient computation of deflection
        angle cosines.
        A nonlinear transformation is used to compute linear lookup table
        (floating point) index from a uniform random number F:

        .. code-block:: python

            index = 0.5*(params[0]/(F - params[2]) - params[1] + 1.0)*data.size

        Parameters
        ----------
        params: list, tuple or ndarray vector of 3 float
            Nonlinear index transformation function parameters provided as
            an array of size 3.
        lut: np.ndarray vector
            Lookup table data.

        Note
        ----
            The offset of the first lookup table entry is stored in
            the 45th user-defined scattering phase function parameter.
            The lookup table size is stored in the 4th user-defined
            scattering phase function parameter.
        '''
        super().__init__()
        self._offset = 0
        self._lut = np.asarray(lut, dtype=np.float64)
        self._params = np.zeros((3,))
        self._params[:] = params

    def _get_params(self) -> list:
        return self._params
    def _set_params(self, params: list):
        self._params[:] = params
    params = property(_get_params, _set_params, None,
                      'Lookup table index transformation parameters.')

    def _get_lut(self) -> np.ndarray:
        return self._lut
    def _set_lut(self, lut: np.ndarray):
        self._lut = lut
    lut = property(_get_lut, _set_lut, None, 'Lookup table object or None.')

    def _get_size(self) -> int:
        return self._lut.size
    size = property(_get_size, None, None, 'Lookup table size.')

    def _get_offset(self) -> int:
        return self._offset
    def _set_offset(self, offset: int):
        self._offset = int(offset)
    offset = property(_get_offset, _set_offset, None,
                      'Offset of the first element in the lookup table.')

    def cl_pack(self, mc: McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the an OpenCL cltypes.Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~Lut.cl_type` method
        for a detailed list of fields.

        Note
        ----
        Calls the :py:meth:`Mc.pack_fp_lut` method to pack the lookup table
        array data.

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
            instance of ClLut if the input argument target is None.
        '''
        if target is None:
            target_type = self.fetch_cl_type(mc)
            target = target_type()

        lut_entry = mc.append_r_lut(self._lut)
        self.offset = lut_entry.offset

        target.a = self._params[0]
        target.b = self._params[1]
        target.c = self._params[2]
        target.offset = self.offset
        target.size = self._lut.size

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'params': self._params, 'lut':self._lut,
               'type': self.__class__.__name__}

    def __str__(self):
        return 'Lut(params=({}, {}, {}), data={})'.format(
            self._params[0], self._params[1], self._params[2], self._lut)


class LutEx(Lut):
    def __init__(self, pftype: PfBase or str, pfargs: list,
                 lutsize: int = 2000, **kwargs):
        '''
        Creates a Monte Carlo lookup table-based scattering phase function
        instance from a standard scattering phase function instance
        (:py:class:`xopto.pf.PfBase` subclass).

        Parameters
        ----------
        pftype: PfBase subclass or str
            Scattering phase function type. Use one defined in the
            :py:mod:`xopto.pf` module.
        pfargs: list, tuple
            A list of scattering phase function parameters passed to the pftype.
        lutsize: int
            Number of elements in the lookup table.
        kwargs: dict
            Additional keyword arguments passed to the lut method of the
            scattering phase function as
            :code:`pftype(*pfargs).lut(lutsize, **kwargs)`
        '''
        if isinstance(pftype, str):
            pftype = getattr(pf, pftype)
            if pftype is None:
                raise ValueError(
                    'Scattering Phase function type "{}" not found '\
                    'in the xopto.pf.pf module!'.format(pftype)
                )

        self._pfargs = pfargs
        self._kwargs = kwargs
        self._pf_obj = pftype(*pfargs)
        lut_args = self._pf_obj.mclut(lutsize, **kwargs)
        super().__init__(*lut_args)

    def _get_pfargs(self) -> tuple:
        return tuple(self._pfargs)
    pfargs = property(_get_pfargs, None, None,
                      'A tuple of scattering phase function arguments.')

    def _get_pf(self) -> PfBase:
        return self._pf_obj
    pf = property(_get_pf, None, None,
                  'Scattering phase function object that was used '\
                  'to create a lookup table.')

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'pftype': self._pf_obj.__class__.__name__, 
                'pfargs':self._pfargs, 'lutsize':self.lut.size,
                'kwargs':self._kwargs, 'type': self.__class__.__name__}

    def __str__(self):
        return 'LutEx(pftype="{:s}", pfargs={}, lutsize={:d}, '\
               '**kwargs={})'.format(
                   self._pf_obj.__class__.__name__, self._pfargs, self.lut.size,
                   self._kwargs)
