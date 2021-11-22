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

from xopto.mcml.mcdetector.base import Detector
from xopto.mcml import cltypes, mctypes, mcobject, mcoptions
from xopto.mcml.mcutil import axis
from xopto.mcml.mcutil.lut import CollectionLut
from xopto.mcml.mcutil import axis


class TotalPl(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClTotalPl(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            direction: mc_point3f_t
                Reference direction/orientation of the source.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle relative to teh direction.
            pl_min: mc_fp_t
                The leftmost edge of the first optical path length accumulator.
            inv_dpl: mc_fp_t
                Inverse of the width of the optical path length accumulators.
            n_pl: mc_size_t
                The number of path length accumulators.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            pl_log_scale: mc_int_t
                A flag indicating logarithmic scale of the path length axis.
            '''
            _fields_ = [
                ('direction', T.mc_point3f_t),
                ('cos_min', T.mc_fp_t),
                ('pl_min', T.mc_fp_t),
                ('inv_dpl', T.mc_fp_t),
                ('n_pl', T.mc_size_t),
                ('offset', T.mc_size_t),
                ('pl_log_scale', T.mc_int_t),
            ]
        return ClTotalPl

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_point3f_t direction;'
            '	mc_fp_t cos_min;',
            '	mc_fp_t pl_min;',
            '	mc_fp_t inv_dpl;',
            '	mc_size_t n_pl;',
            '	mc_size_t offset;',
            '	mc_int_t pl_log_scale;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the detector accumulator in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_detector(__mc_detector_mem const Mc{}Detector *detector){{'.format(loc, Loc),
            '	dbg_print("Mc{}Detector - TotalPl detector:");'.format(Loc),
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_float(INDENT "pl_min (um)", detector->pl_min*1e6f);',
            '	dbg_print_float(INDENT "inv_dpl (1/um)", detector->inv_dpl*1e-6f);',
            '	dbg_print_size_t(INDENT "n_pl:", detector->n_pl);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '	dbg_print_int(INDENT "pl_log_scale:", detector->pl_log_scale);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} TotalPl detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const struct Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	mc_fp_t pl = mcsim_optical_pathlength(mcsim);',
            '	if (detector->pl_log_scale)',
            '	    pl = mc_log(mc_fmax(pl, FP_PLMIN));',
            '	mc_int_t pl_index = mc_int((pl - detector->pl_min)*detector->inv_dpl);',
            '	pl_index = mc_clip(pl_index, 0, detector->n_pl - 1);',
            '',
            '	address = mcsim_accumulator_buffer_ex(mcsim, detector->offset + pl_index);',
            '',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(dot3f(dir, &detector->direction)));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} TotalPl detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def cl_options(self, mc, target=None) -> mcoptions.RawOptions:
        '''
        OpenCL kernel options defined by this object.
        '''
        return [('MC_TRACK_OPTICAL_PATHLENGTH', True)]

    def __init__(self, plaxis: axis.Axis or 'TotalPl' = None, 
                 cosmin: float = 0.0,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Total reflectance-transmittance and optical pathe length detector.

        Parameters
        ----------
        plaxis: axis.Axis or TotalPl
            Object that defines the accumulators along the optical path length
            axis (this axis supports log-scale). If an instance of TotalPl
            a new copy of TotalPl is created.
        cosmin: float
            Cosine of the maximum acceptance angle (relative to the direction)
            of the detector.
        direction: (float, float, float)
            Reference direction/orientation of the detector.
        '''
        if isinstance(plaxis, TotalPl):
            totalpl = plaxis
            cosmin = totalpl.cosmin
            nphotons = totalpl.nphotons
            direction = totalpl.direction
            raw_data = np.copy(totalpl.raw)
            plaxis = totalpl.plaxis
        else:
            if plaxis is None:
                plaxis = axis.Axis(0.0, 1.0, 1)
            raw_data = np.zeros((plaxis.n, 1))
            nphotons = 0

        super().__init__(raw_data, nphotons)

        self._cosmin = 0.0
        self._direction = np.zeros((3,))
        self._pl_axis = plaxis

        self._set_cosmin(cosmin)
        self._set_direction(direction)

    def _get_cosmin(self) -> Tuple[float, float]:
        return self._cosmin
    def _set_cosmin(self, value: float or Tuple[float, float]):
        self._cosmin = min(max(float(value), 0.0), 1.0)
    cosmin = property(_get_cosmin, _set_cosmin, None,
                      'Cosine of the maximum acceptance angle.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('Direction vector norm/length must not be 0!')
        self._direction *= 1.0/norm
    direction = property(_get_direction, _set_direction, None,
                         'Detector reference direction.')

    def _get_normalized(self) -> np.ndarray:
        return self.raw*(1.0/max(self.nphotons, 1.0))
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def _get_plaxis(self) -> axis.Axis:
        return self._pl_axis
    plaxis = property(_get_plaxis, None, None, 'Path length axis object.')

    def _get_pl(self):
        return self._pl_axis.centers
    pl = property(_get_pl, None, None,
                  'Centers of the optical pathlength axis accumulators.')

    def _get_pledges(self):
        return self._pl_axis.edges
    pledges = property(_get_pledges, None, None,
                       'Edges of the optical pathlength axis accumulators.')

    def _get_npl(self):
        return self._pl_axis.n
    npl = property(_get_npl, None, None,
                   'Number of accumulators in the optical pathlength axis.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`TotalPl.cl_type` method for a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        allocation = mc.cl_allocate_rw_accumulator_buffer(self, self.shape)
        target.offset = allocation.offset

        target.cos_min = self._cosmin
        target.direction.fromarray(self._direction)

        target.pl_min = self._pl_axis.scaled_start
        if self._pl_axis.step != 0.0:
            target.inv_dpl = 1.0/self._pl_axis.step
        else:
            target.inv_dpl = 0.0
        target.pl_log_scale = self._pl_axis.logscale
        target.n_pl = self._pl_axis.n

        return target

    def todict(self) -> dict:
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`TotalPl.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'TotalPl',
            'cosmin':self._cosmin,
            'direction': self._direction.tolist(),
            'pl_axis': self.plaxis.todict(),
        }

    @staticmethod
    def fromdict(data: dict) -> 'TotalPl':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`TotalPl.todict` method.
        '''
        detector_type = data.pop('type')
        if detector_type != 'TotalPl':
            raise TypeError(
                'Expected "TotalPl" type bot got "{}"!'.format(detector_type))

        pl_axis_data = data.pop('pl_axis')
        pl_axis_type = pl_axis_data.pop('type')
        plaxis = getattr(axis, pl_axis_type)(**pl_axis_data)

        return TotalPl(plaxis=plaxis, **data)

    def __str__(self):
        return 'TotalPl(plaxis={}, cosmin={}, direction=({}, {}, {}))'.format(
            self._pl_axis, self._cosmin, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))


class TotalLutPl(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClTotalLutPl(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            lut: mc_fp_lut_t
                Detector angular sensitivity lookup table. The lookup table
                is sampled with the absolute value of the incidence angle
                cosine compute relative to the detector reference direction.
            direction: mc_point3f_t
                Reference direction/orientation of the source.
            pl_min: mc_fp_t
                The leftmost edge of the first optical path length accumulator.
            inv_dpl: mc_fp_t
                Inverse of the width of the optical path length accumulators.
            n_pl: mc_size_t
                The number of path length accumulators.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            pl_log_scale: mc_int_t
                A flag indicating logarithmic scale of the path length axis.
            '''
            _fields_ = [
                ('lut', CollectionLut.cl_type(mc)),
                ('direction', T.mc_point3f_t),
                ('pl_min', T.mc_fp_t),
                ('inv_dpl', T.mc_fp_t),
                ('n_pl', T.mc_size_t),
                ('offset', T.mc_size_t),
                ('pl_log_scale', T.mc_int_t),
            ]
        return ClTotalLutPl

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_fp_lut_t lut;',
            '	mc_point3f_t direction;',
            '	mc_fp_t pl_min;',
            '	mc_fp_t inv_dpl;',
            '	mc_size_t n_pl;',
            '	mc_size_t offset;',
            '	mc_int_t pl_log_scale;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the detector accumulator in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_detector(__mc_detector_mem const Mc{}Detector *detector){{'.format(loc, Loc),
            '	dbg_print("Mc{}Detector - TotalLutPl detector:");'.format(Loc),
            '	dbg_print_fp_lut(INDENT "lut:", &detector->lut);',
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_float(INDENT "pl_min (um)", detector->pl_min*1e6f);',
            '	dbg_print_float(INDENT "inv_dpl (1/um)", detector->inv_dpl*1e-6f);',
            '	dbg_print_size_t(INDENT "n_pl:", detector->n_pl);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '	dbg_print_int(INDENT "pl_log_scale:", detector->pl_log_scale);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} TotalLutPl detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const struct Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	mc_fp_t pl = mcsim_optical_pathlength(mcsim);',
            '	if (detector->pl_log_scale)',
            '	    pl = mc_log(mc_fmax(pl, FP_PLMIN));',
            '	mc_int_t pl_index = mc_int((pl - detector->pl_min)*detector->inv_dpl);',
            '	pl_index = mc_clip(pl_index, 0, detector->n_pl - 1);',
            '',
            '	address = mcsim_accumulator_buffer_ex(mcsim, detector->offset + pl_index);',
            '',
            '	mc_fp_t sensitivity = FP_0;',
            '	fp_linear_lut_rel_sample(mcsim_fp_lut_array(mcsim),',
            '		&detector->lut, mc_fabs(dot3f(dir, &detector->direction)),',
            '		&sensitivity)',
            '	dbg_print_float("{} TotalLutPl sensitivity:", sensitivity);'.format(Loc),
            '',
            '	uint32_t ui32w = weight_to_int(weight*sensitivity);'
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} TotalLutPl detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def cl_options(self, mc, target=None) -> mcoptions.RawOptions:
        '''
        OpenCL kernel options defined by this object.
        '''
        return [('MC_TRACK_OPTICAL_PATHLENGTH', True)]

    def __init__(self, lut: CollectionLut, plaxis: axis.Axis = None,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Total reflectance-transmittance and optical path length detector with
        lookup table-based collection sensitivity.

        Parameters
        ----------
        lut: CollectionLut
            Lookup table of the detector angular sensitivity. The lookup table
                is sampled with the absolute value of the incidence angle
                cosine compute relative to the detector reference direction.
        plaxis: axis.Axis
            Object that defines the accumulators along the optical path length
            axis (this axis supports log-scale).
        direction: (float, float, float)
            Reference direction/orientation of the detector.

        Note
        ----
        The detector sensitivity must be valid for packets inside the detector,
        since the photon packets are handeled by the detector after entering
        (refracting into) the detector itself. 
        '''
        if isinstance(lut, TotalLutPl):
            totallutpl = lut
            lut = totallutpl.lut
            nphotons = totallutpl.nphotons
            direction = totallutpl.direction
            raw_data = np.copy(totallutpl.raw)
            plaxis = totallutpl.plaxis
        else:
            if plaxis is None:
                plaxis = axis.Axis(0.0, 1.0, 1)
            raw_data = np.zeros((plaxis.n, 1))
            nphotons = 0

        super().__init__(raw_data, nphotons)

        self._pl_axis = plaxis
        self._lut = None
        self._set_lut(lut)
        self._direction = np.zeros((3,))
        self._set_direction(direction)

    def _get_lut(self) -> CollectionLut:
        return self._lut
    def _set_lut(self, lut: CollectionLut):
        if not isinstance(lut, CollectionLut):
            raise TypeError(
                'Expected a CollectionLut but got {}!'.format(type(lut)))
        self._lut = lut
    lut = property(_get_lut, _set_lut, None,
                   'Lookup table of the detector angular sensitivity.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('Direction vector norm/length must not be 0!')
        self._direction *= 1.0/norm
    direction = property(_get_direction, _set_direction, None,
                        'Detector reference direction.')

    def _get_normalized(self) -> np.ndarray:
        return self.raw*(1.0/max(self.nphotons, 1.0))
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def _get_plaxis(self) -> axis.Axis:
        return self._pl_axis
    plaxis = property(_get_plaxis, None, None, 'Path length axis object.')

    def _get_pl(self):
        return self._pl_axis.centers
    pl = property(_get_pl, None, None,
                  'Centers of the optical pathlength axis accumulators.')

    def _get_pledges(self):
        return self._pl_axis.edges
    pledges = property(_get_pledges, None, None,
                       'Edges of the optical pathlength axis accumulators.')

    def _get_npl(self):
        return self._pl_axis.n
    npl = property(_get_npl, None, None,
                   'Number of accumulators in the optical pathlength axis.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`TotalLutPl.cl_type` method for a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        allocation = mc.cl_allocate_rw_accumulator_buffer(self, self.shape)
        target.offset = allocation.offset

        target.lut = self._lut.cl_pack(mc, target.lut)
        target.direction.fromarray(self._direction)

        target.pl_min = self._pl_axis.scaled_start
        if self._pl_axis.step != 0.0:
            target.inv_dpl = 1.0/self._pl_axis.step
        else:
            target.inv_dpl = 0.0
        target.pl_log_scale = self._pl_axis.logscale
        target.n_pl = self._pl_axis.n

        return target

    def todict(self) -> dict:
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`TotalLutPl.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'TotalLutPl',
            'lut':self._lut.todict(),
            'direction': self._direction.tolist(),
            'pl_axis': self.plaxis.todict(),
        }

    @staticmethod
    def fromdict(data: dict) -> 'TotalLutPl':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`TotalLutPl.todict` method.
        '''
        detector_type = data.pop('type')
        if detector_type != 'TotalLutPl':
            raise TypeError(
                'Expected "TotalLutPl" type bot got "{}"!'.format(detector_type))

        lut = CollectionLut(data.pop('lut'))

        pl_axis_data = data.pop('pl_axis')
        pl_axis_type = pl_axis_data.pop('type')
        plaxis = getattr(axis, pl_axis_type)(**pl_axis_data)

        return TotalLutPl(lut, plaxis=plaxis, **data)

    def __str__(self):
        return 'TotalLutPl(lut={}, plaxis={}, direction=({}, {}, {}))'.format(
            self._lut, self._pl_axis, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
