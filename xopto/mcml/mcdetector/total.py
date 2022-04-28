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
from xopto.mcml import cltypes, mctypes, mcobject
from xopto.mcml.mcutil import axis
from xopto.mcml.mcutil.lut import CollectionLut


class Total(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClTotal(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            direction: mc_point3f_t
                Reference direction/orientation of the source.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle relative to teh direction.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            '''
            _fields_ = [
                ('direction', T.mc_point3f_t),
                ('cos_min', T.mc_fp_t),
                ('offset', T.mc_size_t),
            ]
        return ClTotal

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
            '	mc_size_t offset;',
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
            '	dbg_print("Mc{}Detector - Total detector:");'.format(Loc),
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} Total detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const struct Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	address = mcsim_accumulator_buffer_ex(mcsim, detector->offset);',
            '',
            '	mc_point3f_t detector_direction = detector->direction;',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(mc_dot_point3f(dir, &detector_direction)));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} Total detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def __init__(self, cosmin: float = 0.0,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Total reflectance-transmittance detector.

        Parameters
        ----------
        cosmin: float
            Cosine of the maximum acceptance angle (relative to the direction)
            of the detector.
        direction: (float, float, float)
            Reference direction/orientation of the detector.
        '''
        if isinstance(cosmin, Total):
            total = cosmin
            cosmin = total.cosmin
            nphotons = total.nphotons
            direction = total.direction
            raw_data = np.copy(total.raw)
        else:
            raw_data = np.zeros((1,))
            nphotons = 0

        super().__init__(raw_data, nphotons)

        self._cosmin = 0.0
        self._direction = np.zeros((3,))
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

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`Total.cl_type` method for a detailed list of fields.

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

        return target

    def todict(self) -> dict:
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`Total.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'Total',
            'cosmin':self._cosmin,
            'direction': self._direction.tolist()
        }

    @staticmethod
    def fromdict(data: dict):
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`Total.todict` method.
        '''
        data_ = dict(data)
        detector_type = data_.pop('type')
        if detector_type != 'Total':
            raise TypeError(
                'Expected "Total" type bot got "{}"!'.format(detector_type))

        return Total(**data_)

    def __str__(self):
        return 'Total(cosmin={}, direction=({}, {}, {}))'.format(
            self._cosmin, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))


class TotalLut(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClTotal(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            direction: mc_point3f_t
                Reference direction/orientation of the source.
            lut: mc_fp_lut_t
                Detector angular sensitivity lookup table. The lookup table
                is sampled with the absolute value of the incidence angle
                cosine computed relative to the detector reference direction.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            '''
            _fields_ = [
                ('lut', CollectionLut.cl_type(mc)),
                ('direction', T.mc_point3f_t),
                ('offset', T.mc_size_t),
            ]
        return ClTotal

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_fp_lut_t lut;',
            '	mc_point3f_t direction;'
            '	mc_size_t offset;',
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
            '	dbg_print("Mc{}Detector - TotalLut detector:");'.format(Loc),
            '	dbg_print_fp_lut(INDENT "lut:", &detector->lut);',
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} TotalLut detector hit");'.format(Loc),
            '',
            '	__mc_detector_mem const struct Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	address = mcsim_accumulator_buffer_ex(mcsim, detector->offset);',
            '',
            '	mc_fp_t sensitivity = FP_0;',
            '	mc_point3f_t detector_direction = detector->direction;',
            '	fp_linear_lut_rel_sample(mcsim_fp_lut_array(mcsim),',
            '		&detector->lut, mc_fabs(mc_dot_point3f(dir, &detector_direction)),',
            '		&sensitivity);',
            '	dbg_print_float("{} TotalLut sensitivity:", sensitivity);'.format(Loc),
            '',
            '	uint32_t ui32w = weight_to_int(weight*sensitivity);'
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} TotalLut detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def __init__(self, lut: CollectionLut,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        Total reflectance-transmittance detector.

        Parameters
        ----------
        lut: CollectionLut
            Lookup table of the detector angular sensitivity. The lookup table
                is sampled with the absolute value of the incidence angle
                cosine compute relative to the detector reference direction.
        direction: (float, float, float)
            Reference direction/orientation of the detector.

        Note
        ----
        The detector sensitivity must be valid for packets inside the detector,
        since the photon packets are handeled by the detector after entering
        (refracting into) the detector itself. 
        '''
        if isinstance(lut, TotalLut):
            total = lut
            lut = total.lut
            nphotons = total.nphotons
            direction = total.direction
            raw_data = np.copy(total.raw)
        else:
            raw_data = np.zeros((1,))
            nphotons = 0

        super().__init__(raw_data, nphotons)

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

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`TotalLut.cl_type` method for a detailed list of fields.

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

        return target

    def todict(self) -> dict:
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`Total.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'TotalLut',
            'lut':self._lut.todict(),
            'direction': self._direction.tolist()
        }

    @staticmethod
    def fromdict(data: dict):
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`Total.todict` method.
        '''
        data_ = dict(data)
        detector_type = data_.pop('type')
        if detector_type != 'TotalLut':
            raise TypeError(
                'Expected "TotalLut" type bot got "{}"!'.format(detector_type))
        lut = CollectionLut(data_.pop('lut'))
        return TotalLut(lut, **data_)

    def __str__(self):
        return 'TotalLut(lut={}, direction=({}, {}, {}))'.format(
            self._lut, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
