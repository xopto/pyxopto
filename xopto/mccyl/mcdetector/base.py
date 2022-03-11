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

from typing import List, Dict

import numpy as np

from xopto.mcml import cltypes
from xopto.mcml import mcobject
from xopto.mcml import mctypes
from xopto.mcml import mcoptions


NONE = 'none'
OUTER = 'outer'
SPECULAR = 'specular'


class DetectorBase(mcobject.McObject):
    '''
    Base class of all detectors.
    '''
    def __init__(self, location: NONE or OUTER or SPECULAR):
        '''
        Detector at the sample surface.

        Parameters
        ----------
        location: OUTER or SPECULAR
            Location of the detector.
        '''
        super().__init__()
        if location not in (NONE, OUTER, SPECULAR):
            raise ValueError(
                'Detector location must be "{}", "{}" or "{}"!'.format(
                    NONE, OUTER, SPECULAR)
            )
        self._location = location

    def _get_location(self) -> str:
        return self._location
    def _set_location(self, location: OUTER or SPECULAR):
        if location not in (OUTER, SPECULAR):
            raise ValueError(
                'Surface detector location must be "{}" or "{}"!'.format(
                    OUTER, SPECULAR)
            )
        if location != self._location and self._location != NONE:
            raise RuntimeError('Detector location cannot be changed!')
        self._location = str(location)
    location = property(_get_location, _set_location, None,
                        'Location of the detector.')


class DetectorOuter(DetectorBase):
    '''
    Base class of detectors at the outer sample surface.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create a detector at the outer sample surface.
        '''
        super().__init__(OUTER)


class DetectorSpecular(DetectorBase):
    '''
    Base class of the detectors for specular reflections.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create a detector for specular reflections.
        '''
        super().__init__(SPECULAR)


class DetectorAny(DetectorBase):
    '''
    Base class of the detectors at the outer sample surface or
    for specular reflections.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create a detector for any location.
        '''
        super().__init__(NONE)

class DetectorDefault(DetectorAny):
    '''
    Default / dummy detector for the outer sample surface or for
    specular reflections.
    '''
    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Structure that is passed to the Monte carlo simulator kernel.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents the detectors in
            the Monte Carlo kernel.

            The returned structure type implements the following fields:

            - dummy: mc_int_t
                Dummy field of the dummy detector.
        '''
        T = mc.types
        class ClDetectorDefault(cltypes.Structure):
            _fields_ = [('dummy', T.mc_int_t)]
        return ClDetectorDefault

    @staticmethod
    def cl_options(mc: mcobject.McObject) -> str:
        return ''

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the default detector in the Monte Carlo simulator.
        '''
        Loc = self.location.capitalize()
        return 'struct Mc{}Detector{{mc_int_t dummy;}};'.format(Loc)

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL implementation of the default detector
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_detector(__mc_detector_mem const Mc{}Detector *detector){{'.format(loc, Loc),
            '	dbg_print("Mc{}Detector - default detector:");'.format(Loc),
            '};',
        ))

    def cl_pack(self, mc: mcobject.McObject, \
                target : cltypes.Structure = None) -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`DetectorDefault.cl_type` method for a detailed list
        of fields.

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

        target.dummy = 0

        return target

    def todict(self) -> dict:
        '''
        Export object to dict.
        '''
        return {'type': self.__class__.__name__}

    @classmethod
    def fromdict(cls, data: dict) -> 'DetectorDefault':
        '''
        Create an instance of :py:class:`DetectorDefault` from a
        dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`DetectorDefault.todict` method.
        '''
        detector_type = data.pop('type')
        if detector_type != cls.__name__:
            raise TypeError(
                'Expected a "{}" type bot got "{}"!'.format(
                    cls.__name__, detector_type))

        return cls(**data)


class Detector(DetectorAny):
    def __init__(self, raw_data: np.ndarray, nphotons: int):
        '''
        Detector base class.

        Parameters
        ----------
        raw_data: np.ndarray
            Raw data accumulator.
        nphotons: int
            Number of photon packets used in the simulations.
        '''
        super().__init__()
        self._nphotons = int(nphotons)
        self._raw_data = raw_data

    def _get_raw(self) -> np.ndarray:
        return self._raw_data
    raw = property(_get_raw, None, None,
                   'Raw accumulator data array.')

    def _get_nphotons(self) -> int:
        return self._nphotons
    nphotons = property(_get_nphotons, None, None,
                        'The number of photon packets that produced '\
                        'the raw data accumulator content.')

    def _get_shape(self) -> tuple:
        return self._raw_data.shape
    shape = property(_get_shape, None, None,
                     'Shape of the raw data accumulator array.')

    def _get_raw_total(self):
        return self._raw_data.sum()
    total = property(_get_raw_total, None, None,
                     'Total weight of the accumulated photon packets.')

    def set_raw_data(self, data: np.ndarray, nphotons: int):
        self._raw_data[:] = np.asarray(data, dtype=self._raw_data.dtype)
        self._nphotons = int(nphotons)

    def update_data(self, mc: mcobject.McObject,
                    accumulators: List[np.ndarray], nphotons: int, **kwargs):
        '''
        Update the content of raw data accumulators with simulation results.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that produced the data.
        accumulators: List[np.ndarray]
            List of allocated accumulators (this implementation uses only one).
        nphotons: int
            The number of photon packets that produced the raw data
            accumulator content.
        kwargs: dict
            Additional keyword arguments not used by this implementation.
        '''
        new_data = np.reshape(accumulators[0], self.shape)
        if self._raw_data is not None:
            self._raw_data += new_data*(1.0/mc.types.mc_accu_k)
            self._nphotons += int(nphotons)
        else:
            self._raw_data[:] = new_data*(1.0/mc.types.mc_accu_k)
            self._nphotons = int(nphotons)


class Detectors(mcobject.McObject):
    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Structure that is passed to the Monte carlo simulator kernel.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents surface detectors in the
            Monte Carlo kernel.

            The returned structure type implements the following fields:

            - outer: cltypes.Structure
                Detector at the outer sample surface.
            - specular: cltypes.Structure
                Detector of specular reflections.
        '''
        T = mc.types
        class ClDetectors(cltypes.Structure):
            _fields_ = [
                ('outer', self._outer.fetch_cl_type(mc)),
                ('specular', self._specular.fetch_cl_type(mc))
            ]
        return ClDetectors

    def cl_options(self, mc: mcobject.McObject) -> str:
        '''
        Returns the OpenCL options of the detectors.

        If the detector for the outer sample surface or the detector
        of specular reflection are specified, this function activates
        the OpenCL options that enable the corresponding outer surface
        detector, i.e. MC_USE_OUTER_DETECTOR for the detector at the outer
        sample surface and MC_USE_SPECULAR_DETECTOR for detector of specular
        reflections.
        '''
        options = []
        use_detectors = False

        if type(self._outer) != DetectorDefault:
            options.append(('MC_USE_OUTER_DETECTOR', True))
            options.extend(self._outer.fetch_cl_options(mc))
            use_detectors = True
        if type(self._specular) != DetectorDefault:
            options.append(('MC_USE_SPECULAR_DETECTOR', True))
            options.extend(self._specular.fetch_cl_options(mc))
            use_detectors = True
        if use_detectors:
            options.append(('MC_USE_DETECTORS', True))

        return options

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Declarations of detectors in OpenCL.
        '''
        return '\n'.join((
            self._outer.fetch_cl_declaration(mc),
            'typedef struct McOuterDetector McOuterDetector;',
            '',
            self._specular.fetch_cl_declaration(mc),
            'typedef struct McSpecularDetector McSpecularDetector;',
            '',
            'struct MC_STRUCT_ATTRIBUTES McDetectors{',
            '	McOuterDetector outer;',
            '	McSpecularDetector specular;',
            '};',
            '',
            'inline void mcsim_outer_detector_deposit(',
            '	McSim *mcsim, mc_point3f_t const *pos, mc_point3f_t const *dir,',
            '	mc_fp_t weight);',
            'inline void mcsim_specular_detector_deposit(',
            '	McSim *mcsim, mc_point3f_t const *pos, mc_point3f_t const *dir,',
            '	mc_fp_t weight);',
            '',
            '#define mcsim_outer_detector(psim) (&mcsim_detectors(psim)->outer)',
            '#define mcsim_specular_detector(psim) (&mcsim_detectors(psim)->specular)',
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the detectors.
        '''
        return '\n'.join((
            self._outer.fetch_cl_implementation(mc),
            '',
            self._specular.fetch_cl_implementation(mc),
            '',
            'void dbg_print_detectors(__mc_detector_mem const McDetectors *detectors){',
            '	dbg_print("McDetectors");',
            '	dbg_print_outer_detector(&detectors->outer);',
            '	dbg_print_specular_detector(&detectors->specular);',
            '};',
        ))

    def __init__(self, outer: Detector or 'Detectors' = None,
                 specular: Detector = None):
        '''
        Create an instance of simulator detectors.

        Parameters
        ----------
        outer: DetectorBase or Detectors
            Detector at the outer sample surface, None, or an existing Detectors
            instance. If an existing Detectors instance, a new copy will
            be made.
        specular: DetectorBase
            Detector of specular reflections or None.
        '''
        super().__init__()
        if isinstance(outer, Detectors):
            detectors = outer

            outer = detectors.outer
            outer = type(outer)(outer)

            specular = detectors.specular
            specular = type(specular)(specular)
        else:
            if outer is None:
                outer = DetectorDefault()
            if specular is None:
                specular = DetectorDefault()

        outer.location = OUTER
        specular.location = SPECULAR

        self._outer = outer
        self._specular = specular

    def update_data(self, mc: mcobject.McObject, detector: Detector or str,
                    data: Dict[np.dtype, List[np.ndarray]],
                    nphotons: int = 0):
        '''
        Update the detector data with simulation results.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that produced the results.
        detector: Detector or str
            Detector or location to update.
        data: Dict[np.dtype: List[np.ndarray]]
            A dict of list of numpy data buffers downloaded from the kernel.
        nphotons: int
            The number of photon packets that produced the results.
        '''
        if isinstance(detector, Detector):
            location = detector.location
        else:
            location = detector

        if location not in (OUTER, SPECULAR):
            raise ValueError(
                'Detector location must be one of "{}" or "{}" '
                'but got "{}"!'.format(
                    OUTER, SPECULAR, location
                )
            )
        accumulators = data.get(np.dtype(mc.types.np_accu))
        float_buffers = data.get(np.dtype(mc.types.np_float))
        integer_buffers = data.get(np.dtype(mc.types.np_int))
        getattr(self, location).update_data(
            mc, accumulators=accumulators, float_buffers=float_buffers,
            integer_buffers=integer_buffers, nphotons=nphotons)

    def _get_outer(self) -> Detector:
        return self._outer
    outer = property(_get_outer, None, None,
                   'Detector at the outer sample surface.')

    def _get_specular(self) -> Detector:
        return self._specular
    specular = property(_get_specular, None, None,
                        'Detector of specular reflections.')

    def __iter__(self):
        return iter([self._outer, self._specular])


    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`Detectors.cl_type` method for a detailed list of
        fields.

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

        target.outer = self._outer.cl_pack(mc, target.outer)
        target.specular = self._specular.cl_pack(mc, target.specular)

        return target

    def types(self) -> tuple:
        '''
        Returns a tuple of detector types assigned to this instance.
        '''
        return type(self._outer), type(self._specular)

    def todict(self) -> dict:
        '''
        Export object to dict.
        '''

        return {
            'type': self.__class__.__name__,
            'outer': self._outer.todict(),
            'specular': self._specular.todict()
        }

    def __str__(self):
        return 'Detectors(outer={}, specular={})'.format(
            self._outer, self._specular)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
