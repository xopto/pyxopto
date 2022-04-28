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

from xopto.mcbase.mcutil.lut import EmissionLut, CollectionLut
from xopto.mcbase import cltypes
from xopto.mcbase import mcobject
from xopto.mcbase import mctypes


class MultimodeFiber:
    '''
    Class representing Multimode step-index fibers with an ideal
    emission and collection characteristics defined by the numerical
    aperture. 
    '''

    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClMultimodeFiber(cltypes.Structure):
            _field_ = [
                ('ncore', T.mc_fp_t),
                ('ncladding', T.mc_fp_t),
                ('rcore', T.mc_fp_t),
                ('rcladding', T.mc_fp_t),
                ('na', T.mc_fp_t),
            ]
        return ClMultimodeFiber

    @staticmethod
    def compute_na(ncore: float, ncladding: float) -> float:
        '''
        Computes the numerical aperture of the fiber core given the refractive
        indices of the fiber core and fiber cladding.

        Parameters
        ----------
        ncore: float
            Refractive index of the fiber core.
        ncladding: float
            Refractive index of the fiber cladding.

        Returns
        -------
        na: float
            Numerical aperture of the fiber core.
        '''
        return (ncore**2 - ncladding**2)**0.5

    @staticmethod
    def compute_ncladding(ncore: float, na: float) -> float:
        '''
        Computes the refractive index of the fiber cladding given the
        refractive index and numerical aperture of the fiber core.

        Parameters
        ----------
        ncore: float
            Refractive index of the fiber core.
        na: float
            Numerical aperture of the fiber core.

        Returns
        -------
        ncladding: float
            Refractive index of the fiber cladding.
        '''
        return (ncore**2 - na**2)**0.5

    def __init__(self, dcore: float, dcladding: float, ncore: float, na: float):
        '''
        Multimode fiber constructor.

        Parameters
        ----------
        dcore: float
            Diameter (m) of the fiber core.
        dcladding: float
            Diameter (m) of the fiber cladding.
        ncore: float
            Refractive index of the fiber core.
        na: float
            Numerical aperture of the fiber core.
        '''
        self._dcore = float(dcore)
        self._dcladding = float(dcladding)
        self._ncore = float(ncore)
        self._na = float(na)

    def validate(self):
        '''
        Check the integrity of the object parameters. Raises ValueError.
        '''
        if self._dcore > self._dcladding:
            raise ValueError(
            'Fiber core diameter must be smaller than the fiber '
            'cladding diameter!')

    def _get_dcore(self) -> float:
        return self._dcore
    def _set_dcore(self, value: float):
        self._dcore = float(value)
    dcore = property(_get_dcore, _set_dcore, None,
                     'Diameter (m) of the fiber core.')

    def _get_dcladding(self) -> float:
        return self._dcladding
    def _set_dcladding(self, value: float):
        self._dcladding = float(value)
    dcladding = property(_get_dcladding, _set_dcladding, None,
                         'Diameter (m) of the fiber cladding.')

    def _get_ncore(self) -> float:
        return self._ncore
    def _set_ncore(self, value: float):
        self._ncore = float(value)
    ncore = property(_get_ncore, _set_ncore, None,
                     'Refractive index of the fiber core.')

    def _get_ncladding(self) -> float:
        return MultimodeFiber.compute_ncladding(self._ncore, self._na)
    ncladding = property(_get_ncladding, None, None,
                         'Refractive index of the fiber cladding computed as '
                         '$\\sqrt(ncore^2 - na^2)$')

    def _get_na(self) -> float:
        return self._na
    def _set_na(self, value: float):
        self._na = float(value)
    na = property(_get_na, _set_na, None,
                  'Numerical aperture (NA) of the fiber core.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the OpenCL Structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`~MultimodeFiber.cl_type`
        method for a detailed list of fields.

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
            instance of ClMultimodeFiber if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type()
            target = target_type()

        target.ncore = self.ncore
        target.ncladding = self.ncladding
        target.rcore = self.dcore*0.5
        target.dcladding = self.dcladding*0.5
        target.na = self.na

        return target

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'dcore': self._dcore, 
                'dcladding': self._dcladding,
                'ncore': self._ncore,
                'na': self._na,
                'type': self.__class__.__name__}

    @staticmethod
    def fromdict(data: dict) -> 'MultimodeFiber':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the py:meth:`MultimodeFiber.todict` method.
        '''
        obj_type = data.pop('type')
        if obj_type != 'MultimodeFiber':
            raise TypeError('Expected "MultimodeFiber" type but '
                            'got "{}"!'.format(obj_type))

        return MultimodeFiber(**data)

    def __repr__(self):
        return 'MultimodeFiber(dcore={:f}, dcladding={:f}, ' \
               'ncore={:f}, na={:f})'.format(
                   self._dcore, self._dcladding, self._ncore, self._na)

    def __str__(self):
        return self.__repr__()


class MultimodeFiberLut:
    '''
    Class representing multimode fibers with arbitrary emission and
    collection characteristics defined by a linear lookup table.
    '''

    @staticmethod
    def compute_na(ncore: float, ncladding: float) -> float:
        '''
        Computes the numerical aperture of the fiber core given the refractive
        indices of the fiber core and fiber cladding.

        Parameters
        ----------
        ncore: float
            Refractive index of the fiber core.
        ncladding: float
            Refractive index of the fiber cladding.

        Returns
        -------
        na: float
            Numerical aperture of the fiber core.
        '''
        return (ncore**2 - ncladding**2)**0.5

    @staticmethod
    def compute_ncladding(ncore: float, na: float) -> float:
        '''
        Computes the refractive index of the fiber cladding given the
        refractive index and numerical aperture of the fiber core.

        Parameters
        ----------
        ncore: float
            Refractive index of the fiber core.
        na: float
            Numerical aperture of the fiber core.

        Returns
        -------
        ncladding: float
            Refractive index of the fiber cladding.
        '''
        return (ncore**2 - na**2)**0.5

    def __init__(self, dcore: float, dcladding: float,
                 ncore: float, ncladding: float,
                 emission: EmissionLut=None, collection: CollectionLut=None):
        '''
        Multimode fiber constructor.

        Parameters
        ----------
        dcore: float
            Diameter (m) of the fiber core.
        dcladding: float
            Diameter (m) of the fiber cladding.
        ncore: float
            Refractive index of the fiber core.
        ncladding: float
            Refractive index of the fiber cladding.
        emission: EmissionLut or np.ndarray or str
            Lookup table for numerical sampling of the emission angles. Can
            be an instance of EmissionLut or a file from which to load the
            EmissionLut.
        collection: CollectionLut or str
            Collection sensitivity/efficiency of the fiber defined by a linear
            lookup table array.
            Can be an instance of CollectionLut or a file from which to load
            the CollectionLut.
        '''
        self._dcore = float(dcore)
        self._dcladding = float(dcladding)
        self._ncore = float(ncore)
        if ncladding is None:
            self._ncladding = self._ncore

        if isinstance(emission, str):
            emission = EmissionLut.fromfile(emission)
        elif isinstance(emission, np.ndarray):
            emission = EmissionLut(emission)

        if isinstance(collection, str):
            collection = CollectionLut.fromfile(collection)

        self._emission_lut = emission
        self._collection_lut = collection

    def validate(self):
        '''
        Check the integrity of parameters. Raises a ValueError.
        '''
        if self._dcore > self._dcladding:
            raise ValueError(
            'Fiber core diameter must be smaller than the fiber '
            'cladding diameter!')

    def _get_dcore(self) -> float:
        return self._dcore
    def _set_dcore(self, value:float):
        self._dcore = float(value)
    dcore = property(_get_dcore, _set_dcore, None,
                     'Diameter (m) of the fiber core.')

    def _get_dcladding(self) -> float:
        return self._dcladding
    def _set_dcladding(self, value: float):
        self._dcladding = float(value)
    dcladding = property(_get_dcladding, _set_dcladding, None,
                         'Diameter (m) of the fiber cladding.')

    def _get_ncore(self) -> float:
        return self._ncore
    def _set_ncore(self, value: float):
        self._ncore = float(value)
    ncore = property(_get_ncore, _set_ncore, None,
                     'Refractive index of the fiber core.')

    def _get_ncladding(self) -> float:
        return self._ncladding
    ncladding = property(_get_ncladding, None, None,
                         'Refractive index of the fiber cladding computed as '
                         '$\\sqrt(ncore^2 - na^2)$')

    def _get_emission_lut(self) -> EmissionLut:
        return self._emission_lut
    def _set_emission_lut(self, value: EmissionLut):
        self._nemission_lut = EmissionLut(value)
    emission = property(_get_emission_lut, _set_emission_lut, None,
                        'Emission lookup table.')

    def _get_collection_lut(self) -> CollectionLut:
        return self._collection_lut
    def _set_collection_lut(self, value: CollectionLut):
        self._collection_lut = CollectionLut(value)
    collection = property(_get_collection_lut, _set_collection_lut, None,
                         'Collection lookup table.')

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'dcore': self._dcore, 
                'dcladding': self._dcladding,
                'ncore': self._ncore,
                'emission': self._emission_lut,
                'collection': self._collection_lut,
                'type': self.__class__.__name__}

    def __repr__(self):
        return 'MultimodeFiberLut(dcore={:f}, dcladding={:f}, ' \
               'ncore={:f}, emission={}, collection={})'.format(
                   self._dcore, self._dcladding, self._ncore,
                   self._emission_lut, self._collection_lut)

    def __str__(self):
        return self.__repr__()


class FiberLayout:
    def __init__(self,
                 fiber: MultimodeFiber or MultimodeFiberLut or 'FiberLayout',
                 position: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 direction: Tuple[float, float, float] = (0.0, 0.0 , 1.0)):
        '''
        Optical fiber layout constructor.

        Parameters
        ----------
        fiber: MultimodeFiber or MultimodeFiberLut or FiberLayout
            A multimode optical fiber instance or a fiber layout instance.
            If a fiber layout instance, a new copy is created.
        position: tuple(float, float, float)
            Position of the fiber tip center.
        direction:
            Direction/orientation of the fiber. 
        --------
        '''
        if isinstance(fiber, FiberLayout):
            af = fiber
            fiber = af.fiber
            position = af.position
            direction = af.direction

        self._fiber = fiber
        self._position = np.zeros((3,))
        self._direction = np.array((0.0, 0.0, 1.0))

        self._set_position(position)
        self._set_direction(direction)

    def _get_fiber(self) -> MultimodeFiber or MultimodeFiberLut:
        return self._fiber
    def _set_fiber(self, value: MultimodeFiber or MultimodeFiberLut):
        self._fiber = value
    fiber = property(_get_fiber, _set_fiber, None,
                     'Properties of the optical fiber.')

    def _get_position(self) -> Tuple[float, float, float]:
        return self._position
    def _set_position(self, value: float or Tuple[float, float] or Tuple[float, float, float]):
        if hasattr(value, '__len__'):
            self._position[:len(value)] = value
        else:
            self._position[:] = value
    position = property(_get_position, _set_position, None,
                       'Position of the fiber as a tuple (x, y, z).')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError(
                'Fiber direction vector norm/length must not be 0!')
        self._direction *= 1.0/norm
    direction = property(_get_direction, _set_direction, None,
                         'Fiber reference direction.')

    def todict(self) -> dict:
        '''
        Export object to a dict.
        '''
        return {'fiber': self._fiber.todict(), 
                'position': self._position.tolist(),
                'direction': self._direction.tolist(),
                'type': self.__class__.__name__}

    @staticmethod
    def fromdict(data: dict) -> 'MultimodeFiber':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the py:meth:`FiberLayout.todict` method.
        '''
        obj_type = data.pop('type')
        if obj_type != 'FiberLayout':
            raise TypeError('Expected "FiberLayout" type but '
                            'got "{}"!'.format(obj_type))

        return FiberLayout(fiber=MultimodeFiber(data.pop('fiber')), **data)

    def __str__(self):
        return 'FiberLayout(fiber={}, position=({}, {}, {}), ' \
               'direction=({}, {}, {}))'.format(
                   self._fiber, *self._position, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
