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

from typing import List, Tuple

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from xopto.mcbase import mcobject
from xopto.mcbase import cltypes


class LinearLut(mcobject.McObject):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Create a Structure type that represents LinearLut in the OpenCL
        kernel.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance.

        Returns
        -------
        struct: Structure
            OpenCL structure type that is used to represents a LinearLut
            instance in the OpenCL kernel.

            The returned structure type implements the following fields:

            - first: mc_fp_t
                The value of the independent variable for the first element in
                the lookup table.
            - inv_span: mc_fp_t
                Inverse of the span of the independent variable (the span
                is defined as the difference between the value of the
                independent variable of the last and first element in the
                lookup table).
            - n: mc_size_t
                The number of elements in the lookup table.
            - offset: mc_size_t
                Offset/index of the first element in the lookup table
                from the start of the linear array that holds data of all
                the lookup tables.
        '''
        T = mc.types
        class ClLinearLut(cltypes.Structure):
            _fields_ = [
                ('first', T.mc_fp_t),
                ('inv_span', T.mc_fp_t),
                ('n', T.mc_size_t),
                ('offset', T.mc_size_t)
            ]

        return ClLinearLut

    @staticmethod
    def fromfile(filename: str) -> 'LinearLut':
        '''
        Load LinearLut instance from a npz file.

        Parameters
        ----------
        filename: str
            Data file.

        Returns
        -------
        lut: LinearLut
            A new LinearLut instance created from the file.
        '''
        data = np.load(filename)
        return LinearLut(data['lut_data'], data['first'], data['last'])

    def __init__(self, lut_data: np.ndarray or 'LinearLut',
                 first: float = 0.0, last: float = 1.0):
        '''
        Creates a linear lookup table instance. Samples in the lookup table
        are uniformly spread between the positions of the first and last
        lookup table element that are passed to the constructor.

        Parameters
        ----------
        lut_data: LinearLut, np.ndarray vector or str
            If a numpy array, lookup table data as an array or array-like object.
            If a str, a filename from which to load the instance.
            If a LinearLut, a weak copy is created.
        first: float
            Position of the first element/entry in the lookup table.
        last: float
            Position of the last element/entry in the lookup table.
        '''
        if isinstance(lut_data, LinearLut):
            lut = lut_data
            self._first = lut.first
            self._last = lut.last
            self._lut_data = np.copy(lut.data)
        elif isinstance(lut_data, str):
            np_data = np.load(lut_data)
            self._first = np_data['first']
            self._last = np_data['last']
            self._lut_data = np_data['lut_data']
        else:
            self._first = float(first)
            self._last = float(last)
            self._lut_data = np.asarray(lut_data, dtype=np.float64)

    def _get_first(self) -> float:
        return self._first
    first = property(
        _get_first, None, None,
        'Position of the first point in the lookup table.')

    def _get_last(self) -> float:
        return self._last
    last = property(
        _get_last, None, None,
        'Position of the last point in the lookup table.')

    def _get_span(self) -> Tuple[float, float]:
        return self._last - self._first
    span = property(
        _get_span, None, None,
        'Distance between the positions of the first and last point in the '
        'lookup table.')

    def _get_lut_data(self) -> np.ndarray:
        return self._lut_data
    data = property(
        _get_lut_data, None, None,
        'Lookup table data as a flat numpy array.')

    def __call__(self, x: np.ndarray) -> np.ndarray:
        n = self._lut_data.size
        fp_ind = (np.asarray(x) - self._first)/(self._last - self._first)*(n - 1)
        ind_1 = np.clip(np.floor(fp_ind), 0, n - 1).astype(np.intp)
        ind_2 = np.clip(ind_1 + 1, 0, n - 1)
        w = fp_ind - ind_1

        return (1.0 - w)*self._lut_data[ind_1] + w*self._lut_data[ind_2]

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator.

        Parameters
        ----------
        mc: McObject
            Simulator instance.
        target: cltypes.Struct
            Ctypes structure to be filled with the lut parameters.

        Returns
        -------
        lut: cltypes.Struct
            Filled ctypes structured received as an input arguments.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        lut_entry = mc.append_r_lut(self._lut_data)

        target.offset = lut_entry.offset
        target.first = self.first
        target.inv_span = 1.0/(self.last - self.first)
        target.n = self.data.size

        return target

    def todict(self, np2list: bool = True) -> dict:
        '''
        Save the lookup table configuration and data to a dictionary.

        Parameters
        ----------
        np2list: bool
            Convert numpy data array (lookup table) to a python list.

        Returns
        -------
        data: dict
            Linear lookup table configuration as a dictionary.
        '''
        return {
            'type': 'LinearLut',
            'first': self._first,
            'last': self._last,
            'lut_data': self._lut_data.tolist() if np2list else self._lut_data
        }

    @staticmethod
    def fromdict(data: dict) -> 'LinearLut':
        '''
        Create a linear lookup table instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`LinearLut.todict` method.
        '''
        rt_type = data.pop('type')
        if rt_type != 'LinearLut':
            raise TypeError('Expected "LinearLut" type bot got "{}"!'.\
                format(rt_type))

        return LinearLut(**data)

    def save(self, filename: str):
        '''
        Save instance to a compressed numpy file with extension .npz. 
        '''
        np.savez_compressed(filename, **self.todict(np2list=False))


class EmissionLut(LinearLut):
    def __init__(self, radiance: np.ndarray or 'EmissionLut' or str,
                 costheta: np.ndarray=None,
                 n: int=2000, npts: int=10000, meth: str='simps'):
        '''
        Constructs a lookup table for numerical sampling of the emission
        angles of a source with the given angular emission
        characteristics. Azimuthal symmetry of the emission characteristics
        is assumed.

        The constructed lookup table should be sampled by a uniform random
        variable :math:`r` from [0, 1] and linear interpolation that will
        yield the emission angle cosine (costheta):
        .. math:

            f_{i} = r*n
            i_1 = int(fp_index)
            i_2 = min(i_1 + 1, n - 1)
            \\delta = f_i - i_1
            costheta = lut[i_1]*(1.0 - \\delta) + lut[i_2 + 1]*\\delta

        Parameters
        ----------
        radiance: np.ndarray, EmissionLut or str
            Relative source radiance as a function of the angle of incidence
            cosines given in the costheta array. If the value is an instance
            of EmissionLut, a new weak copy of the instance is created. If the
            value is a str, an instance is loaded from the file represented
            by the string.
        costheta: np.ndarray
            Cosines of the angle of incidence at which the values in the
            radiance array are defined.
            If None, a uniformly distributed vector of values from 0 to 1
            is used.
        n: int
            Size of the lookup table that can be used to numerically sample
            the emission angles.
        npts: int
            Number of uniformly distributed control points from
            min(costheta) to max(costheta) used for simpson approximation
            of the radiance CDF.
        meth: str
            Method that is used to compute the radiance CDF. Use "simps"
            for Simpson or "quad" for adaptive-step integration.
        '''
        lut_random = np.linspace(0.0, 1.0, n)

        if isinstance(radiance, (str, EmissionLut)):
            super().__init__(radiance)

        else: 
            radiance = np.asarray(radiance, dtype=np.float64)

            if costheta is None:
                costheta = np.linspace(0.0, 1.0, radiance.size)
            else:
                costheta = np.asarray(costheta, dtype=np.float64)

            if radiance.size != costheta.size:
                raise ValueError(
                    'The sizes of the radiance and costheta array must be equal!')

            radiance_f = interp1d(costheta, radiance*np.sqrt(1.0 - costheta**2))
            ct_range = costheta.min(), costheta.max()

            if meth == 'simps':
                npts = max(int((max(2*n, npts)//2))*2 + 1, 3)
                ct = np.linspace(ct_range[0], ct_range[1], max(n, npts))
                cdf = np.zeros([int(ct.size//2) + 1])
                radiance_ct = radiance_f(ct)
                dx = (ct[-1] - ct[0])/(ct.size - 1)
                cdf[1:] = dx/3.0*(
                    radiance_ct[:-2:2] + 4.0*radiance_ct[1:-1:2] + radiance_ct[2::2])
                cdf = cdf.cumsum()
                cdf /= cdf[-1]
                lut = interp1d(cdf, ct[::2])(lut_random)
            elif meth == 'quad':
                ct = np.linspace(ct_range[0], ct_range[1], max(n, npts))
                cdf = np.zeros_like(ct)
                for index, ct_item in enumerate(ct):
                    cdf[index] = quad(radiance_f, ct_range[0], ct_item)[0]
                cdf /= cdf[-1]
                lut = interp1d(cdf, ct)(lut_random)

            super().__init__(lut, 0.0, 1.0)


class CollectionLut(LinearLut):
    def __init__(self, sensitivity: np.ndarray or str or 'CollectionLut',
                 costheta: np.ndarray=None, n: int = 1000):
        '''
        A lookup table of detector sensitivity (collection efficiency)
        as a function of the angle of incidence. Azimuthal symmetry of the
        collection sensitivity is assumed.

        Parameters
        ----------
        sensitivity: np.ndarray or str or CollectionLut
            Detector sensitivity (efficiency) as a function of the
            angle of incidence. If the value is a str, an instance is loaded
            from the file represented by the string.
        costheta: np.ndarray
            Angle of incidence cosines at which the detector
            sensitivities are defined. If None, a vector of uniformly
            distributed values from 0 to 1 is used.
        n: int
            The size of lookup table that will be used in Monte Carlo
            simulations.
        '''
        if isinstance(sensitivity, (str, CollectionLut)):
            super().__init__(sensitivity)
        else:
            sensitivity = np.asarray(sensitivity, dtype=np.float64)
            if costheta is None:
                costheta = np.linspace(0.0, 1.0, sensitivity.size)
            else:
                costheta = np.asarray(costheta)

            ct = np.linspace(costheta.min(), costheta.max(), n)
            lut = interp1d(costheta, sensitivity)(ct)
            first, last = ct[0], ct[-1]

            super().__init__(lut, first, last)


class LutEntry:
    def __init__(self, manager: 'LutManager', data: np.ndarray, offset: int):
        '''
        Lookup table entry.

        Parameters
        ----------
        manager: LutManager
            Manager of this lookup table entry.
        data: np.ndarray
            Lookup table data as a numpy array.
        offset: int
            Offset of the first element of this lookup table in the common
            lookup table buffer.
        '''
        self._manager = manager
        self._data = data
        self._offset = int(offset)

    def _get_offset(self) -> int:
        return self._offset
    offset = property(_get_offset, None, None,
                      'Offset of the first element of this lookup table '
                      'entry from the beginning of the data buffer.')

    def _get_manager(self) -> 'LutManager':
        return self._manager
    manager = property(_get_manager, None, None,
                       'Lookup table manager that manages this entry.')

    def _get_data(self) -> np.ndarray:
        return self._data
    data = property(_get_manager, None, None, 'Data array of the entry.')

    def __eq__(self, other):
        if isinstance(other, LutEntry):
            return id(self.data) == id(other.data)
        else:
            return id(self.data) == id(other)

    def __ne__(self, other):
        return not self.__eq__(other)

class LutManager:
    '''
    Manages the data of multiple lookup tables. The lookup tables can be
    easily packed into a single contiguous array that is passed to an
    OpenCL kernel.
    '''

    def __init__(self, dtype):
        '''
        Creates a new instance of the lookup table manager.

        Parameters
        ----------
        dtype: np.dtype
            Data type of the array that will contain all the managed
            lookup table arrays. 
        '''
        self._data = []
        self._offsets = {}
        self._size = 0
        self._dtype = np.dtype(dtype)

    def _get_dtype(self) -> np.dtype:
        return self._dtype
    dtype = property(_get_dtype, None, None, 'Lookup table data type.')

    def append(self, lut: np.ndarray, force: bool = False) \
            -> Tuple[LutEntry, bool]:
        '''
        Appends a new lookup table to the list of managed lookup tables and
        returns the offset from the first element of the contiguous array of
        all managed lookup tables.

        A lookup table is appended only if an existing lookup table with the
        same values is not found in the list of managed lookup tables.

        Parameters
        ----------
        lut: np.ndarray vector
            A new lookup table to append to the list of managed lookup tables.
            The new lookup table is appended only if the same lookup table i
            not found in the list of existing lookup tables.
        force: bool
            If True, the lookup table is appended without comparing its content
            to the list of existing managed lookup tables. This improves the
            performance but can grow the total size of the contiguous
            data array due to duplicate lookup tables in the managed list. 

        Returns
        -------
        lut: LutEntry
            Lookup table entry representing the given data.
        added: bool
            True if a lookup table entry with the same data was not
            found among the existing managed lookup table entries.
            If the value of the force argument is True,
            this value is always returned as True.
        '''
        offset = 0
        existing_lut_found = False

        if not force:
            for existing_lut in self._data:
                if np.array_equal(existing_lut, lut):
                    existing_lut_found = True
                    break
                offset += existing_lut.size
        else:
            offset = self._size

        if not existing_lut_found:
            self._data.append(lut)
            self._size += lut.size

        lut_id = id(lut)
        if lut_id not in self._offsets:
            self._offsets[lut_id] = offset

        return LutEntry(manager=self, data=lut, offset=offset), not existing_lut_found

    def __contains__(self, lut: np.ndarray or LutEntry) -> bool:
        '''
        Returns True, if the lookup table has been already added to the managed
        list of lookup tables.
        
        Note that the Python id is used to identify lookup table object
        and not the lookup table content which is used wen appending a new
        lookup table.

        Parameters
        ----------
        lut: np.ndarray
            Lookup table.

        Returns
        -------
        contains: bool
            Returns True, if the lookup table has been already added to the
            object.
        '''
        if isinstance(lut, LutEntry):
            return id(lut.data) in self._offsets
        else:
            return id(lut) in self._offsets

    def offset(self, lut: np.ndarray or LutEntry) -> int:
        '''
        Get the offset of the first element in the lookup table.

        Parameters
        ----------
        lut: np.ndarray
            Lookup table array as passed to the append method.

        Returns
        -------
        offset:int
            Offset of the first element or None if the lookup table
            is not found.
        '''
        if isinstance(lut, LutEntry):
            return self._offsets.get(id(lut.data))
        else:
            return self._offsets.get(id(lut))

    def clear(self):
        '''
        Clear the lookup table manager - remove all the managed data entries.
        '''
        self._data = []
        self._offsets = {}
        self._size = 0

    def __len__(self):
        return self._size

    def pack_into(self, target: np.ndarray = None, dtype: np.dtype = None) \
                  -> np.ndarray:
        '''
        Move the lookup tables into a contiguous numpy array.
        If target is None, a new target array is created using the data type
        specified in dtype. If the target array is too small to fit all the
        data, a new target of the same data type is allocated.

        Parameters
        ----------
        target: np.ndarray
            A contiguous numpy array that will hold the lookup tables
            on return. One of the input arguments target and dtype must
            not be None!
        dtype: np.dtype
            This data type overloads the data type passed to the constructor.

        Returns
        -------
        data: np.ndarray
            A contiguous numpy array with all the lookup tables.
        '''
        if dtype is None:
            dtype = self.dtype

        if target is None or target.dtype != dtype or target.size < self._size:
            target = np.empty((self._size,), dtype=dtype)

        offset = 0
        for item in self._data:
            target[offset:offset + item.size] = item
            offset += item.size

        return target


class LutManagers:
    '''
    Manages multiple lookup table managers of type
    :py:class:`LutManager`.
    '''
    def __init__(self):
        self._lut_managers = {}

    def manager(self, dtype: np.dtype) -> LutManager:
        '''
        Return lookup table manager for the given data type.
        Lookup table managers are created on the first demand/use.

        Parameters
        ----------
        dtype: np.dtype
            Numpy data type of the lookup table manager.

        Returns
        -------
        allocator: LutManager
            Lookup table manager.

        Note
        ----
        Lookup table managers can be retrieved by the [] operator that takes
        the numpy data type of the lookup table manager. 
        '''
        dtype = np.dtype(dtype)
        manager = self._lut_managers.get(dtype)
        if manager is None:
            self._lut_managers[dtype] = manager = LutManager(dtype)
        return manager

    def clear(self):
        '''
        Clear the data (managed data arrays) of all the managers.
        '''
        for _, manager in self._lut_managers.items():
            manager.clear()

    def __getitem__(self, dtype: np.dtype) -> LutManager:
        return self.manager(dtype)

    def __str_(self):
        return 'LutManagers()'

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))


class RestrictedLutManagers(LutManagers):
    def __init__(self, dtypes: List[np.dtype] or Tuple[np.dtype]):
        '''
        Manages multiple lookup table managers of type
        :py:class:`LutManager`. The range of data types is restricted to
        list of dtypes.

        Parameters
        ----------
        dtypes: list or tuple
            List of data types as numpy.dtype that will be supported by the
            instance.
        '''
        dtypes = [np.dtype(item) for item in dtypes]
        self._dtypes = tuple(dtypes)
        if len(self._dtypes) != len(set(self._dtypes)):
            raise ValueError('Duplicate data types are not allowed in '
                             'restricted lookup table managers.')
        super().__init__()

    def manager(self, dtype: np.dtype) -> LutManager:
        '''
        Return lookup table manager for the given data type.
        Lookup table managers are created on the first demand/use.

        Parameters
        ----------
        dtype: np.dtype
            Numpy data type of the lookup table manager.

        Returns
        -------
        allocator: LutManager
            Lookup table manager.

        Note
        ----
        Lookup table managers can be retrieved by the [] operator that takes
        the numpy data type of the lookup table manager. 
        '''
        dtype = np.dtype(dtype)
        if dtype not in self._dtypes:
            raise TypeError(
                'Lookup table manager for the given data type is not available!')
        return super().manager(dtype)
        dtype = np.dtype(dtype)

    def __getitem__(self, dtype: np.dtype) -> LutManager:
        return self.manager(dtype)

    def __str_(self):
        return 'RestrictedLutManagers()'

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
