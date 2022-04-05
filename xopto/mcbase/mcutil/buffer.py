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

from typing import Tuple, List

import numpy as np

from xopto.mcbase import mcobject


class BufferAllocation:
    def __init__(self, offset: int, size: int, shape: tuple,
                 dtype: np.dtype, owner: any, initializer: np.ndarray = None,
                 allocator: 'BufferAllocator' = None, download=True):
        '''
        Buffer allocation instance.

        Parameters
        ----------
        offset: int
            Offset of the first element of the allocated buffer from the
            start of the common memory buffer (in number of elements).
        size: int
            Size of the allocated buffer in number of items.
        shape: tuple
            Shape of the allocated buffer.
        dtype: np.dtype
            Numpy data type of the allocated buffer.
        owner: any
            Owner of the allocation.
        initializer: np.ndarray or None
            Initializer for the allocation.
        allocator: BufferAllocator
            Buffer allocator that created this buffer.
        download: bool
            Set to True if the data buffer should be downloaded after
            executing the kernel.
            The downloaded data will be passed to the :py:meth:`update_data`
            of the owner.
        '''
        self._offset = int(offset)
        self._shape = tuple(shape)
        self._size = int(size)
        self._dtype = np.dtype(dtype)
        self._owner = owner
        self._allocator = allocator
        self._download = bool(download)
        if initializer is not None:
            initializer = np.ascontiguousarray(initializer, dtype=self._dtype)
        self._initializer = initializer

    def free(self):
        '''
        Clear the allocation.
        '''
        self._allocator.free(self)
        self._allocator = None

    def _get_owner(self) -> any:
        return self._owner
    owner = property(_get_owner, None, None, 'Allocation owner.')

    def _get_allocator(self) -> 'BufferAllocator':
        return self._allocator
    allocator = property(_get_allocator, None, None, 'Allocation allocator.')

    def _get_initializer(self) -> np.ndarray or None:
        return self._initializer
    initializer = property(_get_initializer, None, None, 'Buffer initializer.')

    def _get_offset(self) -> int:
        return self._offset
    offset = property(_get_offset, None, None, 'Offset of the buffer start.')

    def _get_size(self) -> int:
        return self._size
    size = property(_get_size, None, None, 'Size of the allocated buffer.')

    def _get_dtype(self) -> np.dtype:
        return self._dtype
    dtype = property(_get_dtype, None, None, 'Numpy data type as dtype.')

    def _get_shape(self) -> tuple:
        return self._shape
    shape = property(_get_shape, None, None, 'Buffer shape.')

    def _get_download(self) -> bool:
        return self._download
    download = property(_get_download, None, None,
                        'Download required after execution.')

    def __str__(self):
        return 'BufferAllocation(offset={}, size={}, shape={}, dtype={}, '\
               'owner={}, initializer={}, allocator={}, download={})'.format(
                   self.offset, self._size, self.shape, self.dtype,
                   self.owner, self.initializer, self.allocator, self.download)

    def __repr__(self):
        return '{:s} # id={}'.format(self.__str__(), id(self))


class BufferAllocator:
    def __init__(self, dtype=np.float64):
        self._dtype = np.dtype(dtype)
        self._allocations = {}
        self._size = 0

    def _get_dtype(self) -> np.dtype:
        return self._dtype
    dtype = property(_get_dtype, None, None, 'Allocator data type.')

    def allocate(self, owner: any, shape: Tuple[int], download=True) -> int:
        '''
        Allocate a new data buffer for the owner. The allocations
        cannot be cleared individually. Each owner can allocate
        multiple buffers.

        Parameters
        ----------
        owner: any
            Owner of the allocated buffer.
        shape: tuple
            Shape of the data buffer.
        download: bool
            Set to True if the data buffer should be downloaded after
            executing the kernel.
            The downloaded data will be passed to the :py:meth:`update_data`
            of the owner.

        Returns
        -------
        allocation: BufferAllocation
            Buffer allocation object.
        '''
        allocations = self._allocations.get(owner)
        if allocations is None:
            allocations = []
            self._allocations[owner] = allocations

        # create a new allocation
        size = np.prod(shape)
        allocation = BufferAllocation(
            offset=self._size, shape=tuple(shape),
            size=int(size), dtype=self._dtype,
            owner=owner, allocator=self,
            download=download
        )

        allocations.append(allocation)
        self._size += size

        return allocation

    def allocations(self, owner: any) -> List[dict]:
        '''
        Get a list of existing allocations for the owner.

        Parameters
        ----------
        owner: any
            Allocation owner.

        Returns
        -------
        allocation: List[dict]
            List of allocation dicts with the following keys:
                - size: int
                    Allocation buffer size.
                - shape: tuple
                    Allocation buffer shape
                - offset: int
                    Offset of the first buffer element in the common data
                    buffer array.
                - owner: any
                    Allocation owner.

            The returned value is an empty list if no allocations for
            the given owner are found.
        '''
        return self._allocations.get(owner, [])

    def extract(self, data: np.ndarray, owner: any) -> List[np.ndarray]:
        '''
        Extract buffers of the owner from the full memory buffer.

        Parameters
        ----------
        data: np.ndarray
            Common data buffer.
        owner: any
            Find and return buffers that were allocated by this owner.

        Returns
        -------
        owner_buffers: List[np.ndarray]
            Buffers that were allocated by the owner.
        '''
        buffers = []
        for allocation in self.allocations(owner):
            buffer = data[allocation.offset: allocation.offset + allocation.size]
            buffer.shape = allocation.shape
            buffers.append(buffer)

        return buffers

    def free(self, allocation: BufferAllocation):
        '''
        Free a buffer allocation. This will not change the total
        buffer size. The memory space will be left unused.
        '''
        if allocation in self:
            self._allocations.get(allocation.owner).remove(allocation)
            allocation.free()

    def __contains__(self, allocation: BufferAllocation):
        return allocation in self._allocations.get(allocation.owner, [])

    def clear(self):
        '''
        Remove all allocations.
        '''
        for owner in self._allocations:
            for allocation in self._allocations[owner]:
                allocation.free()

        self._allocations = {}
        self._size = 0

    def create_buffer(self, mc: mcobject.McObject,
                      out: np.ndarray = None) -> np.ndarray:
        '''
        Create a new numpy buffer or return the existing buffer if the size
        and type match the requirements.

        Parameters
        ----------
        mc: McObject
            Monte Carlo simulator instance.
        out: np.ndarray
            Optional existing buffer array.

        Returns
        -------
        buffer: np.ndarray
            Existing buffer if passed as an input argument and matches the
            required size and data type, else a new buffer
        '''
        if out is None or out.size != self._size:
            out = np.zeros((self._size), dtype=self._dtype)
        return out

    def defragment(self):
        '''
        Defragment allocations. Note that this will change the location
        of buffers.
        '''
        offset = 0
        for owner, items in self._allocations.items():
            for allocation in items:
                allocation._offset = offset
                offset += allocation.size
        self._size = offset

    def _get_size(self) -> int:
        return self._size
    size = property(_get_size, None, None, 'Current buffer size.')


class BufferAllocators:
    '''
    Manages multiple general Buffer allocators of type
    :py:class:`BufferAllocator`.
    '''
    def __init__(self):
        self._allocators = {}

    def allocator(self, dtype: np.dtype) -> BufferAllocator:
        '''
        Return allocator for the given buffer data type.
        Allocators are created on the first demand/use.

        Parameters
        ----------
        dtype: np.dtype
            Numpy data type of the buffer allocator.

        Returns
        -------
        allocator: ContiguousNumpyAllocator
            Buffer allocator.

        Note
        ----
        Allocators can be retrieved by the [] operator that takes the
        numpy data type of the allocator. 
        '''
        dtype = np.dtype(dtype)
        allocator = self._allocators.get(dtype)
        if allocator is None:
            self._allocators[dtype] = allocator = BufferAllocator(dtype)
        return allocator

    def __iter__(self):
        return iter(self._allocators.values())

    def __contains__(self, item: BufferAllocator):
        return item in self._allocators.values()

    def __getitem__(self, dtype: np.dtype) -> BufferAllocator:
        return self.allocator(dtype)

    def __str_(self):
        return 'BufferAllocators()'

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))


class RestrictedBufferAllocators(BufferAllocators):
    '''
    Manages multiple general Buffer allocators of type
    :py:class:`BufferAllocator`. The buffer types are restricted to a list
    of types passed to the constructor.
    '''
    def __init__(self, dtypes: List[np.dtype]):
        dtypes = [np.dtype(item) for item in dtypes]
        self._dtypes = tuple(dtypes)
        if len(self._dtypes) != len(set(self._dtypes)):
            raise ValueError('Duplicate data types are not allowed in '
                             'restricted buffer allocators.')
        super().__init__()

    def allocator(self, dtype: np.dtype) -> BufferAllocator:
        '''
        Return allocator for the given buffer data type.
        Allocators are created on the first demand/use.

        Parameters
        ----------
        dtype: np.dtype
            Numpy data type of the buffer allocator.

        Returns
        -------
        allocator: ContiguousNumpyAllocator
            Buffer allocator.

        Note
        ----
        Allocators can be retrieved by the [] operator that takes the
        numpy data type of the allocator. 
        '''
        dtype = np.dtype(dtype)
        if dtype not in self._dtypes:
            raise TypeError(
                'Allocator for the given data type is not available!')

        return super().allocator(dtype)

    def __getitem__(self, dtype: np.dtype) -> BufferAllocator:
        return self.allocator(dtype)

    def __str_(self):
        return 'RestrictedBufferAllocators()'

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))


class NumpyAllocator:
    def __init__(self, dtype, size=1000000):
        '''
        Manages temporary allocations of contiguous numpy-based buffers.
        The initial size of the buffer grows as required by the allocations.
        The internal buffer is never released.

        Parameters
        ----------
        dtype: np.dtype
            Data type of the allocator.
        size: int
            Initial minimum size of the internal buffer created on the first
            allocation.
        '''
        self._dtype = np.dtype(dtype)
        self._buffer = None
        self._initial_size = int(size)
        self._state_positions = []
        self._pos = 0

    def allocate(self, size: int) -> np.ndarray:
        '''
        Allocates a buffer with the given size.

        Parameters
        ----------
        size: int
            Size of the allocated buffer in number of buffer items (NOT bytes).

        Returns
        -------
        buffer: np.ndarray
            Allocated numpy buffer
        '''
        size = int(size)

        if self._buffer is None:
            n = max(size, self._initial_size)
            self._buffer = np.empty((n,), dtype=self._dtype)
        elif self._buffer.size < self._pos + size:
            n = self._pos + size
            self._buffer = np.empty((n,), dtype=self._dtype)

        buffer = self._buffer[self._pos: self._pos + size]
        self._pos += size

        return buffer

    def clear(self):
        '''
        Release all allocations made by :py:meth:`allocate`.

        Note
        -----
        The internal numpy buffer is not released, only the internal buffer
        manager will clear all allocations.
        '''
        if self._state_positions:
            raise RuntimeError(
                'Cannot clear the buffer allocations within a with statement!')
        self._pos = 0

    def available(self) -> int:
        '''
        Returns the size available for allocation before a new internal
        buffer allocation is required.
        '''
        return self._size() - self._pos

    def _size(self) -> int:
        '''
        Returns the internal size of the allocated buffer
        '''
        if self._buffer is None:
            return 0
        return self._buffer.size

    def _get_size(self) -> int:
        return self._pos
    size = property(_get_size, None, None, 'Currently allocated size.')

    def _push_state(self):
        '''
        All allocations after this point will be released after calling
        th :py:meth:`_pop_state` method.
        '''
        self._state_positions.append(self._pos)
        return self

    def _pop_state(self):
        '''
        Release all buffers allocated since the last call to
        the :py:meth:`_push_state` method.
        '''
        if self._state_positions:
            self._pos = self._state_positions.pop()
        else:
            self._pos = 0

    def __enter__(self):
        self._push_state()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._pop_state()

    def __str__(self):
        return 'NumpyAllocator(dtype={}, size={})'.format(
            self._dtype, self._size)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))


class NumpyAllocators:
    def __init__(self, size: int = 0,
                 dtypes: List[np.dtype] or Tuple[np.dtype] = None):
        '''
        Manages temporary allocations of contiguous numpy-based buffers.
        The initial size of the internal numpy buffers grows as required by the
        allocations but is never released (reduced in size).

        If a list of data types (dtypes) is passed to the constructor, only
        buffers of the listed types can be allocated (restricted allocator).

        Parameters
        ----------
        size: int
            Initial minimum size of the internal buffers created on the
            first allocation.
        dtype: np.dtype
            Data types of the allocator. If None, any data type is allowed.
        '''
        self._dtypes = dtypes
        self._initial_size = int(size)
        self._pos_by_dtype = {}

        self._states = []

        self._is_dtype_allowed = None
        if dtypes is None:
            self._is_dtype_allowed = self._allow_all_dtypes
        else:
            self._is_dtype_allowed = self._allow_listed_dtypes

        self._buffer_by_dtype = {}
        if dtypes is not None:
            for dtype in dtypes:
                self._buffer_by_dtype[np.dtype(dtype)] = None 
                self._pos_by_dtype[dtype] = 0
            self._dtypes = tuple(self._buffer_by_dtype.keys())

    def _allow_all_dtypes(self, dtype=np.dtype) -> bool:
        return True

    def _allow_listed_dtypes(self, dtype=np.dtype) -> bool:
        return dtype in self._buffer_by_dtype

    def allocate(self, dtype: np.dtype, size: int) -> np.ndarray:
        '''
        Allocates a buffer of the given data type and size.
        Raises TypeError if allocating a data type that is not allowed.

        Parameters
        ----------
        dtype: np.dtype
            Data type of the buffer.
        size: int
            Size of the allocated buffer in number of buffer items (NOT bytes).

        Returns
        -------
        buffer: np.ndarray
            Allocated numpy buffer
        '''
        dtype = np.dtype(dtype)
        out = None
        if self._is_dtype_allowed(dtype):
            buffer = self._buffer_by_dtype.get(dtype)
            pos = self._pos_by_dtype.get(dtype, 0)
            if buffer is None:
                n = max(size, self._initial_size)
                buffer = np.empty((n,), dtype=dtype)
                self._buffer_by_dtype[dtype] = buffer
                pos = 0
            if buffer.size < pos + size:
                buffer = np.empty((pos + size), dtype=dtype)
                self._buffer_by_dtype[dtype] = buffer
                # pos = 0
            out = buffer[pos: pos + size]
            pos += size
            self._pos_by_dtype[dtype] = pos
        else:
            raise TypeError('The requested data type cannot be allocated!')

        return out

    def release(self, dtype: np.dtype = None):
        '''
        Release all allocations made by :py:meth:`allocate`.

        Parameters
        ----------
        dtype: np.dtype
            If not None, only the allocations of the specified data type
            will be cleared. If None, all allocations are cleared.

        Note
        ----
        The internal numpy buffers are not released, only the internal buffer
        manager will clear all allocations.
        '''
        if self._states:
            raise RuntimeError(
                'Cannot clear the buffer allocations within a with statement!')

        if dtype is None:
            for key in self._pos_by_dtype:
                self._pos_by_dtype[key] = 0
        else:
            self._pos_by_dtype[np.dtype(dtype)] = 0

    def size(self, dtype: np.dtype) -> int:
        '''
        Total size of all allocations of the given data type.

        Parameters
        ----------
        dtype: np.dtype
            Data type of the allocations.

        Returns
        -------
        size: int
            Total size of the allocations for the given data type.
        '''
        return self._pos_by_dtype.get(np.dtype(dtype), 0)

    def available(self, dtype=np.dtype) -> int:
        '''
        Returns the size available for allocation before a new internal
        buffer allocation is required for the given data type.

        Parameters
        ----------
        dtype: np.dtype
            Data type of the allocations.

        Returns
        -------
        size: int
            Total buffer size available before a new internal buffer allocation
            is required.
        '''
        dtype = np.dtype(dtype)
        if dtype in self._buffer_by_dtype:
            return self._buffer_by_dtype[dtype].size - self._pos_by_dtype[dtype]

        return 0

    def _push_state(self):
        '''
        All allocations after this point will be released after calling
        th :py:meth:`_pop_state` method.
        '''
        self._states.append(dict(self._pos_by_dtype))
        return self

    def _pop_state(self):
        '''
        Release all buffers allocated since the last call to
        the :py:meth:`_push_state` method.
        '''
        self._pos_by_dtype.update(self._states.pop())

    def __enter__(self):
        self._push_state()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._pop_state()

    def __str__(self):
        return 'NumpyAllocator(dtypes={}, size={})'.format(
            self._dtypes, self._initial_size)

    def __rpr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
