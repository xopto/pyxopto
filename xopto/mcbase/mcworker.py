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

import time

import numpy as np

from xopto.mcbase import cltypes
from xopto.mcbase import mcobject
from xopto.mcbase import mctypes
from xopto.mcbase import mcoptions
from xopto.cl import clinfo
from xopto.cl import clrng

from xopto.mcbase.mcutil.buffer import \
    RestrictedBufferAllocators, BufferAllocator, BufferAllocation, \
    NumpyAllocators, NumpyAllocator

from xopto.mcbase.mcutil.lut import \
    RestrictedLutManagers, LutManager, LutEntry

import pyopencl as cl


class ClWorker(mcobject.McObject):
    def __init__(self,
                 types: mctypes.McDataTypesBase = mctypes.McDataTypesSingle,
                 cl_devices: str or cl.Device or List['cl.Device'] 
                             or cl.Context or cl.CommandQueue = None,
                 cl_build_options: List[str] = None,
                 cl_profiling: bool = False,
                 **kwargs):
        '''
        OpenCL worker with support for allocation of read-write data buffers
        of arbitrary types and management of read-only lookup tables.

        Parameters
        ----------
        types: mctypes.McDataTypes
            A class that defines all the simulator data types.
            Use one of the following predefined type classes derived from
            mctypes.McDataTypes:

            - mctypes.McDataTypesSingle

                - 32-bit size type
                - 32-bit default integers,
                - 64-bit detector accumulators,
                - single precision floating-point arithmetics,
                - 32-bit photon packet counter (maximum number of photon
                  packets per OpenCL kernel call limited to 4,294,967,295)

            - mctypes.McDataTypesDouble

                - 32-bit size type,
                - 32-bit default integers,
                - 64-bit detector accumulators,
                - double precision floating-point arithmetics,
                - 32-bit photon packet counter (maximum number of photon
                  packets per OpenCL kernel call limited to 4,294,967,295)

            - mctypes.McDataTypesSingleCnt64

                - 64-bit size type,
                - 32-bit default integers,
                - 64-bit detector accumulators,
                - single precision floating-point arithmetics,
                - 64-bit photon packet counter (maximum number of photon
                  packets per OpenCL kernel call virtually unlimited)

            - mctypes.McDataTypesDoubleCnt64

                - 64-bit size type,
                - 32-bit default integers,
                - 64-bit detector accumulators,
                - double precision floating-point arithmetics,
                - 64-bit photon packet counter (maximum number of photon
                  packets per OpenCL kernel call virtually unlimited)

        cl_devices: str or cl.Device or List[cl.Device] or cl.Context or cl.CommnadQueue
            A python list of OpenCL devices that are used for
            conducting the simulation. See the clGpuDevices and clCpuDevices
            functions of the :py:mod:`xopto.clinfo` module for details on
            obtaining a list of OpenCl capable devices. If None is provided,
            the first available device is used (GPU devices have priority over
            CPU devices). Use function :py:func:`xopto.clinfo.device` to get
            a desired device by simple keywords, e.g. a call
            clinfo.device(['amd', 'nvidia', 'hd', 'cpu'], that will search
            for an AMD GPU, Nvidia GPU, Intel Hd GPU, any CPU and return
            the first device found.
            The value of this input argument can also be an
            instance of OpenCL Context or an instance of OpenCL CommandQueue.
            Note that in case an instance of CommandQueue is passed, the value
            of parameter cl_profiling is ignored since it is not possible to
            enable or disable profiling on an existing OpenCL CommandQueue.

        cl_build_options: List[str]
            A list of OpenCL build option as specified by the OpenCl manuals at
            https://www.khronos.org/.
            An example of commonly used build options:
            cloptions=['-cl-opt-disable', '-Werror', '-cl-fast-relaxed-math', '-cl-mad-enable'].

        cl_profiling: bool
            Enables OpenCL command queue profiling.
            Note that in case an instance of CommandQueue is passed as the
            cl_devices parameter, the value of parameter cl_profiling is
            ignored since it is not possible to enable or disable profiling
            on an existing OpenCL CommandQueue.

        kwargs: dict
            Keyword arguments passed to the Mixin classes.

        Note
        ----
        Use the :py:meth:`create_allocator` method to create buffer allocators
        and the :py:meth:`create_lut_manager` method to create managers of
        lookup tables.
        '''
        # Data types that will be used by this simulator instance.
        if not issubclass(types, mctypes.McDataTypesBase):
            raise TypeError(
                'The types argument must be a subclass of McDataTypes!')
        self._types = types

        # Allocators for read-write OpenCL buffers that are passed to the
        # OpenCL kernels as a single buffer.
        self._cl_rw_allocators = RestrictedBufferAllocators(
            dtypes=(types.np_accu, types.np_int, types.np_float))
        # OpenCL buffers that are produced by the allocators of read-write
        # OpenCL buffers.
        self._cl_rw_allocators_buffers = {}

        # Allocators for temporary numpy data buffers.
        self._np_allocators = NumpyAllocators()

        # Storage for all the read-only lookup table managers used by this object.
        self._r_lut_managers= RestrictedLutManagers(
            dtypes=(types.np_int, types.np_float))

        # Numpy data buffers.
        self._np_buffers = {}
        # OpenCL data buffers.
        self._cl_buffers = {}

        self._cl_devices = self._cl_queue = self._cl_context = None
        # Save the list of OpenCL devices
        if isinstance(cl_devices, cl.Device):
            self._cl_devices = [cl_devices]
        elif isinstance(cl_devices, str):
            self._cl_devices = [clinfo.device(cl_devices)]
        elif isinstance(cl_devices, cl.Context):
            self._cl_context = cl_devices
            self._cl_devices = self._cl_context.devices
        elif isinstance(cl_devices, cl.CommandQueue):
            self._cl_queue = cl_devices
            self._cl_context = self._cl_queue.context
            self._cl_devices = self._cl_context.devices
            cl_profiling = bool(
                self._cl_queue.properties & 
                cl.command_queue_properties.PROFILING_ENABLE)

        self._cl_profiling = bool(cl_profiling)
        cl_cq_properties = None
        if self._cl_profiling:
            cl_cq_properties = cl.command_queue_properties.PROFILING_ENABLE

        # OpenCL build options
        if cl_build_options is None:
            cl_build_options = []
        self._cl_build_options = [str(item) for item in cl_build_options]

        # The OpenCL context used by the worker.
        if self._cl_context is None:
            self._cl_context = cl.Context(self._cl_devices)
        # The OpenCL queue used by the worker.
        if self._cl_queue is None:
            self._cl_queue = cl.CommandQueue(
                self._cl_context, properties=cl_cq_properties)
        # The latest executable build with this the worker.
        self._cl_exec = None

        # a dict lookup for type names from numpy dtype
        self._dtype_to_type_str = {
            np.dtype('uint64'): 'uint64', np.dtype('int64'): 'int64',
            np.dtype('uint32'): 'uint32', np.dtype('int32'): 'int32',
            np.dtype('uint16'): 'uint16', np.dtype('int16'): 'int16',
            np.dtype('uint8'): 'uint8', np.dtype('int8'): 'int8',
            np.dtype('float32'): 'float', np.dtype('float64'): 'double'
        }

    def event_timing(self, event: cl.Event) -> Tuple[float, float]:
        '''
        If OpenCL profiling is enabled, the timing related to an OpenCL event
        can be retrieved through the profile property of the event:

        - ev.profile.queued   - nanosecond counter captured on event queued
        - ev.profile.submit   - nanosecond counter captured on event submitted
        - ev.profile.start    - nanosecond counter captured on start of execution
        - ev.profile.complete - nanosecond counter captured on end of execution

        This method returns the time (s) that was required to start executing
        of the command and the time (s) required to execute the command.

        Parameters
        ----------
        event: cl.Event
            OpenCL event instance.

        Returns
        -------
        dt_delay: float
            Time (s) that was required to start the execution of the command
            (profile.start - profile.submit).
        dt_exec: float
            Time (s) that was required to execute the command
            (profile.complete - profile.start).
        '''
        if self.cl_profiling:
            return  (event.profile.start - event.profile.submit)*1e-9, \
                    (event.profile.complete - event.profile.start)*1e-9

    def cl_build(self, cl_src: str, verbose: bool = False) -> cl.Program:
        '''
        Build OpenCL source code. A context and command queue are created on
        the first run.

        Parameters
        ----------
        cl_src: str
            OpenCL source code as a string.
        verbose: bool
            Turns on verbose reporting.

        Returns
        -------
        program: cl.Program
            OpenCL executable.
        '''
        if not isinstance(cl_src, str):
            raise TypeError('The openCl source code must be a string!')

        if self._cl_context is None:
            self._cl_context = cl.Context(self._cl_devices)
        if self._cl_queue is None:
            properties = None
            if self._cl_profiling:
                properties = cl.command_queue_properties.PROFILING_ENABLE
            self._cl_queue = cl.CommandQueue(
                self._cl_context, properties=properties)

        if verbose:
            print('Executing OpenCL code on: {}'.format(
                self._cl_context.devices))
            print('OpenCL build options:', self._cl_build_options)

        tb = time.perf_counter()
        cl_exec = cl.Program(self._cl_context, cl_src).build(
            options=self._cl_build_options)
        # options=['-cl-opt-disable', '-Werror']
        buildtime = time.perf_counter() - tb

        if verbose:
            print('Source code built in {:.3f} ms.'.format(buildtime*1000.0))

        return cl_exec

    def _get_cl_device(self) -> List[cl.Device]:
        return self._cl_context.devices
    cl_device = property(_get_cl_device, None, None,
                         'OpenCL device that is used to run this simulator '
                         'instance. ')

    def _get_cl_context(self) -> cl.Context:
        return self._cl_context
    cl_context = property(_get_cl_context, None, None, 'OpenCL context.')

    def _get_cl_build_options(self) -> Tuple[str]:
        return self._cl_build_options
    cl_build_options = property(_get_cl_build_options, None, None,
                                'Returns a tuple of OpenCL build options '
                                'that were passed to the constructor.')

    def _get_cl_profiling(self) -> bool:
        return self._cl_profiling
    cl_profiling = property(_get_cl_profiling, None, None,
                            'Returns True if OpenCL command queue '
                            'allows profiling.')

    def _get_cl_queue(self) -> cl.CommandQueue:
        return self._cl_queue
    cl_queue = property(_get_cl_queue, None, None, 'OpenCL command queue.')

    def _get_cl_exec(self) -> cl.Program:
        return self._cl_exec
    cl_exec = property(_get_cl_exec, None, None,
                       'The latest OpenCL program built with this worker.')

    def dtype_to_typename(self, dtype: np.dtype) -> str:
        '''
        Return a standard short type name for the given numpy data type.

        Parameters
        ----------
        dtype: np.dtype
            Numpy data type.

        Returns
        -------
        typename: str
            A short standard type name for the given numpy data type.
        '''
        name = self._dtype_to_type_str.get(np.dtype(dtype))

        if name is None:
            raise TypeError('The provided numpy data type is not supported!')

        return name

    def np_allocators(self) -> NumpyAllocators:
        '''
        Get allocators of temporary numpy buffers.

        Returns
        -------
        allocators: NumpyAllocators
            Numpy allocators of temporary buffers.
        '''
        return self._np_allocators

    def cl_rw_allocator(self, dtype:np.dtype) -> BufferAllocator:
        '''
        Returns OpenCL read-write buffer allocator for the give data type.

        Parameters
        ----------
        dtype: np.dtype
            Numpy data type used with the allocator.

        Returns
        -------
        allocator: BufferAllocator
            Buffer allocator for the given data type.
        '''
        return self._cl_rw_allocators[dtype]

    def cl_allocate_rw_buffer(self, dtype: np.dtype, owner: any, shape:
                              tuple, download=True):
        '''
        Allocate a read-write OpenCL buffer of the given data type using
        the related OpenCL buffer allocator.

        Parameters
        ----------
        dtype: np.dtype
            Numpy data type of the buffer allocator.
        owner: any
            Object that will own the allocated buffer.
        shape: Tuple[int]
            Shape of the buffer array to allocate.
        download: bool
            Set to True if the data buffer should be downloaded after
            executing the kernel.
            The downloaded data will be passed to the :py:meth:`update_data`
            of the owner.

        Returns
        -------
        allocation: BufferAllocation
            Buffer allocation object with information on the allocation.
        '''
        return self._cl_rw_allocators[dtype].allocate(
            owner=owner, shape=shape, download=download)

    def cl_rw_allocator_buffer(self, dtype: np.dtype, fill=None) -> cl.Buffer:
        '''
        Fetch the OpenCL read-write buffer of the allocator for the given
        data type.
        A new OpenCL buffer is created on the first call. A new OpenCL buffer
        is also allocated if the size of the existing buffer is too small for
        the allocations.
        The fill argument can be used to initialize the buffer with the given
        value.

        Parameters
        ----------
        dtype: np.dtype
            Buffer allocator data type.
        fill: np.dtype.type
            If not None, initialize the buffer with the given value / fill.

        Returns
        -------
        buffer: cl.Buffer
            OpenCL read-write buffer of the allocator.
        '''
        allocator = self._cl_rw_allocators[dtype]
        cl_buffer = self._cl_rw_allocators_buffers.get(allocator)
        nbytes = allocator.dtype.itemsize*allocator.size
        if cl_buffer is None or cl_buffer.size < nbytes:
            # an OpenCL buffer does not exist yet - create one
            if nbytes > 0:
                cl_buffer = cl.Buffer(
                        self._cl_context, cl.mem_flags.READ_WRITE, nbytes)
                self._cl_rw_allocators_buffers[allocator] = cl_buffer

        if fill is not None and nbytes > 0:
            self.cl_w_buffer_fill(cl_buffer, dtype, fill)

        return cl_buffer

    def cl_w_buffer_fill(self, cl_buffer:cl.Buffer, dtype:np.dtype,
                         fill: int or float, offset: int = 0, size: int = None,
                         cl_kernel=None):
        '''
        Fast fill of writable OpenCL buffers with a given scalar value. 

        The OpenCL buffer must be writable, since the initialization is
        performed in an OpenCL kernel!

        Parameters
        ----------
        cl_buffer: cl.Buffer
            The OpenCL buffer to fill.
        dtype: np.dtype
            Numpy dtype of a buffer element/item.
        fill: int or float
            Scalar value used as a buffer fill. Must be convertible to
            the given data type (dtype).
        size: int
            The last item of the buffer relative to the offset that will be
            filled. If None, the entire buffer from the offset will be filled.
            Note that the size is given in buffer items not bytes!
        offset: int
            First item of the buffer that will be filled.
            Note that the size is given in buffer items not bytes!
        cl_kernel: cl.Kernel
            OpenCL kernel that will be used to fill the buffer. If None,
            a matching kernel will be searched in the current executable
            :py:attr:`cl_exec`. The fill kernels must follow the
            following footprint:

            .. code-block::

                fill_<dtype>(__global T *buffer, T fill_value,
                             mc_size_t size, mc_size_t offset)

            where dtype should be one of (double, float, int64, uint64,
            int32 or uint32).

        Note
        ----
        Note that the size and offset are given in buffer items not bytes!
        '''
        dtype = np.dtype(dtype)
        if cl_kernel is None:
            type_name = self.dtype_to_typename(dtype)
            kernel_name = 'fill_{}'.format(type_name)
            cl_kernel = getattr(self.cl_exec, kernel_name, None)
            if cl_kernel is None:
                raise TypeError('No fill kernel is defined for the given '
                                'data type {}!'.format(dtype))

        cl_size = int(cl_buffer.size/dtype.itemsize)

        if offset >= cl_size:
            raise IndexError('Fill offset exceeds the OpenCL buffer size!')
        if size is None:
            size = cl_size - offset
        elif size + offset > cl_size:
            raise IndexError('Fill range (offset: offset + size) exceeds '
                             'the OpenCL buffer size!')

        dtype_sizet = np.dtype(self.types.np_size)
        cl_kernel(self.cl_queue, (cl_size, ), None, cl_buffer, dtype.type(fill),
                  dtype_sizet.type(size), dtype_sizet.type(offset)).wait()


    def r_lut_manager(self, dtype: np.dtype) -> LutManager:
        '''
        Return a read-only lookup table manager that matches the given
        data type.

        Parameters
        ----------
        dtype: np.dtype
            Data type managed by the lookup table manager.
        '''
        return self._r_lut_managers[dtype]

    def append_r_lut_data(
        self, dtype: np.dtype, data:np.ndarray,
        force: bool = False) -> LutEntry:
        '''
        Add a read-only lookup table data array to the list of managed
        lookup table arrays.

        Parameters
        ----------
        dtype: np.dtype
            Data type of the lookup table manager.
        data: np.ndarray
            Lookup table data to be added to the managed list. The numpy
            array must be of a type that can be converted to the data
            type of the lookup table manager.
        force: bool
            Add the array to the managed list even if an existing entry for the
            given data array is found.

        Returns
        -------
        lutentry: LutEntry
            Lookup table entry. Use the :py:attr:LutEntry.offset
            to get the offset of the first element in the common data array.
        '''
        return self._r_lut_managers[dtype].append(data, force=force)

    def cl_r_lut_buffer(self, dtype: np.dtype, update=True):
        '''
        Returns an OpenCL buffer (create if not exist) representing the
        read-only data arrays managed by the related read-only lookup
        table manager.
        The value of update argument controls if the content of the OpenCL
        buffer is updated with the content of the managed arrays.
        If an OpenCL buffer does not exist, a new one is created and updated
        regardless of the value of the update argument.

        Parameters
        ----------
        dtype: np.dtype
            Data type of the lookup table manager.
        update: bool
            If True, update the flat numpy array  with the
            current content of the managed lookup tables.
            If an OpenCL buffer does not exist, a new is created and updated
            regardless of the value of the update argument.

        Returns
        -------
        cl_float_lut: cl.Buffer
            OpenCL buffer of the floating-point lookup table array.
        '''
        lut_manager = self._r_lut_managers[dtype]
        cl_buffer = self._cl_buffers.get(lut_manager)
        if cl_buffer is None or update:
            np_buffer = self._np_buffers.get(lut_manager)
            np_buffer = lut_manager.pack_into(np_buffer)
            self._np_buffers[lut_manager] = np_buffer 

            cl_buffer = self._get_cl_lut_buffer(
                lut_manager, np_buffer=np_buffer, access='r')

        return cl_buffer

    def _get_cl_lut_buffer(
            self, lut_manager: BufferAllocator,
            np_buffer: np.ndarray = None, access: str = 'r') -> cl.Buffer:
        '''
        Internal method that returns an existing or creates a new numpy data
        array-based read-write or read-only OpenCL buffer.

        Parameters
        ----------
        lut_manager: str
            Lookup table manager.
        np_buffer: np.ndarray
            Numpy data buffer used as initializer.
        access: 'r' or 'rw'
            Access flags for the OpenCL buffer.

        Returns
        -------
        buffer: cl.Buffer
            OpenCL buffer.
        '''
        # allocate and initialize the buffer
        cl_buffer = self._cl_buffers.get(lut_manager)
        
        mf = {'rw': cl.mem_flags.READ_WRITE,
              'r': cl.mem_flags.READ_ONLY}.get(access, 'rw')
        mf_cp_host = cl.mem_flags.COPY_HOST_PTR
    
        if cl_buffer is None:
            # the OpenCL buffer does not exist yet - create one
            if np_buffer is not None and np_buffer.size > 0:
                self._cl_buffers[lut_manager] = cl.Buffer(
                        self._cl_context, mf | mf_cp_host,
                        hostbuf=np_buffer)
            else:
                self._cl_buffers[lut_manager] = cl.Buffer(
                    self._cl_context, mf,
                    lut_manager.dtype.itemsize)
        else:
            if np_buffer is not None and np_buffer.size > 0:
                if np_buffer.nbytes == cl_buffer.size:
                    cl.enqueue_copy(
                        self._cl_queue,
                        cl_buffer,
                        np_buffer).wait()
                else:
                    print('Allocating a new lut OpenCL Buffer!', cl_buffer)
                    self._cl_buffers[lut_manager] = cl.Buffer(
                        self._cl_context, mf | mf_cp_host,
                        hostbuf=np_buffer)

        return self._cl_buffers[lut_manager]

    def cl_rw_buffer(self, name:str, data: np.ndarray = None,
                     size: int = None) -> cl.Buffer:
        '''
        Create a read-write OpenCL buffer and optionally initialize the buffer
        with data.

        Parameters
        ----------
        name: str
            A unique name for the buffer.
        data: any object with a buffer interface
            Initializer for the OpenCl buffer.
        size: int
            Size of the buffer in bytes.

        Returns
        -------
        buffer: cl.Buffer
            The OpenCL buffer.

        Note
        ----
        One of the size or data input arguments must be defined. If both
        are defined, data will be used to derive the size and initialize
        the OpenCL buffer.
        '''
        cl_buffer = self._cl_buffers.get(name)
        if cl_buffer is None:
            if size is None and data is None:
                raise ValueError(
                    'Cannot create an OpenCL buffer if the initialization '
                    'data and the buffer size are unknown!'
                )
            if data is not None:
                cl_buffer = cl.Buffer(
                    self._cl_context,
                    cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR,
                    hostbuf=data)
            else:
                cl_buffer = cl.Buffer(
                    self._cl_context, cl.mem_flags.READ_WRITE, size)
            self._cl_buffers[name] = cl_buffer

        elif data is not None and data.size > 0:
            if cl_buffer.size < data.nbytes:
                cl_buffer = cl.Buffer(
                    self._cl_context,
                    cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR,
                    hostbuf=data)
                self._cl_buffers[name] = cl_buffer
            else:
                cl.enqueue_copy(
                        self._cl_queue,
                        cl_buffer,
                        data).wait()

        return cl_buffer

    def cl_w_buffer(self, name:str, size: int = None) -> cl.Buffer:
        '''
        Create a write only OpenCL buffer and optionally initialize the buffer
        with data. If a buffer of equal or larger size is found for the
        given name, a new buffer is not created.

        Parameters
        ----------
        name: str
            A unique name for the buffer.
        size: int
            Size of the buffer in bytes.

        Returns
        -------
        buffer: cl.Buffer
            The opencl buffer.
        '''
        cl_buffer = self._cl_buffers.get(name)
        if cl_buffer is None:
            if size is None:
                raise ValueError(
                    'Cannot create an OpenCL buffer if the size of the buffer'
                    'is unknown!'
                )
            else:
                cl_buffer = cl.Buffer(
                    self._cl_context, cl.mem_flags.WRITE_ONLY, size)
                self._cl_buffers[name] = cl_buffer

        else:
            if cl_buffer.size < size:
                cl_buffer = cl.Buffer(
                    self._cl_context, cl.mem_flags.WRITE_ONLY, size)
                self._cl_buffers[name] = cl_buffer

        return cl_buffer

    def cl_r_buffer(self, name: str, data: np.ndarray,
                    size: int = None) -> cl.Buffer:
        '''
        Create a read-only OpenCL buffer and optionally initialize the buffer
        with data. If a buffer of equal or larger size is found for the
        given name, a new buffer is not created.

        Parameters
        ----------
        name: str
            A unique name for the buffer.
        data: any object with a buffer interface
            Initializer for the OpenCl buffer.
        size: int
            Size of the buffer in bytes.
        buffer: cl.Buffer
            The OpenCL buffer.

        Note
        ----
        One of the size or data input arguments must be defined. If both
        are defined, data will be used to derive the size and initialize
        the OpenCL buffer.
        '''
        cl_buffer = self._cl_buffers.get(name)
        if cl_buffer is None:
            if size is None and data is None:
                raise ValueError(
                    'Cannot create an OpenCL buffer if the initialization '
                    'data and the buffer size are unknown!'
                )
            if data is not None:
                cl_buffer = cl.Buffer(
                    self._cl_context,
                    cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                    hostbuf=data
                )
            else:
                cl_buffer = cl.Buffer(
                    self._cl_context, cl.mem_flags.READ_ONLY, size
                )
            self._cl_buffers[name] = cl_buffer

        elif data is not None:
            size = self._sizeof(data)
            if size > 0:
                if size != cl_buffer.size:
                    cl_buffer = cl.Buffer(
                        self._cl_context,
                        cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR,
                        hostbuf=data)
                    self._cl_buffers[name] = cl_buffer
                else:
                    cl.enqueue_copy(
                            self._cl_queue,
                            cl_buffer,
                            data).wait()

        return cl_buffer

    def _sizeof(self, data):
        if isinstance(data, np.ndarray):
            return data.nbytes
        else:
            return cltypes.sizeof(data)

    def cl_allocation_download(self, allocation: BufferAllocation,
                               out: np.ndarray = None) -> np.ndarray:
        '''
        Download the content of the allocated buffer from the OpenCL device.

        Parameters
        ----------
        allocation: BufferAllocation
            Buffer allocation as returned by the
            :py:meth:`allocate_buffer` methods.
        out: np.ndarray
            Optional target numpy data array that will be filled with the
            content of the OpenCL buffer.
            The type and size of the numpy array must match the allocation.

        Returns
        -------
        out: np.ndarray
            Numpy data array that is filled with the OpenCL buffer content.
        '''
        if allocation.allocator not in self._cl_rw_allocators:
            raise ValueError('Buffer allocation was not made by this object!')

        offset = allocation.offset*allocation.dtype.itemsize
        if out is None:
            out = np.empty(allocation.shape, dtype=allocation.dtype)
        # a permissive check
        elif out.nbytes != allocation.size*allocation.dtype.itemsize:
            raise ValueError(
                'The output buffer does not match the type and/or size of '
                'the allocation!')
        
        cl.enqueue_copy(
            self._cl_queue, out, self._cl_buffer_from_allocation(allocation),
            device_offset=offset
        ).wait()

        return out

    def cl_allocation_upload(self, allocation: BufferAllocation,
                             data: np.ndarray):
        '''
        Upload an allocated OpenCL buffer from host to the OpenCL device.

        Parameters
        ----------
        allocation: BufferAllocation
            Buffer allocation as returned by the
            :py:meth:`cl_allocate_rw_int_buffer` or
            :py:meth:`cl_allocate_rw_accumulator_buffer` or
            :py:meth:`cl_allocate_rw_accumulator_buffer` methods.
        data: np.ndarray
            Data array to be uploaded. The type and size of the numpy array
            must match the allocation.
        '''
        if allocation.allocator not in self._cl_rw_allocators:
            raise ValueError('Buffer allocation was not made by this object!')

        offset = allocation.offset*allocation.dtype.itemsize
        if allocation.dtype != data.dtype or allocation.size != data.size:
            raise ValueError(
                'The data buffer does not match the type and/or size of '
                'the allocation!')
        
        cl.enqueue_copy(
            self._cl_queue, self._cl_buffer_from_allocation(allocation), data,
            device_offset=offset
        ).wait()

    def _cl_buffer_from_allocation(
            self, allocation: BufferAllocation) -> cl.Buffer:
        '''
        Internal method that fetches an OpenCl buffer related to the
        given allocation.

        Parameters
        ----------
        allocation: BufferAllocation
            Buffer allocation object.

        Returns
        -------
        cl_buffer: cl.Buffer
            OpenCL buffer related to the allocation.
        '''
        return self._cl_rw_allocators_buffers[allocation.allocator]

    def _get_cl_rw_allocators(self) -> RestrictedBufferAllocators:
        return self._cl_rw_allocators
    cl_rw_allocators = property(_get_cl_rw_allocators, None, None,
                                'Allocators of OpenCL buffers.')

    def _get_np_allocators(self) -> NumpyAllocators:
        return self._np_allocators
    np_allocators = property(_get_np_allocators, None, None,
                             'Allocators of temporary numpy buffers.')

    def _get_types(self) -> mctypes.McDataTypesBase:
        return self._types
    types = property(_get_types, None, None,
                     'Data types used by the OpenCL kernel.')

    def _get_np_buffers(self) -> dict:
        return self._np_buffers
    np_buffers = property(_get_np_buffers, None, None, 'Numpy data buffers.')

    def _get_cl_buffers(self) -> dict:
        return self._cl_buffers
    cl_buffers = property(_get_cl_buffers, None, None, 'OpenCL data buffers.')


class ClWorkerStandardBufferLutMixin:
    '''
    ClWorker class mixin that creates interfaces for a
    standard set of OpenCL read-write buffer allocators and read-ony lookup
    table managers.
    '''

    def _get_float_r_lut_manager(self) -> LutManager:
        return self.r_lut_manager(self.types.np_float)
    float_r_lut_manager = property(
        _get_float_r_lut_manager, None, None,
        'Floating-point read-only data lookup table manager.')

    def _get_int_r_lut_manager(self) -> LutManager:
        return self.r_lut_manager(self.types.np_int)
    int_r_lut_manager = property(
        _get_int_r_lut_manager, None, None,
        'Integer read-only data lookup table manager.')

    def _get_cl_rw_accumulator_allocator(self) -> BufferAllocator:
        return self.cl_rw_allocator(self.types.np_accu)
    cl_rw_accumulator_allocator = property(
        _get_cl_rw_accumulator_allocator, None, None,
        'Allocator of read-write accumulator type OpenCL buffers.')

    def _get_cl_rw_float_allocator(self) -> BufferAllocator:
        return self.cl_rw_allocator(self.types.np_float)
    cl_rw_float_allocator = property(
        _get_cl_rw_float_allocator, None, None,
        'Allocator of read-write floating-point type OpenCL buffers.')

    def _get_cl_rw_int_allocator(self) -> BufferAllocator:
        return self.cl_rw_allocator(self.types.np_int)
    cl_rw_int_allocator = property(
        _get_cl_rw_int_allocator, None, None,
        'Allocator of read-write integer type OpenCL buffers.')

    def cl_r_float_lut(self, update: bool = True) -> cl.Buffer:
        '''
        Returns an OpenCL buffer (create if not exist) representing the
        read-only data arrays managed by the floating-point lookup table
        manager.
        The value of update argument controls if the content of the OpenCL
        buffer is updated with the content of the managed arrays.
        If an OpenCL buffer does not exist, a new is created and updated
        regardless of the value of the update argument.

        Parameters
        ----------
        update: bool
            If True, update the flat numpy array of lookup tables with the
            current content of the managed lookup tables.
            If an OpenCL buffer does not exist, a new one is created and
            updated regardless of the value of the update argument.

        Returns
        -------
        cl_float_lut: cl.Buffer
            OpenCL buffer of the floating-point lookup table array.
        '''
        return self.cl_r_lut_buffer(self.types.np_float, update=update)

    def cl_r_int_lut(self, update: bool = True) -> cl.Buffer:
        '''
        Returns an OpenCL buffer (create if not exist) representing the
        read-only integer arrays managed by the integer lookup table
        manager.
        The value of update argument controls if the content of the OpenCL
        buffer is updated with the content of the managed arrays.
        If an OpenCL buffer does not exist, a new one is created and updated
        regardless of the value of the update argument.

        Parameters
        ----------
        update: bool
            If True, update the flat numpy array of lookup tables with the
            current content of the managed lookup tables.
            If an OpenCL buffer does not exist, a new one is created and
            updated regardless of the value of the update argument.

        Returns
        -------
        cl_float_lut: cl.Buffer
            OpenCL buffer of the floating-point lookup table array.
        '''
        return self.cl_r_lut_buffer(self.types.np_int, update=update)

    def cl_rw_accumulator_buffer(self, fill=None) -> cl.Buffer:
        '''
        Fetch the OpenCL buffer of the accumulator allocator.
        A new OpenCL buffer is created on the first call.

        Parameters
        ----------
        fill: int
            Scalar value used to fill the buffer. Must be convertible to
            the given data type (dtype).

        Returns
        -------
        buffer: cl.Buffer
            OpenCL buffer.
        '''
        return self.cl_rw_allocator_buffer(self.types.np_accu, fill=fill)

    def cl_rw_float_buffer(self, fill: float = None) -> cl.Buffer:
        '''
        Fetch the OpenCL buffer of the floating-point allocator.
        A new OpenCL buffer is created on the first call.

        Parameters
        ----------
        fill: float
            Scalar value used to fill the buffer. Must be convertible to
            the given data type (dtype).

        Returns
        -------
        buffer: cl.Buffer
            OpenCL buffer.
        '''
        return self.cl_rw_allocator_buffer(self.types.np_float, fill=fill)

    def cl_rw_int_buffer(self, fill: int = None) -> cl.Buffer:
        '''
        Fetch the OpenCL buffer of the integer allocator.
        A new OpenCL buffer is created on the first call.

        Parameters
        ----------
        fill: int
            Scalar value used to fill the buffer. Must be convertible to
            the given data type (dtype).

        Returns
        -------
        buffer: cl.Buffer
            OpenCL buffer.
        '''
        return self.cl_rw_allocator_buffer(self.types.np_int, fill=fill)

    def cl_allocate_rw_accumulator_buffer(self, owner: any, shape: Tuple[int],
                                          download: bool = True) \
            -> BufferAllocation:
        '''
        Allocate a new read-write accumulator buffer with OpenCL type mc_accu_t.

        Parameters
        ----------
        owner: any
            Allocation owner / caller.
        shape: Tuple[int]
            Allocation buffer shape.
        download: bool
            Set to True if download is required after executing the kernel.
            The downloaded data will be passed to the :py:meth:`update_data`
            method of the owner.

        Returns
        -------
        allocation: BufferAllocation
            The buffer allocation.
        '''
        return self.cl_allocate_rw_buffer(self.types.np_accu, owner, shape,
                                          download=download)

    def cl_allocate_rw_float_buffer(self, owner: any, shape: tuple,
                                    download: bool = True) -> BufferAllocation:
        '''
        Allocate a new read-write floating-point buffer with OpenCL type mc_fp_t.

        Parameters
        ----------
        owner: any
            Allocation owner / caller.
        shape: Tuple[int]
            Allocation buffer shape.
        download: bool
            Set to True if download is required after executing the kernel.
            The downloaded data will be passed to the :py:meth:`update_data`
            method of the owner.

        Returns
        -------
        allocation: BufferAllocation
            The buffer allocation.
        '''
        return self.cl_allocate_rw_buffer(self.types.np_float, owner, shape,
                                          download=download)

    def cl_allocate_rw_int_buffer(self, owner: any, shape: tuple,
                                  download: bool = True) -> BufferAllocation:
        '''
        Allocate a new read-write integer buffer with OpenCL type mc_int_t.

        Parameters
        ----------
        owner: any
            Allocation owner / caller.
        shape: Tuple[int]
            Allocation buffer shape.
        download: bool
            Set to True if download is required after executing the kernel.
            The downloaded data will be passed to the :py:meth:`update_data`
            method of the owner.

        Returns
        -------
        allocation: BufferAllocation
            The buffer allocation.
        '''
        return self.cl_allocate_rw_buffer(self.types.np_int, owner, shape,
                                          download=download)

    def append_r_int_lut(self, data:np.ndarray, force: bool = False) -> LutEntry:
        '''
        Append data to the integer type read-only lookuptable with
        The OpenCL data type of the read-only lookup table is mc_int_t.

        Parameters
        ----------
        data: np.ndarray
            Lookup table data to manage.
        force: bool
            If True, append the data array even if it is already managed.

        Returns
        -------
        lut_entry: LutEntry
            Lookup table entry. Use the offset property to locate the
            first element of the entry in the lookup table buffer.
        new: bool
            True if a lookup table entry with the same data was not
            found among the existing managed lookup table entries.
            If the value of the force argument is True,
            this value is always returned as True.
        '''
        return self.append_r_lut_data(self.types.np_int, data, force=force)

    def append_r_float_lut(self, data:np.ndarray, force: bool = False) \
            -> LutEntry:
        '''
        Append data to the floating-point type read-only lookuptable.
        The OpenCL data type of the read-only lookup table is mc_float_t.

        Parameters
        ----------
        data: np.ndarray
            Lookup table data to manage.
        force: bool
            If True, append the data array even if it is already managed.

        Returns
        -------
        lut_entry: LutEntry
            Lookup table entry. Use the offset property to locate the
            first element of the entry in the lookup table buffer.
        new: bool
            True if a lookup table entry with the same data was not
            found among the existing managed lookup table entries.
            If the value of the force argument is True,
            this value is always returned as True.
        '''
        return self.append_r_lut_data(self.types.np_float, data, force=force)

    def clear_r_float_lut(self):
        '''
        Clear the content of the floating point read-only lookup table manager.
        '''
        self.float_r_lut_manager.clear()

    def clear_r_int_lut(self):
        '''
        Clear the content of the integer read-only lookup table manager.
        '''
        self.int_r_lut_manager.clear()

    def _get_np_r_float_lut(self) -> np.ndarray:
        return self.np_buffers[self.r_lut_manager(self.types.np_float)]
    np_r_float_lut = property(
        _get_np_r_float_lut, None, None,
        'Numpy data buffer of the read-only floating-point type lookup table.')

    def _get_np_r_int_lut(self) -> np.ndarray:
        return self.np_buffers[self.r_lut_manager(self.types.np_int)]
    np_r_int_lut = property(
        _get_np_r_int_lut, None, None,
        'Numpy data buffer of the read-only integer type lookup table.')


class ClWorkerRngMixin:
    '''
    ClWorker class mixin that can be used to create an OpenCl Random number
    generator.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Random number generator mixin for an OpenCL worker. Use the
        :py:attr:`rng` to access the random number generator and
        :py:attr:`cl_max_threads` to determine the maximum number
        of threads that can be run concurrently with this random number
        generator. The mutable seeeds can be accessed through the
        :py:attr:`rng_seeds_x` property and the immutable seeds through
        the :py:attr:`rng_seeds_a`.

        Parameters
        ----------
        rnginit: np.uint64
            OpenCL random number generator initializer as a 64-bit unsigned
            integer.
            Use this initializer if there is a need to put the random, number
            generator into a known initial state. 
        '''
        super().__init__(*args, **kwargs)

        rnginit = kwargs.get('rnginit')
        # Initialization of the Monte Carlo random number generator.
        self._rng = rng = clrng.Random()
        if rnginit is not None:
            rnginit = np.uint64(rnginit)
        rng_seeds_x, rng_seeds_a = rng.seeds(rng.maxseeds, xinit=rnginit)
        self.np_buffers['rng_seeds_x'] = rng_seeds_x
        self.np_buffers['rng_seeds_a'] = rng_seeds_a

        # Maximum number of concurrent OpenCL threads that can be run with the
        # random number generator (total number of available seeds).
        self._cl_max_threads = int(rng.maxseeds//512)*512

    def _get_rng_seeds_x(self):
        return self.np_buffers['rng_seeds_x']
    rng_seeds_x = property(_get_rng_seeds_x, None, None,
                           'Mutable random number generator seeds X.')

    def _get_rng_seeds_a(self):
        return self.np_buffers['rng_seeds_a']
    rng_seeds_a = property(_get_rng_seeds_a, None, None,
                           'Immutable random number generator seeds A.')

    def _get_cl_max_threads(self) -> int:
        return self._cl_max_threads
    cl_max_threads = property(_get_cl_max_threads, None, None,
                             'The maximum number of OpenCl threads that can '
                             'be run concurrently with this random number '
                             'generator.')

    def _get_rng(self) -> clrng.Random:
        return self._rng
    rng = property(_get_rng, None, None, 'Random number generator for OpenCL.')


if __name__ == '__main__':
    TEST_CODE = '\n'.join((
        '__kernel void test(int a, int b){',
        '   int c = a + b;',
        '}',
    ))
    class Worker(ClWorkerStandardBufferLutMixin, ClWorker):
        def __init__(self,
                 types: mctypes.McDataTypesBase = mctypes.McDataTypesSingle,
                 options: List[mcoptions.McOption] = None,
                 cl_devices=None, cl_build_options=None):
            super().__init__(types, options, cl_devices, cl_build_options)


    cw = Worker(cl_devices='hd')
    prog = cw.cl_build(TEST_CODE, verbose=True)
