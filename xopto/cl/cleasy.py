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

import time
import os.path

from typing import List, Tuple

import pyopencl as cl
import numpy as np

from xopto.cl import clinfo
from xopto import ROOT_PATH, USER_TMP_PATH


class ArgScalar:
    def __init__(self, argtype: np.dtype,
                 initializer: bool or int or float or np.ndarray = None):
        '''
        Used to declare a scalar kernel argument with the
        :py:meth:`Program.declare` method.

        Parameters
        ----------
        argtype: numpy.dtype
            Argument data type.
        initializer: int, float, scalar numpy types
            Initializer for the OpenCL buffer.
        '''
        self._type = argtype
        self._access = ''
        self._buffer = None
        if initializer is not None:
            initializer = argtype(initializer)
        self._initializer = initializer

    def initializer(self) -> bool or int or float or np.ndarray:
        '''
        Returns the initializer passed to the constructor.
        '''
        return self._initializer

    def access(self) -> str:
        '''
        Access mode of the argument ("r", "w", "rw" or "") as seen by the
        kernel (NOT by the host/python).
        '''
        return self._access

    def todevice(self, ctx: cl.Context, queue: cl.CommandQueue, value=None,
                 buffer: cl.Buffer = None, argind: int = None,
                 verbose: bool = False):
        '''
        Transfer data to OpenCL device.

        This method is called by the :py:class:`Program` instance before
        executing an OpenCL kernel and should not be called directly.

        Parameters
        ----------
        ctx: cl.Context
            OpenCL context.
        queue: cl.CommandQueue
            OpenCL command queue.
        value: bool or int or float or np.ndarray
            Value that is used to update the argument before it is passed to
            the kernel. Use existing value if None.
        buffer: cl,Buffer
            Target OpenCL device buffer or None to create a new one.
        argind: int
            Index of this argument.
        verbose:
            Turn on verbose reporting.

        Returns
        -------
        buffer: cl.Buffer
            The input buffer or a new buffer that is properly initialized
            with the data.
        '''
        if value is not None:
            self._buffer = self._type(value)
        buffer = self._buffer

        return buffer

class ArgArray:
    @staticmethod
    def _clMemFlags(rwstr):
        mf = cl.mem_flags
        return {'r':mf.READ_ONLY, 'w':mf.WRITE_ONLY,
                'rw':mf.READ_WRITE, 'wr':mf.READ_WRITE}[rwstr]

    def __init__(self, argtype: np.dtype, access: str = "rw",
                 initializer: np.ndarray = None):
        '''
        Used to declare an array kernel argument with the 
        :py:meth:`Program.declare` method.

        Parameters
        ----------
        argtype: np.dtype
            Argument data type as a numpy dtype.
        access: str
            Access specifier for the argument. Must be one of "r", "w" or "rw".
            Note that the access flags are defined as seen by the OpenCL kernel
            (not as seen by the host/python).
        initializer: np.ndarray
            Initializer for the OpenCL buffer.
        '''
        self._type = np.dtype(argtype)
        self._access = str(access).lower()
        if self._access not in ("r", "w", "rw"):
            raise ValueError('Argument access flags mus be one of '
                             '"r", "w" or "rw"!')
        if initializer is not None:
            initializer = np.asarray(initializer, dtype=self._type)
        self._initializer = initializer

        self._clbuffer = None

    def access(self) -> str:
        '''
        Access flags of the variable ("r", "w", "rw" or "") as seen by the
        OpenCL kernel (NOT as seen by the host/python)
        '''
        return self._access

    def initializer(self) -> np.ndarray:
        '''
        Returns the initializer passed to the constructor.
        '''
        return self._initializer


    def todevice(self, ctx: cl.Context, queue: cl.CommandQueue,
                 npdata: np.ndarray, clbuffer=None, argind: int =None,
                 verbose: bool = False) -> cl.Buffer:
        '''
        Transfer data to OpenCL device. 
        
        This method is called by the :py:class:`Program` instance before
        executing an OpenCL kernel and should not be called directly

        Parameters
        ----------
        ctx: cl.Context
            OpenCL context.
        queue: cl.CommandQueue
            OpenCL queue.
        npdata: np.ndarray
            Value that is used to update the argument before it is passed to
            the kernel. Use existing value if None.
        buffer: cl.Buffer
            Target OpenCL device buffer or None to create a new one.
        argind: int
            Index of this argument.
        verbose: bool
            Turn on verbose reporting.

        Returns
        -------
        buffer: cl.Buffer
            The input buffer or a new buffer that is properly initialized
            with the data.
        '''
        if clbuffer is None:
            clbuffer = self._clbuffer
            # if we have an existing buffer and array is None, we are done
            if clbuffer is not None and npdata is None:
                return clbuffer

        copyflag = 0
        if 'r' in self._access and npdata is not None:
            copyflag = cl.mem_flags.COPY_HOST_PTR

        if npdata is not None and self._type != npdata.dtype:
            raise TypeError('Kernel argument no. {} should be an array '\
                'of type {}!'.format(argind + 1, self._type))

        if clbuffer is None or \
                npdata is not None and clbuffer.size != npdata.nbytes:
            if verbose:
                print('Kernel argument no. {} requires a new '\
                      'device buffer allocation.'.format(argind + 1))

            if clbuffer is not None:
                clbuffer.release()

            if copyflag:
                clbuffer = cl.Buffer(ctx, self._clMemFlags(self._access) |
                                     copyflag, hostbuf=npdata)
            else:
                clbuffer = cl.Buffer(ctx, self._clMemFlags(self._access),
                                     size=npdata.nbytes)

        else:
            if copyflag:
                if verbose:
                    print('Kernel argument no. {} is '\
                          'being initialized with data.'.format(argind + 1))
                cl.enqueue_copy(queue, clbuffer, npdata)

        # update the local buffer
        self._clbuffer = clbuffer

        return clbuffer

    def fromdevice(self, queue: cl.CommandQueue, npdata: np.ndarray,
                   clbuffer: cl.Buffer) -> np.ndarray:
        '''
        Transfer data from OpenCL device into a numpy buffer.

        This method is called by the :py:class:`Program` instance after
        executing an OpenCL kernel and should not be called directly.

        Parameters
        ----------
        queue: cl.CommandQueue
            OpenCL queue.
        npdata: np.ndarray
            Target numpy array.
        buffer: cl.Buffer
            OpenCL buffer that will be transferred to the npdata numpy array.

        Returns
        -------
        npdata: np.ndarray
            The input numpy array filled with the contents of the OpenCL buffer. 

        Note
        ----
        Data will be transferred only if the OpenCL buffer has a write flag
        ("w"). OpenCL buffers with no write flag cannot be modified and there
        is no point in transfering the data.
        '''
        if 'w' in self._access:
            cl.enqueue_copy(queue, npdata, clbuffer)

class ArgLocalMemory:
    def __init__(self, byte_size: int):
        '''
        Used to declare a local memory buffer kernel argument
        with the Program.declare method.

        Parameters
        ----------
        byte_size: int
            Size of the local buffer in bytes.
        '''
        self._byte_size = int(byte_size)
        self._access = ''
        self._initializer = None
        self._clbuffer = None

    def access(self) -> str:
        '''
        Returns the argument access mode passed to the constructor.
        Note that local memory buffers do not support access attributes.
        '''
        return self._access

    def initializer(self) -> np.ndarray:
        '''
        Returns the initializer passed to the constructor.
        Note that initializers are not supported by local memory buffers.
        '''
        return self._initializer

    def todevice(self, ctx: cl.Context = None, queue: cl.CommandQueue=None,
                 nparray: np.ndarray = None, clbuffer: cl.Buffer = None,
                 argind: int = None, verbose: bool = False) -> cl.Buffer:
        '''
        Transfer data to the OpenCL device.

        This method is called by the :py:class:`Program` instance before
        executing an OpenCL kernel and should not be called directly

        Parameters
        ----------
        ctx: cl.Context
            OpenCL context.
        queue: cl.CommandQueue
            OpenCL queue.
        nparray: np.ndarray
            Value that is used to update the argument before it is passed to
            the kernel. Use existing value if None.
            Cannot be used with local OpenCL memory type.
        clbuffer: cl.Buffer
            OpenCL device buffer that will be used or None.
        verbose: bool
            Turn on verbose reporting.
        argind: int
            Index of the kernel argument.

        Returns
        -------
        buffer: cl.Buffer
            OpenCL device buffer.
        '''
        if self._clbuffer is None:
            self._clbuffer = cl.LocalMemory(self._byte_size)

        return self._clbuffer


class Program:
    DEFAULT_CL_DIR = os.path.join(ROOT_PATH, 'cl', 'kernel')
    ''' Place to look for the OpenCL kernels related to this class. '''

    DEFAULT_EXPORT_FILE = os.path.join(USER_TMP_PATH, 'clprogram.c')
    ''' Default filename suffix for auto-generated OpenCL source code ''' 

    with open(os.path.join(DEFAULT_CL_DIR, 'clbase.h'), 'rt') as fid:
        BASE_LIB = fid.read() + '\n'
    ''' Load the required OpenCL code on import. '''

    @staticmethod
    def clbase(doubleprecision: bool = False) -> str:
        '''
        Return the class-related OpenCL source code in the requested precision.

        Parameters
        ----------
        doubleprecision: bool
            Turns on double precision if set to True.

        Returns
        -------
        clcode: str
            OpenCL code in requested precision.
        '''
        if doubleprecision:
            return '#define USE_DOUBLE_PRECISION\n\n' + Program.BASE_LIB
        else:
            return Program.BASE_LIB

    def __init__(self, code: str,
                 device: str or list or tuple or cl.Device = None,
                 buildopts: list = [], verbose: bool = False,
                 exportsrc: bool = False):
        '''
        Creates a new instance of an OpenCL program.

        Parameters
        ----------
        code: str
            OpenCL code as a string.
        device: str or list or tuple or cl.Device
            A list of devices on which the code is to be executed. Can use
            the same arguments as with the :py:func:`xopto.cl.clinfo.device`
            function.
        buildopts: list
            A list of build options passed to the OpenCL compiler.
        verbose: bool
            If set to True, additional information is displayed during the
            build and kernel calls.
        exportsrc, bool, src
            A path (default is used if True) where the target OpnCL source
            code is saved.
        '''
        if device is None:
            self._ctx = cl.create_some_context()
        else:
            if isinstance(device, str) or \
                    (isinstance(device, tuple) or  isinstance(device, list)) and \
                    isinstance(device[0], str):
                device = clinfo.device(device)
            if not isinstance(device, tuple) and not isinstance(device, list):
                device = [device]
            self._ctx = cl.Context(device)
            if verbose:
                print('Using OpenCL device(s): {}'.format(self._ctx.devices))

        self._queue = cl.CommandQueue(self._ctx)

        if exportsrc:
            if isinstance(exportsrc, bool):
                exportsrc = self.DEFAULT_EXPORT_FILE
            with open(exportsrc, 'wt') as fid:
                fid.write(code)

        if verbose:
            t1 = time.perf_counter()
        self._clprg = cl.Program(self._ctx, code).build(options=buildopts)
        self._kernels = self._clprg.kernel_names.split(';')
        if verbose:
            print('Code build in {:.1f} ms.'.format(
                (time.perf_counter() - t1)*1000.0))
            print('Kernels in the compiled source: {}.'.format(self._kernels))

        self._declarations = {}

        self._verbose = bool(verbose)

    def kernelinfo(self, name: str) -> List[dict]:
        '''
        Returns information about the arguments of the specified kernel.

        Parameters
        ----------
        kernel: str
            Kernel name as in the OpenCL source file.

        Returns
        -------
        info: List[dict]
            A list of dicts with information on the kernel arguments.
        '''
        kernel = getattr(self._clprg, name)

        info = []
        for ind in range(kernel.num_args):
            aq = kernel.get_arg_info(ind, cl.kernel_arg_info.ADDRESS_QUALIFIER)
            aq = cl.kernel_arg_address_qualifier.to_string(aq)

            acq = kernel.get_arg_info(ind, cl.kernel_arg_info.ACCESS_QUALIFIER)
            acq = cl.kernel_arg_access_qualifier.to_string(acq)

            tn = kernel.get_arg_info(ind, cl.kernel_arg_info.TYPE_NAME)

            an = kernel.get_arg_info(ind, cl.kernel_arg_info.NAME)

            tq = kernel.get_arg_info(ind, cl.kernel_arg_info.TYPE_QUALIFIER)
            tq = cl.kernel_arg_type_qualifier.to_string(tq)

            info.append({'address':aq, 'access':acq, 'type':tn,
                         'qualifier':tq, 'name':an})

        return info

    def device(self) -> cl.Device:
        '''
        Return the OpenCL context device.
        '''
        return self._ctx.devices[0]

    def declare(self, kernel: str,
                args: List[ArgScalar or ArgArray or ArgLocalMemory]) \
                    -> List[ArgScalar or ArgArray or ArgLocalMemory]:
        '''
        Use this function to declare a kernel by name. Use the same name as in
        the OpenCL code.

        Parameters
        ----------
        kernel: str
            Kernel name as in the OpenCL file.
        args: list/tuple of arguments
            Define a list of arguments in terms of :py:class:`ArgScalar`,
            :py:class:`ArgArray` or :py:class:`ArgLocalMemory` instances.

        Returns
        -------
        decl: list
            Declaration of kernel arguments.
        '''
        if kernel not in self._kernels:
            raise ValueError('Kernel "{}" does not exist in the '\
                'compiled source!'.format(kernel))
        clBuffers = [None]*len(args)
        for ind, arg in enumerate(args):
            if arg.initializer() is not None:
                clBuffers[ind] = arg.todevice(
                    self._ctx, self._queue, arg.initializer(),
                    clBuffers[ind], verbose=self._verbose, argind=ind)
        self._declarations[kernel] = [kernel, args, clBuffers]

        return self._declarations[kernel]

    def __getattr__(self, name):
        if name in self._declarations:
            return lambda *args, **kwargs: self._exec(name, *args, **kwargs)
        elif name in self._kernels:
            raise RuntimeError(
                'Kernel "{}" was not declared yet! '
                'Declare a kernel by invoking the "declare" method!'.format(
                    name)
            )
        else:
            raise AttributeError(
                "'Program' object has no attribute '{}'".format(str(name)))

        return None

    def _exec(self, kernel: str,
              args: List[ArgScalar or ArgArray or ArgLocalMemory],
              globalwg: int or Tuple[int] or List[int],
              localwg: int or Tuple[int] or List[int] = None):
        '''
        Executes an OpenCL kernel.

        Parameters
        ----------
        kernel: str
            OpenCL kernel name that will be executed.
        args: List[ArgScalar or ArgArray or ArgLocalMemory]
            A list of kernel arguments. If None is used for a particular
            argument, the existing buffer content is used by the kernel. If
            an OpenCL buffer cannot be created an error is raised.

            Numpy arrays of arguments with "r" or "rw" access passed in the
            list will be transferred to the OpenCL device before executing
            the kernel.

            Numpy arrays of arguments with "w" or "rw" access passed in the
            list will be updated with data from the OpenCL after executing the
            kernel.
        globalwg: int or Tuple[int] or List[int]
            Global work group size.
        localwg: int or Tuple[int] or List[int]
            Local work group size. If None, an optimal value is selected by
            the OpenCL runtime.
        '''
        if self._verbose:
            t1_exec = time.perf_counter()

        if isinstance(globalwg, int):
            globalwg = [globalwg]

        if isinstance(localwg, int):
            localwg = [localwg]

        kernel, declaredArgs, clBuffers = self._declarations[kernel]

        if self._verbose:
            t1 = time.perf_counter()
        # transfer data to the device
        for ind, declaredArg in enumerate(declaredArgs):

            clBuffers[ind] = declaredArg.todevice(
                self._ctx, self._queue, args[ind], clBuffers[ind],
                verbose=self._verbose, argind=ind)

            if clBuffers[ind] is None:
                raise RuntimeError('Kernel argument no. {} '\
                                   'is undefined'.format(ind + 1))
        if self._verbose:
            print('Arguments of kernel "{}" uploaded to OpenCL device '
                  'in {:.1f} ms.'.format(
                      kernel, (time.perf_counter() - t1)*1000.0))

        if self._verbose:
            t1 = time.perf_counter()
        getattr(self._clprg, kernel)(self._queue, globalwg, localwg,
                                     *clBuffers)
        if self._verbose:
            print('Kernel "{}" executed in {:.1f} ms.'.format(
                kernel, (time.perf_counter() - t1)*1000.0))

        if self._verbose:
            t1 = time.perf_counter()
        # read data from the device
        for ind, arg in enumerate(args):
            if 'w' in declaredArgs[ind].access() and arg is not None:
                declaredArgs[ind].fromdevice(self._queue,
                                             arg, clBuffers[ind])
        if self._verbose:
            print('Arguments of kernel "{}" downloaded from OpenCL device '
                  'in {:.1f} ms.'.format(
                      kernel, (time.perf_counter() - t1)*1000.0))

        if self._verbose:
            print('Total python + kernel execution time {:.1f} ms.'.format(
                1000.0*(time.perf_counter() - t1_exec)))

    def upload_kernel_arg(self, kernel: str, arg_index: int,
                          np_data: np.ndarray):
        '''
        Updates the contents of the OpenCL buffer for the given kernel
        argument.

        Parameters
        ----------
        kernel: str
            A declared OpenCL kernel name.
        arg_inex: int
            Index of the kernel argument.
        np_data: npy.ndarray
            Numpy data to transfer to the kernel argument.
        '''
        kernel, declaredArgs, clBuffers = self._declarations[kernel]

        declaredArgs[arg_index].todevice(
            self._ctx, self._queue, np_data, clBuffers[arg_index],
            verbose=self._verbose, argind=arg_index)

if __name__ == "__main__":
    import os

    from xopto.cl import clrng

    os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"

    # Make a kernel code that executes the random number generator in
    # parallel.
    code = \
    '''
    __kernel void rng(__global uint64_t *X, __global uint32_t const *A,
                      __global float *rnd){
      /* load the random number generator state for this thread */
      uint64_t x = X[get_global_id(0)];
      uint32_t a = A[get_global_id(0)];

      rnd[get_global_id(0)] = mlrngccf(&x, a);

      /* save the random number generator state of this thread */
      X[get_global_id(0)] = x;
    };


    __kernel void test(cl_fp_t k, __global cl_fp_t *in_out){
        in_out[get_global_id(0)] *= k;
    };
    '''

    # Load the OpenCL source code, select OpenCL device and build the program
    k = Program(Program.clbase() + code, verbose=1, exportsrc=True,
                device=['nvidia', 'amd', 'hd', 'cpu'])

    # Get some seeds for the Opencl random number generator. 
    rng = clrng.Random()
    x, a = rng.seeds() # get the maximum number of seeeds
    rngdata = np.empty([rng.maxseeds], dtype=np.float32)

    # Declare the kernel arguments.
    # 1st argument - seed X array (the OpenCL kernel reeds and updates these seeds)
    # 2nd argument - seed a array (the OpenCL kernel only reads these seeds)
    # 3rd argument - empty numpy array for the random numbers (the OpenCL kernel will write into this array)
    k.declare('rng', [ArgArray(np.uint64, 'rw', x),
                      ArgArray(np.uint32, 'r', a), ArgArray(np.float32, 'w')])

    # Call the random number generator.
    k.rng([None, None, rngdata], rngdata.size)
    print(rngdata.min(), rngdata.max())

    #%%
    np_in = np.random.rand(1000).astype(np.float32)

    k.declare('test', [ArgScalar(np.float32, 10.0),
                       ArgArray(np.float32, 'rw')])

    k.test([1/10, np_in], np_in.size)
