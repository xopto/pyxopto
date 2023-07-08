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
import ctypes
import sys
import platform
import os
import os.path
import time

import numpy as np
import numba as nb

from xopto import PRIMES_PATH, USER_PRIMES_PATH, USER_BIN_PATH, VERBOSE


def generate_primes(n: int, out: np.ndarray = None,
                    verbose: bool = False) -> np.ndarray:
    '''
    Generates n safe primes used by the random number generator.

    Parameters
    ----------
    n: int
        Number of safe primes to generate.
    out: np.ndarray of type np.uint64
        Optional output array for safe primes that must be of type
        np.uint64 and shape (3, n).
    verbose: bool
        Verbose mode.

    Returns
    -------
    out: np.ndarray of type np.uint64
        Output array with the safe random numbers.

    Note
    ----
    Generating a typical batch of 150,000 primes takes about 5 min. This
    function requires gmpy2 module.
    '''
    import gmpy2

    def is_prime(value):
        for n in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]:
            if not gmpy2.is_prime(gmpy2.mpz(str(value)), n):
                return False
        return True

    n = int(n)

    if out is not None:
        if not isinstance(out, np.ndarray):
            raise TypeError('The output array must be of type np.uint64!')
        if out.dtype != np.uint64:
            raise TypeError('Output array must be of type np.uint64!')
        if out.shape != (n, 3):
            raise ValueError('The shape of the output array must be (3, n)!')
    else:
        out = np.zeros((n, 3), dtype=np.uint64)

    tstart = None
    if verbose:
        print('Searcing for {:d} primes'.format(n))
        tstart = time.perf_counter()

    a = gmpy2.mpz("4294967118")
    i = 0
    while i < n:
        n2 = gmpy2.mpz(gmpy2.mul_2exp(a, 32))
        n2 = gmpy2.mpz(gmpy2.sub(n2, 1))
        if is_prime(n2):
            n1 = gmpy2.mpz(gmpy2.sub(n2, 1))
            n1 = gmpy2.mpz(gmpy2.f_div_2exp(n1, 1))
            if is_prime(n1):
                out[i] = (np.uint64(a), np.uint64(n2), np.uint64(n1))
                i += 1
                if verbose:
                    print('\rFound {}/{} primes ...'.format(i, n), end='')

        a = gmpy2.mpz(gmpy2.sub(a, 1))

    if verbose:
        print('\nDone in {:.3f} s\n'.format(time.perf_counter() - tstart))

    return out

def generate_primes_fast(n):
    cdll = None
    try:
        if platform.system() == 'Linux':
            if sys.maxsize > 0xFFFFFFFF:
                # print('64-bit os')
                dllFullName = os.path.join(USER_BIN_PATH, "rng64.so")
                cdll = ctypes.cdll.LoadLibrary(dllFullName)
            else:
                # print('32-bit os')
                dllFullName = os.path.join(USER_BIN_PATH, "rng32.so")
                cdll = ctypes.cdll.LoadLibrary(dllFullName)

        elif platform.system() == 'Windows':
            if sys.maxsize > 0xFFFFFFFF:
                # print('64-bit os')
                dllFullName = os.path.join(USER_BIN_PATH, "rng64.dll")
                cdll = ctypes.cdll.LoadLibrary(dllFullName)
            else:
                # print('32-bit os')
                dllFullName = os.path.join(USER_BIN_PATH, "rng32.dll")
                cdll = ctypes.cdll.LoadLibrary(dllFullName)

        else:
            raise RuntimeError(
                'Your operating system is not supported '
                'at the moment. \n'
                'Program is tested only on Windows and Linux'
                'operating systems. \n'
                'You must configure rng.dll or rng.so manually '
                'for your system. You can find the random numbe '
                'generator code in xopto/src/rng/rng/rng.cpp.')
    except OSError as err:
        print(err)
        if VERBOSE:
            print('Found no compatible build of the random number library.')
            print('Using pure Python implementation of the random number generator.')
            print('Run xopto.rebuild() to build the library.')

    out = a = n1 = n2 = find_primes = None

    if cdll is not None:
        # int generate_primes(
        #   unsigned int n, uint32_t *vec_a, uint64_t vec_n2, uint64_t vec_n3);
        find_primes = None
        try:
            find_primes = cdll.find_primes
            find_primes.argtypes = [
                ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_uint64)
            ]
            find_primes.restype = int
        except (OSError, AttributeError):
            pass
        if find_primes is not None:
            a = np.zeros((n,), dtype=np.uint32)
            n1 = np.zeros((n,), dtype=np.uint64)
            n2 = np.zeros((n,), dtype=np.uint64)
            n = find_primes(
                ctypes.c_uint32(n),
                a.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                n1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
                n2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64))
            )
            out = np.hstack([
                    np.reshape(a[:n], (n, 1)),
                    np.reshape(n2[:n], (n, 1)),
                    np.reshape(n1[:n], (n, 1)),
            ])
    if out is None:
        out = generate_primes(n)

    return out


@nb.jit(nopython=True)
def _rng_init(x: np.uint64, a: np.uint32, fora: np.uint32, n_rng: int, xinit: np.uint64) -> int:
    begin = fora[0]

    ui32_1 = np.uint32(1)
    ui32_0xFFFFFFFF = np.uint32(0xFFFFFFFF)
    ui64_0 = np.uint64(0)
    ui64_32 = np.uint64(32)
    ui64_0xFFFFFFFF = np.uint64(0xFFFFFFFF)
    ui64_0x100000000 = np.uint64(0x100000000)
    double_0x100000000 = np.double(ui64_0x100000000)

    if xinit == ui64_0 or \
            np.uint32(xinit >> ui64_32) >= begin - ui32_1 or \
            np.uint32(xinit) >= ui64_0xFFFFFFFF:
        return 1

    for i in range(n_rng):
        a[i] = fora[i + 1]
        x[i] = ui64_0

        while x[i] == ui64_0 or \
                np.uint32(x[i] >> ui64_32) >= fora[i + 1] - ui32_1 or \
                np.uint32(x[i]) >= ui32_0xFFFFFFFF:

            # generate a random number
            xinit = (xinit & ui64_0xFFFFFFFF)*begin + (xinit >> ui64_32)

            # Calculate c and store in the upper 32 bits of x[i].
            # Make sure 0<=c<a.
            x[i] = np.uint32(
                # np.floor(
                    np.double(np.uint32(xinit))/
                    double_0x100000000*fora[i + 1]
                # )
            )
            x[i] = x[i] << ui64_32

            # Generate a random number and store in the lower 32 bits of x[i]
            # (as the initial x of the generator)
            # x will be 0<=x<b, where b is the base 2^32
            xinit = (xinit & ui64_0xFFFFFFFF)*begin + (xinit >> ui64_32)

            x[i] += np.uint32(xinit)

    return 0

def _rand_cc(data, n, x, a):
    ui64_0xFFFFFFFF = np.uint64(0xFFFFFFFF)
    ui64_32 = np.uint64(32)
    ui32_0x7FFFFF = np.uint32(0x7FFFFF)
    single_0x7FFFFF = np.single(0x7FFFFF)

    for i in range(n):
        x = (x & ui64_0xFFFFFFFF)*a + (x >> ui64_32)
        data[i] = np.single(
                (np.uint32(x) & ui32_0x7FFFFF)
            )/single_0x7FFFFF
            
    return data, x, a

def _rand_co(data, n, x, a):
    ui64_0xFFFFFFFF = np.uint64(0xFFFFFFFF)
    ui64_32 = np.uint64(32)
    ui32_0x7FFFFF = np.uint32(0x7FFFFF)
    single_0x800000 = np.single(0x800000)

    for i in range(n):
        x = (x & ui64_0xFFFFFFFF)*a + (x >> ui64_32)
        data[i] = np.single(
                (np.uint32(x) & ui32_0x7FFFFF)
            )/single_0x800000
            
    return data, x, a


def _rand_oc(data, n, x, a):
    ui64_0xFFFFFFFF = np.uint64(0xFFFFFFFF)
    ui64_32 = np.uint64(32)
    ui32_0x7FFFFF = np.uint32(0x7FFFFF)
    single_0x7FFFFF = np.single(0x7FFFFF)

    for i in range(n):
        x = (x & ui64_0xFFFFFFFF)*a + (x >> ui64_32)
        data[i] = np.single(
                np.uint32(1) | (np.uint32(x) & ui32_0x7FFFFF)
            )/single_0x7FFFFF
            
    return data, x, a

def _rand_oo(data, n, x, a):
    ui64_0xFFFFFFFF = np.uint64(0xFFFFFFFF)
    ui64_32 = np.uint64(32)
    ui32_0x7FFFFF = np.uint32(0x7FFFFF)
    single_0x800000 = np.single(0x800000)

    for i in range(n):
        x = (x & ui64_0xFFFFFFFF)*a + (x >> ui64_32)
        data[i] = np.single(
                np.uint32(1) | (np.uint32(x) & ui32_0x7FFFFF)
            )/single_0x800000
            
    return data, x, a


class Random:
    PRIMES_FILE = 'safeprimes_base32_500k.npz'
    ''' Default npz file with primes. '''

    NPRIMES = 500000
    ''' Default number of primes to generate. '''

    def __init__(self, primes=None, forcepy=False):
        '''
        Create an instance of a random number generator.

        Class depends on external rng64.dll (rng64.so) and rng32.dll (rng32.so)
        that can be generated/built by calling :py:func:`xopto.rebuild`.
        The built shared library will be placed into
        :py:const:`xopto.USER_BIN_PATH` folder.
        If an external library is not found, a slower python implementation
        of the required functions is used instead.

        Parameters
        ----------

        primes: np.ndarray
            If provided, the primes must be organized into a numpy array
            of shape (nprimes, 2). Type of the array should be uint32 or
            uint64. See the one of the safeprimes_base32_<>.npz files
            for more details.

        forcepy: bool
            If True, use python implementation of the random number
            helper functions even if a shared library is available.

        Note
        ----
        Use the seeds method to generate seed pairs X and A for a
        a GPU-based random number generator. Each GPU thread should take one
        pair of values, e.g. X[i], A[i], to use with the random number
        generator. Note that the value of seed x is changed on each call to the
        random number generator and should be saved before the GPU thread
        returns, e.g. to the global memory. The value of seed a is not changed
        during the calls to the random number generator.

        Examples of different uniform random number generators can be found
        in xopto/cl/kernels/ml.h file.
        '''
        self._cdll = None
        if not forcepy:
            try:
                if platform.system() == 'Linux':
                    if sys.maxsize > 0xFFFFFFFF:
                        # print('64-bit os')
                        dllFullName = os.path.join(USER_BIN_PATH, "rng64.so")
                        self._cdll = ctypes.cdll.LoadLibrary(dllFullName)
                    else:
                        # print('32-bit os')
                        dllFullName = os.path.join(USER_BIN_PATH, "rng32.so")
                        self._cdll = ctypes.cdll.LoadLibrary(dllFullName)

                elif platform.system() == 'Windows':
                    if sys.maxsize > 0xFFFFFFFF:
                        # print('64-bit os')
                        dllFullName = os.path.join(USER_BIN_PATH, "rng64.dll")
                        self._cdll = ctypes.cdll.LoadLibrary(dllFullName)
                    else:
                        # print('32-bit os')
                        dllFullName = os.path.join(USER_BIN_PATH, "rng32.dll")
                        self._cdll = ctypes.cdll.LoadLibrary(dllFullName)

                else:
                    raise RuntimeError(
                        'Your operating system is not supported '
                        'at the moment. \n'
                        'Program is tested only on Windows and Linux'
                        'operating systems. \n'
                        'You must configure rng.dll or rng.so manually '
                        'for your system. You can find the random number '
                        'generator code in xopto/src/rng/rng/rng.cpp.')
            except OSError:
                if VERBOSE:
                    print('Found no compatible build of the random number library.')
                    print('Using python implementation of the random number generator.')
                    print('Run xopto.rebuild() to build the library.')

        self._init_RNG = None
        self._rand_cc = self._rand_co = self._rand_oc = self._rand_oo = None

        if self._cdll is not None:
            # create api
            # int init_RNG(uint64_t *x, uint32_t *a, uint32_t *fora,
            #   const uint32_t n_rng, uint64_t xinit){
            self._init_RNG = self._cdll.init_RNG
            #self._init_RNG.argtypes = [ctypes.POINTER(ctypes.c_uint64),
            #    ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint32),
            #    ctypes.c_uint32, ctypes.c_uint64]
            #self._init_RNG.restype = ctypes.c_int

            # int rand_cc(float *data, uint32_t n, uint64_t* x, uint32_t a)
            self._rand_cc = self._cdll.rand_cc
            #self._rand_co.argtypes = [ctypes.POINTER(ctypes.c_float),
            #    ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint64),
            #    ctypes.c_uint32]
            #self._rand_co.restype = ctypes.c_int

            # int rand_co(float *data, uint32_t n, uint64_t* x, uint32_t a)
            self._rand_co = self._cdll.rand_co
            #self._rand_co.argtypes = [ctypes.POINTER(ctypes.c_float),
            #    ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint64),
            #    ctypes.c_uint32]
            #self._rand_co.restype = ctypes.c_int

            # int rand_oc(float *data, uint32_t n, uint64_t* x, uint32_t a)
            self._rand_oc = self._cdll.rand_oc
            #self._rand_oc.argtypes = [ctypes.POINTER(ctypes.c_float),
            #    ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint64),
            #    ctypes.c_uint32]
            #self._rand_oc.restype = ctypes.c_int

            # int rand_oc(float *data, uint32_t n, uint64_t* x, uint32_t a)
            self._rand_oo = self._cdll.rand_oo
            #self._rand_oo.argtypes = [ctypes.POINTER(ctypes.c_float),
            #    ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint64),
            #    ctypes.c_uint32]
            #self._rand_oo.restype = ctypes.c_int

        if primes is None:
            fullfilename = os.path.join(PRIMES_PATH, self.PRIMES_FILE)
            primes = None
            try:
                f = np.load(fullfilename)
                primes= f['data']
                f.close()
            except Exception:
                pass
            if primes is None:
                print('Failed to load primes from {:s}!'.format(fullfilename))
                print('Generating primes ... this will take a few minutes ...')
                primes = generate_primes(self.NPRIMES, verbose=True)
                try:
                    print('Saving primes to {}'.format(fullfilename))
                    np.savez_compressed(fullfilename, primes)
                except Exception:
                    pass

            if primes.ndim == 1:
                # supporting flat format
                primes.shape = (int(primes.size//2), 2,)

        self._primes = np.asarray(primes, dtype=np.uint32)
        # must be contiguous - passed to a C function
        self._fora = np.ascontiguousarray(np.uint32(self._primes[:, 0]))

    def _get_maxseeds(self) -> int:
        return self._primes.shape[0] - 1

    maxseeds = property(
        _get_maxseeds, None, None,
        'Maximum number of random generator seeds that can be produced '
        'by the available primes.'
    )

    def initializer(self):
        '''
        Creates an initializer (32-bit unsigned integer) required to
        compute seeds (see the seeds method) of the random number generator.
        '''
        while 1:
            tmp = np.uint64(np.random.rand(2)*0xFFFFFFFF)
            xinit = tmp[0] + tmp[1]*np.uint64(0xFFFFFFFF)

            if not  ((xinit == np.uint32(0)) | \
                    (np.uint32(xinit >> np.uint64(32)) >=
                    (self._fora[0] - np.uint32(1))) | \
                    (np.uint32(xinit) >= np.uint32(0xFFFFFFFF))):
                break

        return xinit

    def seeds(self, n: int = -1, xinit: np.uint64 = None) \
            -> Tuple[np.ndarray, np.ndarray]:
        '''
        Generates the requested number of seeds. Each generator requires two
        initializer/seeds (x, a), one 64-bit unsigned integer (x) and one
        32-bit unsigned integer (a). The seed x is updated on each call to
        the random number generator, while seed a remains unchanged.

        Parameters
        ----------
        n: int
            The requested number of seeds - must not exceed the value of
            property :py:attr:`~Random.maxseeds`. Use -1 to
            generate all available seeds.
        xinit: int
            A 64-bit unsigned integer random number as initializer. If
            none is provided the builtin numpy random generator is used
            to create one.

        Returns
        -------
        x: np.uint64 vector
            A vector of n 64-bit unsigned integers representing x seeds.
        a: np.uint32 vector
            A vector of n 32-bit unsigned integers representing a seeds.

        Note
        ----
        Each random number generator should be initialized by a pair
        odf seeds x[i], a[i].
        '''

        n = int(n)
        if n < 0:
            n = self.maxseeds

        if n > self._primes.shape[0] - 1:
            raise ValueError('Cannot generate more than {} seeds ... '\
            'requested {} seeds!'.format(self._primes.shape[0] - 1, n))

        canReinitialize = False
        if xinit is None:
            # tmp = np.uint64(np.random.rand(2)*0xFFFFFFFF)
            # xinit = tmp[0] + tmp[1]*np.uint64(0xFFFFFFFF)
            xinit = self.initializer()
            canReinitialize = True
        else:
            xinit = np.uint64(xinit)

        res = 1
        while res != 0:
            x = np.zeros([n], dtype=np.uint64)
            a = np.zeros([n], dtype=np.uint32)
            if self._init_RNG is not None:
                res = self._init_RNG(
                    x.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
                    a.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    self._fora.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    ctypes.c_uint32(n), ctypes.c_uint64(xinit)
                )
            else:
                res = _rng_init(x, a, self._fora, n, xinit)

            if res != 0:
                if canReinitialize:
                    # the autogenerated initializer is bad
                    # try with a new initializer
                    xinit = self.initializer()
                else:
                    # user provided initializer is invalid ... stop
                    raise ValueError(
                        'The provided initializer {} is not valid. '
                        'Try a different 64-bit unsigned integer!'.format(
                            xinit))

        return x, a

    def cc(self, n: int, x: int = None, a: int = None) \
            -> Tuple[np.ndarray, int, int]:
        '''
        Generate n random numbers from [0.0, 1.0] using seeds x and a.
        If provided, the two seeds should be generated by the
        :py:meth:`~Random.seeds` method function.

        Parameters
        ----------
        n: int
            Generates n random number.
        x: int
            64-bit unsigned integer seed updated and returned on
            completion.
        a: int
            32-bit unsigned integer seed updated and returned on
            completion.

        Returns
        -------
        data: np.float32
            A vector of n random numbers.
        xnew: int
            Updated value of the input seed x.
        a: int
            The input seed a.
        '''
        if x is None or a is None:
            x, a = self.seeds(1)
            x, a = x[0], a[0]

        n = max(0, int(n))
        data = np.zeros([n], dtype=np.float32)
        if self._rand_cc is not None:
            cx = ctypes.c_uint64(x)
            self._rand_cc(
                data.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_uint32(n),
                ctypes.byref(cx),
                ctypes.c_uint32(a)
            )
            x = cx.value
        else:
            data, x, a = _rand_cc(data, n, x, a)

        return data, x, a

    def co(self, n, x: int = None, a: int = None) \
            -> Tuple[np.ndarray, int, int]:
        '''
        Generate n random numbers from [0.0, 1.0) using seeds x and a.
        If provided, the two seeds should be generated by the
        :py:meth:`~Random.seeds` method.

        Parameters
        ----------
        n: int
            Generates n random number.
        x: int
            64-bit unsigned integer seed updated and returned on
            completion.
        a: int
            32-bit unsigned integer seed updated and returned on
            completion.

        Returns
        -------
        data: np.float32
            A vector of n random numbers.
        xnew: int
            Updated value of the input seed x.
        a: int
            The input seed a.
        '''
        if x is None or a is None:
            x, a = self.seeds(1)
            x, a = x[0], a[0]

        n = max(0, int(n))
        data = np.zeros([n], dtype=np.float32)
        if self._rand_co is not None:
            cx = ctypes.c_uint64(x)
            self._rand_co(
                data.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_uint32(n),
                ctypes.byref(cx),
                ctypes.c_uint32(a)
            )
            x = cx.value
        else:
            data, x, a = _rand_co(data, n, x, a)

        return data, x, a

    def oc(self, n, x: int = None, a: int = None) \
            -> Tuple[np.ndarray, int, int]:
        '''
        Generate n random numbers from (0.0, 1.0] using seeds x and a.
        If provided, the two seeds shoud be generated by the
        :py:meth:`~Random.seeds` method.

        Parameters
        ----------
        n: int
            Generates n random number.
        x: int
            64-bit unsigned integer seed updated and returned on
            completion.
        a: int
            32-bit unsigned integer seed updated and returned on
            completion.

        Returns
        -------
        data: np.float32
            A vector of n random numbers.
        xnew: int
            Updated value of the input seed x.
        a: int
            The input seed a.
        '''

        if x is None or a is None:
            x, a = self.seeds(1)
            x, a = x[0], a[0]

        n = max(0, int(n))
        data = np.zeros([n], dtype=np.float32)
        if self._rand_cc is not None:
            cx = ctypes.c_uint64(x)
            self._rand_oc(
                data.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_uint32(n),
                ctypes.byref(cx),
                ctypes.c_uint32(a)
            )
            x = cx.value
        else:
            data, x, a = _rand_cc(data, n, x, a)

        return data, x, a

    def oo(self, n, x: int = None, a: int = None) \
            -> Tuple[np.ndarray, int, int]:
        '''
        Generate n random numbers from (0.0, 1.0) using seeds x and a.
        If provided, the two seeds should be generated by the
        :py:meth:`~Random.seeds` method.

        Parameters
        ----------
        n: int
            Generates n random number.
        x: int
            64-bit unsigned integer seed updated and returned on
            completion.
        a: int
            32-bit unsigned integer seed updated and returned on
            completion.

        Returns
        -------
        data: np.float32
            A vector of n random numbers.
        xnew: int
            Updated value of the input seed x.
        a: int
            The input seed a.
        '''

        if x is None or a is None:
            x, a = self.seeds(1)
            x, a = x[0], a[0]

        n = max(0, int(n))
        data = np.zeros([n], dtype=np.float32)
        if self._rand_oo is not None:
            cx = ctypes.c_uint64(x)
            self._rand_oo(
                data.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                ctypes.c_uint32(n),
                ctypes.byref(cx),
                ctypes.c_uint32(a)
            )
            x = cx.value
        else:
            data, x, a = _rand_oo(data, n, x, a)

        return data, x, a

if __name__ == '__main__':
    import argparse
    import time

    if not os.path.isdir(USER_PRIMES_PATH):
        os.makedirs(USER_PRIMES_PATH)

    default_out = os.path.join(USER_PRIMES_PATH, 'safeprimes_base32.npz')
    n_default = 150000

    parser = argparse.ArgumentParser(description='OpenCL random number generator')
    parser.add_argument('command', metavar='COMMAND',
                        choices=['test', 'primes'],
                        help='Use "primes" to generate primes or "test" to '
                        'test the performance/initialization of the generator.')
    parser.add_argument('-n', '--nprimes', metavar='NUM_PRIMES',
                        type=int, default=n_default,
                        help='Generate the specified number of primes. Use -1 '
                        'for default number of primes ({}).'.format(n_default))
    parser.add_argument('-o', '--output', metavar='OUTPUT_FILE',
                        type=str, default="",
                        help='Output npz file name for the generated seeds! '
                             'Defaults to "{}".'.format(default_out))
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        required=False, help='Turn on verbose mode.')

    parser.add_argument('-p', '--python', action='store_true', default=False,
                        required=False,
                        help='Force python implementation of the random number '
                             'generator initialization!')

    args = parser.parse_args()

    command = args.command
    n = args.nprimes
    verbose = args.verbose
    out = args.output
    forcepy = bool(args.python)

    if n < 0:
        n = n_default
    
    if not out:
        out = default_out

    tstart = time.perf_counter()

    if command == 'primes':
        seeds = generate_primes(n, verbose=verbose)
        np.savez_compressed(out, data=np.array(seeds[:,:2]))
        print('Generated {} primes in {:.3f} s'.format(
            n, time.perf_counter() - tstart))
        print('Primes saved to "{:s}".'.format(out))

    elif command == 'test':
        mcrng = Random(forcepy=forcepy)
        X, A = mcrng.seeds(
            mcrng.maxseeds, xinit=np.uint64(12343546457474))    

        print('Generated {} seeds in {:.3f} ms'.format(
            X.size, 1000.0*(time.perf_counter() - tstart)))


    '''
    ngens = mcrnd.maxseeds
    nvals = 0

    X, A = mcrnd.seeds(ngens, xinit=np.uint64(12343546457474))
    t2 = time.perf_counter()
    print('Generated {} seeds in {:.3f} ms'.format(ngens, 1000.0*(t2-t1)))


    if nvals > 0:
        R = np.zeros([ngens, nvals], dtype=np.float32)
        for i in range(ngens):
            data, xnew, a = mcrnd.oo(nvals, X[i], A[i])
            R[i] = data
    '''
