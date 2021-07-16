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

import pyopencl as cl
from xopto.cl import clinfo

import numpy as np

def test_cl_khr_int64_base_atomics(devices: cl.Device) -> bool:
    for device in devices:
        if not clinfo.device_extension(device, 'cl_khr_int64_base_atomics'):
            return False

    np_n = np.array([0], dtype=np.uint64)

    try:
        src = '\n'.join((
            '#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable',
            '',
            '__kernel void test(__global ulong *n){',
            '	atom_inc(n);',
            '};',
        ))
        cl_ctx = cl.Context(devices)
        cl_queue = cl.CommandQueue(cl_ctx)
        cl_exec = cl.Program(cl_ctx, src).build()

        mf = cl.mem_flags

        cl_n = cl.Buffer(cl_ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=np_n)
        cl_exec.test(cl_queue, [1], [1], cl_n).wait()
        cl.enqueue_copy(cl_queue, np_n, cl_n).wait()
    except cl.RuntimeError:
        pass

    return np_n[0] == 1


if __name__ == '__main__':
    devices = clinfo.gpus()
    result = test_cl_khr_int64_base_atomics(devices)
    print('Devices ', devices, 'test cl_khr_int64_base_atomics:', result)
