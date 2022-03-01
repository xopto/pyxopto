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
import argparse
import os.path

from xopto import USER_TMP_PATH
from xopto.cl import clinfo
from xopto.mcvox.test import validate
from xopto.mcvox import mcoptions


if __name__ == '__main__':
    devices = ['amd', 'nvidia', 'cpu', 'hd']
    wgsize = None
    nphotons_standard = 10e6
    default_math = 'std'
    
    parser = argparse.ArgumentParser(description='MCVOX performance test')
    parser.add_argument('-d', '--device', metavar='DEVICE_NAME',
                        type=str, default="",
                        help='OpenCL device identifier/name')
    parser.add_argument('-i', '--index', metavar='DEVICE_INDEX',
                        type=int, default=0,
                        help='OpenCL device index')
    parser.add_argument('-w', '--wgsize', metavar='WG_SIZE',
                        type=int, default=-1,
                        help='Workgroup size ... use -1 for auto adjust')
    parser.add_argument('-n', '--nphotons', metavar='NUM_PHOTONS',
                        type=float, default=10e6,
                        help='Number of photon packets. Minimum is 1,000,000.')
    parser.add_argument('-l', '--lut', action='store_true', default=False,
                        required=False,
                        help='Use lookup table-based sampling of the '\
                             'scattering angle.')
    parser.add_argument('-m', '--math', metavar='OPENCL_MATH',
                        type=str, default=default_math,
                        choices=('native', 'std'),
                        help='OpenCL Math functions - "native" for native, '\
                             'or "std" (default) for standard Math functions.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        required=False,
                        help='Enable verbose mode.')
    parser.add_argument('--method', metavar='MC_METHOD', type=str,
                        choices = ('albedo_weight', 'aw',
                                   'albedo_rejection', 'ar',
                                   'microscopic_beer_lambert', 'mbl'),
                        default='albedo_weight', required=False,
                        help='Select the Monte Carlo simulation method.')

    mc_options = [
        mcoptions.McFloatLutMemory.constant_mem,
        mcoptions.McMaterialMemory.constant_mem
    ]

    args = parser.parse_args()
    if args.device:
        devices = [args.device]
    device_index = max(int(args.index), 0)
    if args.wgsize > 0:
        wgsize = args.wgsize
    nphotons = max(int(args.nphotons), 100000)
    if args.math == 'native':
        mc_options.append(mcoptions.McUseNativeMath.on)
    if args.method:
        kernel_method = str(args.method)
        kernel_method = {'ar': 'albedo_rejection',
                         'aw': 'albedo_weight',
                         'mbl': 'microscopic_beer_lambert'}.get(
                               kernel_method, kernel_method)
        mc_options.append(getattr(mcoptions.McMethod, kernel_method))

    usepflut = bool(args.lut)

    verbose = bool(args.verbose)

    options = {
        'nphotons': nphotons,
        'usepflut': usepflut,
        'wgsize': wgsize,
        'verbose': verbose,
        'visualize': False,
        'cl_devices': clinfo.device(devices, device_index),
        'cl_build_options': ['-cl-fast-relaxed-math', '-cl-mad-enable'],
        'options': mc_options,
        'exportsrc': os.path.join(USER_TMP_PATH, 'mcvox_performancetest.h')
    }
    spf_mode = 'LOOKUP TABLE-BASED' if usepflut else 'ANALYTICAL'
    print('--------------------------------------------------------------------------')
    print('Staring performance test with:\n'
          '   {:s} sampling of the scattering phase function (SPF)\n'
          '   {} MC kernel method'.format(spf_mode, kernel_method)
    )
    print('...')
    t1 = time.perf_counter()
    test_obj = validate.SingleLayerUniformFiberRadial(**options)
    t_mc = test_obj.run(**validate.run_options(options))
    t2 = time.perf_counter()
    print('All done')
    if nphotons >= nphotons_standard:
        res = 'passed' if test_obj.passed() else 'failed'
        print('Validation {}!'.format(res))

    k = float(nphotons_standard)/float(nphotons)
    dev = test_obj.mc.cl_device[0]

    print('--------------------------------------------------------------------------')
    print('Test completed using {} SPF and {} MC kernel method:'.format(spf_mode, kernel_method))
    print('    {:>6.1f} s - {} {}'.format(k*(t2 - t1), dev.vendor, dev.name))
    print('               ({:.1f} s for the Monte Carlo kernel)'.format(k*t_mc))
    print('--------------------------------------------------------------------------')
    print('Compare to other devices (all for ANALYTICAL SPF and REGULAR MC kernel):')
    print('--------------------------------------------------------------------------')
    print('       8.3 s - NVIDIA Corporation GeForce, RTX A6000\n'
          '               (Ubuntu 20.04, 460.73) GPU @ 1700 MHz')
    print('       9.9 s - NVIDIA Corporation GeForce, ASUS GF RTX 2080 TI STRIX\n'
          '               (Ubuntu 20.04, 460.80) GPU @ 1900 MHz')
    print('      21.7 s - NVIDIA Corporation GeForce, ASUS GF GTX 1080 TI STRIX\n'
          '               (Windows 8.1) GPU @ 1900 MHz')
    print('      38.6 s - NVIDIA Corporation GeForce, ASUS GF GTX 1070 OC\n'
          '               (Ubuntu 20.04, 460.84) GPU @ 1900 MHz')
    print('      48.0 s - Advanced Micro Devices, Inc. Grenada R9 390 \n'
          '               (Ubuntu 16.04 amdgpu-pro-17.50-511655 driver)')
    print('      52.3 s - Advanced Micro Devices, Inc. Polaris20 RX 580\n'
          '               (Win 8.1, 18.3.4 AMD driver) GPU @ 1420 MHZ, MEM @ 1750 MHZ')
    print('      54.0 s - Advanced Micro Devices, Inc. Tahiti R9 280\n'
          '               (Win 8.1, 17.7.1 driver) GPU OC @ 1200 MHz, MEM @ 1500 MHZ')
    print('      58.0 s - Advanced Micro Devices, Inc. Grenada R9 390\n'
          '               (Win 8.1 16.7.2 driver)')
    print('      60.0 s - Advanced Micro Devices, Inc. Polaris20 RX 570\n'
          '               (Ubuntu 16.04, amdgpu-pro-17.50 driver)')
    print('      63.4 s - Advanced Micro Devices, Inc. Tahiti R9 280\n'
          '               (Win 8.1, 17.7.1 driver) GPU @ 1020 MHz, MEM @ 1500 MHZ')
    print('      82.7 s - GTX 770 (Ubuntu 16.04, NVIDIA 390.25 driver)')
    print('     114.0 s - Tesla K40c Linux (Linux - cluster)')
    print('     178.0 s - Advanced Micro Devices, Inc. Grenada R9 390\n'
          '               (Win 8.1, 17.2.1 driver)')
    print('     269.3 s - NVIDIA Corporation GeForce, GT 750M\n'
          '               (Windows 8.1) GPU @ 1162 MHz')
    print('     528.4 s - Intel(R) Corporation Intel(R) HD Graphics 4600')
    print('     496.7 s - Intel(R) Corporation Intel(R) Gen9 HD Graphics NEO')
    print('    4679.6 s - Intel(R) Corporation Intel(R) Core(TM) i7-4700HQ CPU @ 2.40GHz')
