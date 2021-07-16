.. ****************************** Begin license ********************************
.. Copyright (C) Laboratory of Imaging technologies,
..               Faculty of Electrical Engineering,
..               University of Ljubljana.
..
.. This file is part of PyXOpto.
..
.. PyXOpto is free software: you can redistribute it and/or modify
.. it under the terms of the GNU General Public License as published by
.. the Free Software Foundation, either version 3 of the License, or
.. (at your option) any later version.
..
.. PyXOpto is distributed in the hope that it will be useful,
.. but WITHOUT ANY WARRANTY; without even the implied warranty of
.. MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
.. GNU General Public License for more details.
..
.. You should have received a copy of the GNU General Public License
.. along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
.. ******************************* End license *********************************

.. _{{ MC }}-working-with-opencl-devices-label:

.. include:: common.rst

.. _opencl-devices-label:

OpenCL devices
==============

OpenCL devices can be easily discovered and selected with the help of the
:py:mod:`xopto.cl.clinfo` module.

To list all available OpenCL devices use the
:py:func:`~xopto.cl.clinfo.info` function. The human-readable output contains
information on the relevant properties of the installed OpenCL devices.

.. code-block:: python

    from xopto.cl import clinfo

    clinfo.info()

An example output of the above command on a system with a single integrated
Intel GPU:

.. code-block:: bash

    ===============================================================
    Platform name: Intel(R) OpenCL HD Graphics
    Platform profile: FULL_PROFILE
    Platform vendor: Intel(R) Corporation
    Platform version: OpenCL 2.1 
    ---------------------------------------------------------------
            Device name: Intel(R) Gen9 HD Graphics NEO
            Device type: ALL | GPU
            Device available: Yes
            Device maximum clock frequency: 1150 MHz
            Device default address space: 64 bit
            Device maximum number of constant kernel arguments: 8
            Device maximum kernel argument size: 1024 B
            Device global memory: 12559 MB
            Device global memory cache: 512 kB
            Device maximum global memory allocation size: 4095 MB
            Device local memory: 64 kB
            Device constant memory: 4194296.0 kB
            Device max clock speed: 1150 MHz
            Device compute units: 24
            Device maximum work group size: 256
            Maximum number of work-item dimensions (global or local): 3
            Maximum number of work-items in each dimension of the work group: [256, 256, 256]
            Maximum number of image-objects that can be passed to a kernel: 128
            Device profiling timer resolution: 83 ns
            Device supported OpenCL extensions:
                    cl_khr_byte_addressable_store
                    cl_khr_fp16
                    cl_khr_global_int32_base_atomics
                    cl_khr_global_int32_extended_atomics
                    cl_khr_icd
                    cl_khr_local_int32_base_atomics
                    cl_khr_local_int32_extended_atomics
                    cl_intel_subgroups
                    cl_intel_required_subgroup_size
                    cl_intel_subgroups_short
                    cl_khr_spir
                    cl_intel_accelerator
                    cl_intel_driver_diagnostics
                    cl_khr_priority_hints
                    cl_khr_throttle_hints
                    cl_khr_create_command_queue
                    cl_intel_subgroups_char
                    cl_intel_subgroups_long
                    cl_khr_fp64
                    cl_khr_subgroups
                    cl_khr_il_program
                    cl_intel_spirv_device_side_avc_motion_estimation
                    cl_intel_spirv_media_block_io
                    cl_intel_spirv_subgroups
                    cl_khr_spirv_no_integer_wrap_decoration
                    cl_intel_unified_shared_memory_preview
                    cl_khr_mipmap_image
                    cl_khr_mipmap_image_writes
                    cl_intel_planar_yuv
                    cl_intel_packed_yuv
                    cl_intel_motion_estimation
                    cl_intel_device_side_avc_motion_estimation
                    cl_intel_advanced_motion_estimation
                    cl_khr_int64_base_atomics
                    cl_khr_int64_extended_atomics
                    cl_khr_image2d_from_buffer
                    cl_khr_depth_images
                    cl_intel_media_block_io
                    cl_khr_3d_image_writes
                    cl_intel_va_api_media_sharing


To get a list of all available OpenCL GPU devices use the
:py:func:`~xopto.cl.clinfo.gpus` function. Two optional input arguments
allow filtering of the  OpenCL devices by a platform name
(:code:`platform`) and / or device name (:code:`device`). The matching of the
given device and / or platform name is done in a case insensitive way. The
target OpenCL device and / or platform name must contain the given string.

In the following example we fetch a list of installed OpenCL GPU devices:

.. code-block:: python

    from xopto.cl import clinfo

    gpus = clinfo.gpus()

In this example, we fetch only the Nvidia GPU devices:

.. code-block:: python

    from xopto.cl import clinfo

    nvidia_gpus = clinfo.gpus(platform='nvidia')

Note that the :py:func:`~xopto.cl.clinfo.gpus` returns :python:`None` if no
devices are found that match the search criteria.

We can use the :py:func:`xopto.cl.clinfo.info` function to print out
human-readable information on the found OpenCL devices.

.. code-block:: python

    from xopto.cl import clinfo

    nvidia_gpus = clinfo.gpus(platform='nvidia')

    clinfo.info(devices=nvidia_gpus)

Alternatively, the :py:func:`~xopto.cl.clinfo.gpu` function can be used to
retrieve the n-th OpenCL device that matches the search criteria. The following
example returns the second installed Nvidia GPU device. Note that the parameter
:code:`index` that is used to select a particular device is zero-based
(0 - 1 :superscript:`st` device, 1 - 2 :superscript:`nd` device, ...).

.. code-block:: python

    from xopto.cl import clinfo

    nvidia_gpus = clinfo.gpu(platform='nvidia', index=1)

Use the :py:func:`~xopto.cl.clinfo.cpus` and :py:func:`~xopto.cl.clinfo.cpu`
functions to find or select CPU devices. The two functions take the same
parameters as the :py:func:`~xopto.cl.clinfo.gpus` and
:py:func:`~xopto.cl.clinfo.gpu` functions.

To retrieve OpenCL device information into a :python:`dict` instead of
displaying the information in the console, use the 
:py:func:`~xopto.cl.clinfo.device_info` function.

.. code-block:: python

    from xopto.cl import clinfo

    nvidia_gpu = clinfo.gpu(platform='nvidia')
    gpu_info = clinfo.device_info(nvidia_gpu)

    print(gpu_info['name'])
    print(gpu_info['platform'])

To check if a device supports one or more OpenCL extensions, use the
:py:func:`~xopto.cl.clinfo.device_extension` function that returns True if
the extension is supported and False if not (returns a list of boolean values
if multiple extensions are queried). See the 
`full list <https://www.khronos.org/registry/OpenCL/sdk/1.2/docs/man/xhtml/EXTENSION.html>`_
of extensions supported by the OpenCL 1.2 programming language.

.. code-block:: python

    from xopto.cl import clinfo

    nvidia_gpu = clinfo.gpu(platform='nvidia')

    supports_double = clinfo.device_extension(nvidia_gpu, 'cl_khr_fp64')

    supports_64bit_atomic_operations = clinfo.device_extension(
        nvidia_gpu, 'cl_khr_int64_base_atomics')

