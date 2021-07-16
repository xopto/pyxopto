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

'''
This module provides functions for working with OpenCL resources.

- :py:func:`info` - lists the available OpenCL platforms and devices
- :py:func:`all_devices` - lists all the available OpenCL devices
- :py:func:`gpus` - lists the available OpenCL GPU devices
- :py:func:`cpus` - lists the available OpenCL CPU devices
- :py:func:`device_info` - returns a dictionary with information on the specified device
- :py:func:`device` - find a device that matches the given search criteria
'''

import pyopencl as cl

from typing import List


def info(show: bool = True, devices: list = None) -> str:
    '''
    Returns a string descriptor of all installed OpenCl devices.

    Parameters
    ----------
    show: bool
        If set to True, the information is displayed on screen.
    devices: list[cl.Device]
        List of OpenCL devices or None. If None (default), all OpenCL
        devices are listed.

    Returns
    -------
    descriptor: str
        If parameter show is False, a string descriptor is returned,
        otherwise the information is displayed on screen.

    Examples
    --------
    Dispalys information of all available OpenCL devices.

    >>> import clinfo
    >>> clinfo.info()
    '''
    info_str = ''
    for platform in cl.get_platforms():
        for device in platform.get_devices():
            if devices is not None and device not in devices:
                continue
            info_str += '\n===========================' \
                '====================================\n'
            info_str += 'Platform name: {}\n'.format(platform.name)
            info_str += 'Platform profile: {}\n'.format(platform.profile)
            info_str += 'Platform vendor: {}\n'.format(platform.vendor)
            info_str += 'Platform version: {}\n'.format(platform.version)
            info_str += '----------------------------' \
                '-----------------------------------\n'
            info_str += '\tDevice name: {}\n'.format(device.name)
            info_str += '\tDevice type: {}\n'.format(
                cl.device_type.to_string(device.type))
            info_str += '\tDevice available: {}\n'.format(
                {False:'No', True:'Yes'}[bool(device.available)])
            info_str += '\tDevice maximum clock frequency: {} MHz\n'.format(
                device.max_clock_frequency)
            info_str += '\tDevice default address space: {} bit\n'.format(
                device.address_bits)
            info_str += '\tDevice maximum number of constant kernel '\
                'arguments: {}\n'.format(device.max_constant_args)
            info_str += '\tDevice maximum kernel argument size: '\
                '{} B\n'.format(device.max_parameter_size)
            info_str += '\tDevice global memory: {} MB\n'.format(
                device.global_mem_size//1024//1024)
            info_str += '\tDevice global memory cache: {} kB\n'.format(
                device.global_mem_cache_size//1024)
            info_str += '\tDevice maximum global memory allocation size: '\
                '{} MB\n'.format(device.max_mem_alloc_size//1024//1024)
            info_str += '\tDevice local memory: {} kB\n'.format(
                device.local_mem_size//1024)
            info_str += '\tDevice constant memory: {} kB\n'.format(
                device.max_constant_buffer_size/1024)
            info_str += '\tDevice max clock speed: {} MHz\n'.format(
                device.max_clock_frequency)
            info_str += '\tDevice compute units: {}\n'.format(
                device.max_compute_units)
            info_str += '\tDevice maximum work group size: {}\n'.format(
                device.max_work_group_size)
            info_str += '\tMaximum number of work-item dimensions '\
                '(global or local): {}\n'.format(
                    device.max_work_item_dimensions)
            info_str += '\tMaximum number of work-items in each dimension '\
                'of the work group: {}\n'.format(
                    device.max_work_item_sizes)
            info_str += '\tMaximum number of image-objects that can be '\
                'passed to a kernel: {}\n'.format(
                    device.max_write_image_args)
            info_str += '\tDevice profiling timer resolution: {} ns\n'.format(
                device.profiling_timer_resolution)

            info_str += '\tDevice supported OpenCL extensions:\n'
            extensions = device.extensions.strip().split(' ')
            for extension in extensions:
                info_str += '\t\t {}\n'.format(extension)

    if show:
        print(info_str)
    else:
        return info_str

def _find_device(device_type: str = None, platform: str = None,
                 device: str = None) -> list:
    devices = []
    if device_type in ('gpu', 'cpu', 'accelerator', 'default'):
        device_type = {'gpu':cl.device_type.GPU, 'cpu':cl.device_type.CPU,
                       'accelerator':cl.device_type.ACCELERATOR}.get(
                           device_type, cl.device_type.DEFAULT)

    if device is not None:
        device = str(device).lower()
    if platform is not None:
        platform = str(platform).lower()

    for platform_item in cl.get_platforms():
        if platform is not None and platform not in platform_item.name.lower():
            continue
        for device_item in platform_item.get_devices():
            if device is not None and device not in device_item.name.lower():
                continue
            if device_type is not None and device_item.type != device_type:
                continue
            devices.append(device_item)
    return devices

def all_devices(platform: str = None, device: str = None) -> List[cl.Device]:
    '''
    Returns a list of all available OpenCl devices (a list of cl.Device).
    Devices can be filtered by the platform and device name. The matching of
    names is case insensitive.

    Parameters
    ----------
    platform: str
        If not None, the platform name of the returned devices must include
        this string (matching is case insensitive).
    device: str
        If not None, the device name of the returned devices must include
        this string (matching is case insensitive).

    Returns
    -------
    devices: List[cl.Device]
        List of found OpenCL devices.
    '''
    return _find_device(None, platform, device)

def gpus(platform: str = None, device: str = None) -> List[cl.Device]:
    '''
    Returns a list of all available OpenCl GPU devices
    (a list of cl.Device). Devices can be filtered by the platform and
    device name. The matching of names is case insensitive.

    Parameters
    ----------
    platform: str
        If not None, the platform name of the returned devices must include
        this string (matching is case insensitive).
    device: str
        If not None, the device name of the returned devices must include
        this string (matching is case insensitive).

    Returns
    -------
    devices: List[cl.Device]
        List of found OpenCL GPU devices.
    '''
    return _find_device(cl.device_type.GPU, platform, device)

def gpu(platform=None, device=None, index=0) -> cl.Device:
    '''
    Returns an available OpenCl GPU device at the specified index.
    Devices can be filtered by the platform and
    device name. The matching of names is case insensitive.

    Parameters
    ----------
    platform: str
        If not None, the platform name of the returned devices must include
        this string (matching is case insensitive).
    device: str
        If not None, the device name of the returned devices must include
        this string (matching is case insensitive).
    index: int
        Sets the number of matching OpenCL devices that will be skipped during
        the search. Useful for hosts with multiple equal OpenCL devices.
        Default value is 0 (return the first matching OpenCL device).

    Returns
    -------
    device: cl.Device
        The first founfd OpenCL GPU device.
    '''
    return _find_device(cl.device_type.GPU, platform, device)[index]

def cpus(platform=None, device=None) -> List[cl.Device]:
    '''
    Returns a list of all available OpenCL CPU devices
    (a list of pyopencl.Device). Devices can be filtered by the platform and
    device name. The matching of names is case insensitive.

    Parameters
    ----------
    platform: str
        If not None, the platform name of the returned devices must include
        this string (matching is case insensitive).
    device: str
        If not None, the device name of the returned devices must include
        this string (matching is case insensitive).

    Returns
    -------
    devices: List[cl.Device]
        List of found OpenCL CPU devices.
    '''
    return _find_device(cl.device_type.CPU, platform, device)

def cpu(platform: str = None, device: str = None, index: int = 0) -> cl.Device:
    '''
    Returns an available OpenCL CPU device at the specified index.
    Devices can be filtered by the platform and
    device name. The matching of names is case insensitive.

    Parameters
    ----------
    platform: str
        If not None, the platform name of the returned devices must include
        this string (matching is case insensitive).
    device: str
        If not None, the device name of the returned devices must include
        this string (matching is case insensitive).
    index: int
        Sets the number of matching OpenCL devices that will be skipped during
        the search. Useful for hosts with multiple equal OpenCL devices.
        Default value is 0 (return the first matching OpenCL device).

    Returns
    -------
    device: cp.Device
        The first founfd OpenCL CPU device.
    '''
    return _find_device(cl.device_type.CPU, platform, device)[index]

def device_info(device: cl.Device) -> dict:
    '''
    Returns a dict with information on the given OpenCL device. the returned
    dict will have the following fields:

    * name: str
    * vendor: str
    * available: bool
    * extensions: List[str],
    * addressSpace: int,
    * platform: str,
    * maxConstKernelArgs: int,
    * maxKernelArgSize: int,
    * type: str,
    * globalMemory: int,
    * globalMemoryCacheSize: int,
    * maxAllocationSize: int,
    * localMemory: int,
    * constantMemory: int,
    * maxClkFrequency: int,
    * computeUnits: int,
    * maxWorkGroupSize: int

    Parameters
    ----------
    device: cl.Device
        An OpenCL device. Use the all_devices, gpus or cpus functions to get
        a list of all available OpenCl or just GPU or CPU devices, respectively.

    Returns
    -------
    info: dict
        A dictionary with device information.
    '''
    return {
        'name':device.name,
        'vendor':device.vendor,
        'available':device.available,
        'extensions':device.extensions.strip().split(' '),
        'addressSpace':device.address_bits,
        'platform':device.platform,
        'maxConstKernelArgs':device.max_constant_args,
        'maxKernelArgSize':device.max_parameter_size,
        'type':cl.device_type.to_string(device.type),
        'globalMemory':device.global_mem_size,
        'globalMemoryCacheSize':device.global_mem_cache_size,
        'maxAllocationSize':device.max_mem_alloc_size,
        'localMemory':device.local_mem_size,
        'constantMemory':device.max_constant_buffer_size,
        'maxClkFrequency':device.max_clock_frequency,
        'computeUnits':device.max_compute_units,
        'maxWorkGroupSize':device.max_work_group_size
    }

def device(what: str or List[str] or cl.Device = None, index: int = 0) \
        -> cl.Device:
    '''
    Finds an OpenCl device according to the given keyword(s).
    The keywords are compared to the OpenCL device name, device vendor name and
    device platform name of each available OpenCL device.

    Parameters
    ----------
    device: str or list[str] or cl.Device
        A keyword or a list of keywords that are compared to the OpenCL
        device name, device vendor name and device platform name
        (comparison is case insensitive).
        The first OpenCL device that contains the keyword is returned if
        the index argument is set to 0.
        If this parameter is an instance of :class:`pyopencl.Device`,
        the instance is immediately returned without any further action.
    index: int
        Sets the number of matching OpenCL devices that will be skipped during
        the search. Useful for hosts with multiple equal OpenCL devices.
        Default value is 0 (return the first matching OpenCL device).
    what: str or a list/tuple of str

    Returns
    -------
    dev: cl.Device
        An OpenCl device if one that matches the keywords is found.

    Examples
    --------
    Searches through all available OpenCL devices for the "amd" keyword. If
    none is found, the next keyword "nvidia" is used, and finally the "cpu"
    keyword is used if the first two keywords are not found in any of the
    available OpenCL device name, vendor name or platform name.

    >>> import clinfo
    >>> dev = clinfo.device(["amd", "nvidia", "cpu"])
    '''
    if what is None:
        what = ''
    if isinstance(what, str):
        what = [what]
    elif isinstance(what, cl.Device):
        return what

    clgpus = gpus()
    clcpus = cpus()
    for item in what:
        item = str(item).lower()
        for gpu in clgpus:
            if str(gpu.vendor).lower().find(item) >= 0 or \
                    str(gpu.name).lower().find(item) >= 0 or \
                    str(gpu.platform).lower().find(item) >= 0:
                if index <= 0:
                    return gpu
                else:
                    index -= 1

        for cpu in clcpus:
            if str(cpu.vendor).lower().find(item) >= 0 or \
                    str(cpu.name).lower().find(item) >= 0 or \
                    str(cpu.platform).lower().find(item) >= 0:
                if index <= 0:
                    return cpu
                else:
                    index -= 1

def device_extension(device: cl.Device, extension: str or list or tuple) \
        -> bool or List[bool]:
    '''
    Checks if the OpenCL device supports the given extension. See
    https://www.khronos.org/registry/OpenCL/sdk/1.2/docs/man/xhtml/EXTENSION.html
    for a list of extensions supported by the OpenCl 1.2 language version.
    Some frequently used extensions:

        - cl_khr_fp64
            Double precision floating point.

        - cl_khr_fp16
            Half precision floating point.

        - cl_khr_int64_base_atomics
            64-bit integer atomic operations.

        - cl_khr_int64_extended_atomics
            64-bit integer extended atomic operations that include "min", "max"
            "and", "or" and "xor" operations.

    Parameters
    ----------
    device: cl.Device
        OpenCl device instance.
    extension: str or list/tuple of str
        OpenCL extensions.

    Returns
    -------
    available: bool or list of bool
        True, if the extension is supported by the device, else False.
    '''
    if isinstance(extension, str):
        return extension in device.extensions
    else:
        return [item in device.extensions for item in extension]

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        'Find and list capabilities of OpenCL devices')

    parser.add_argument(
        '-g', '--gpus', action='store_true',
        required=False, help='Search for GPUs only.')

    parser.add_argument(
        '-c', '--cpus', action='store_true',
        required=False, help='Search for CPUs only.')

    args = parser.parse_args()

    if args.gpus:
        print(info(False, devices=gpus()))
    elif args.cpus:
        print(info(False, devices=cpus()))
    else:
        print(info(False))
