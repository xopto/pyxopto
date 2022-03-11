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

import os
import os.path
import copy

import numpy as np

from xopto.mcvox import mc
from xopto.mcml.mcutil.fiber import MultimodeFiber
from xopto.cl import clinfo

DEFAULT_CONFIG = {
    'rmax': 15.0e-3,
    'batch_packets': 10000000,
    'num_packets': 1000000000,
    'medium':{
        'mua': 0.0,
        'mus': 0.0,
        'n': 1.5,
        'g': 0.9,
    },
    'sample':{
        'mua': 0.5e2,
        'mus': 350.0e2,
        'n': 1.5,
        'g': 0.9,
    },
    'inclusion':{
        'z': 0.5e-3,
        'd': 0.2e-3,
        'mua': 230.0e2,
        'mus': 100.0e2,
        'n': 1.0,
        'g': 0.9,
    },
    'volume':{
        'size': {'x': 1e-3, 'y': 1e-3, 'z': 1e-3},
        'n': {'x': 201, 'y': 201, 'z': 201}
    },
    'time':{
        'start': 0.0,
        'stop': 6e-12,
        'n': 60
    }
}

def default_configuration() -> dict:
    '''
    Returns a copy of the default simulation configuration that can be
    optionally customized.
    '''
    return copy.deepcopy(DEFAULT_CONFIG)

def process_config(config: dict) -> dict:
    '''
    Process the simulation configuration. This will create objects for
    the Monte Carlo simulator and return the amendded configuration.

    Parameters
    ----------
    config: dict
        Customized default configuration.

    Returns
    -------
    config: dict
        Configuration with the created Monte Carlo simulator objects.
    '''
    medium = config['medium']
    inclusion = config['inclusion']
    sample = config['sample']

    materials = mc.mcmaterial.Materials([
        mc.mcmaterial.Material(n=medium['n'], mua=medium['mua'],
                               mus=medium['mus'], pf=mc.mcpf.Hg(medium['g'])),
        mc.mcmaterial.Material(n=sample['n'], mua=sample['mua'],
                               mus=sample['mus'], pf=mc.mcpf.Hg(sample['g'])),
        mc.mcmaterial.Material(n=inclusion['n'], mua=inclusion['mua'],
                               mus=inclusion['mus'], pf=mc.mcpf.Hg(inclusion['g'])),
    ])

    vol = config['volume']
    dx, nx = vol['size']['x'], vol['n']['x']
    dy, ny = vol['size']['y'], vol['n']['y']
    dz, nz = vol['size']['z'], vol['n']['z']
    t = config['time']
    xaxis = mc.mcgeometry.Axis(-dx*0.5, dx*0.5, nx)
    yaxis = mc.mcgeometry.Axis(-dy*0.5, dy*0.5, ny)
    zaxis = mc.mcgeometry.Axis(0.0, dz, nz)
    taxis = mc.mcfluence.Axis(t['start'], t['stop'], t['n'])

    voxels = mc.mcgeometry.Voxels(xaxis, yaxis, zaxis)
    Z, Y, X = voxels.meshgrid()
    inclusion_mask = X**2 + (Z - inclusion['z'])**2 <= inclusion['d']**2/4.0

    source = mc.mcsource.Line()

    fluence = mc.mcfluence.Fluencet(xaxis, yaxis, zaxis, taxis)

    config.update({'materials': materials, 'voxels': voxels, 'source': source,
                   'fluence': fluence, 'inclusion_mask': inclusion_mask})

    return config

def tr_deposition(filename: str, config: dict = None,
                  cl_device: clinfo.cl.Device = None,
                  overwrite: bool = False, verbose: bool = False):
    '''
    Generates data for visualization of time-resolved energy deposition
    in medium with a cylindrical inclusion.

    Parameters
    ----------
    filename: str
        Output file for the simulation results or None.
    config: dict
        Monte Carlo configuration as returned by
        :py:func:`default_configuration` and optionally customized.
    cl_device:
        OpenCL device on which to run the simulations or None.
    overwrite: bool
        If False, existing output files will not be overwritten.
    verbose:
        If nonzero, progress reports are periodically printed to the console.

    Returns
    -------
    results: dict
        Simulation results and simulation configuration data.
    '''
    if filename and os.path.isfile(filename) and not overwrite:
        raise FileExistsError('The output file "{}" already exists!')

    if cl_device is None:
        cl_device = clinfo.gpu()

    if config is None:
        config = default_configuration()

    config = process_config(config)

    mc_obj = mc.Mc(voxels=config['voxels'], materials=config['materials'],
                   source=config['source'], fluence=config['fluence'],
                   cl_devices=cl_device,
                   cl_build_options=[mc.cloptions.FastMath])
    mc_obj.rmax = config['rmax']

    mc_obj.voxels.material[:] = 1
    mc_obj.voxels.material[config['inclusion_mask']] = 2

    results = None
    packets_launched = 0
    num_packets = config['num_packets']
    batch_packets = config['batch_packets']
    while packets_launched < num_packets:
        n_launch = min(batch_packets, num_packets - packets_launched)
        results = mc_obj.run(n_launch, out=results)
        packets_launched += n_launch
        if verbose:
            print('    Launched {:,d}/{:,d} packets'.format(
                packets_launched, num_packets), end='\r')
    
    _, fluence_res, _ = results
    result = {
        'materials': mc_obj.materials.todict(),
        'voxels': mc_obj.source.todict(),
        'fluence': fluence_res.todict(),
        'fluence_data': fluence_res.data,
        'rmax': mc_obj.rmax,
        'batch_packets': config['batch_packets'],
        'num_packets': packets_launched,
    }

    if filename:
        np.savez_compressed(filename, **result)

    return result

if __name__ == '__main__':
    from xopto.util.animation.tr import common

    common.main(tr_deposition, common.create_deposition_animation_xz,
                default_configuration())
