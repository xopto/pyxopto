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

from xopto.mcml import mc
from xopto.mcml.mcutil.fiber import MultimodeFiber
from xopto.cl import clinfo


DEFAULT_CONFIG = {
    'rmax': 15.0e-3,
    'batch_packets': 100000,
    'sv_num_packets': 1000,
    'sample':{
        'mua': 2.0e2,
        'mus': 125.0e2,
        'n': 1.33,
        'g': 0.85,
        'd': np.inf
    },
    'fiber':{
        'ncore': 1.452,
        'dcore': 200e-6,
        'dcladding': 220e-6,
        'na': 0.22
    },
    'fiber_sds': 0.5e-3,
    'fiber_incidence': np.deg2rad(np.linspace(-45, 60, 22)), # 5 deg resolution
    'trace_maxlen': 512,
    'radial_detector': {'start': 0.0, 'stop': 1.0e-3, 'n': 100},
    'sv':{
        'xaxis': {'start': -0.75e-3, 'stop': 0.75e-3, 'n': 300},
        'yaxis': {'start': -0.75e-3, 'stop': 0.75e-3, 'n': 300},
        'zaxis': {'start':  0.0e-3,  'stop': 1.0e-3,  'n': 200}
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
    fp_eps = np.finfo(np.float32).eps
    sample = config['sample']
    fiber_cfg = config['fiber']
    fiber_cosmin = (1.0 - (fiber_cfg['na']/fiber_cfg['ncore'])**2)**0.5

    fiber = MultimodeFiber(**fiber_cfg)

    config['layers'] = mc.mclayer.Layers([
        mc.mclayer.Layer(d=np.inf, mua=0.0, mus=0.0,
                         n=fiber_cfg['ncore'], pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(d=sample['d'], mua=sample['mua'], mus=sample['mus'],
                         n=sample['n'], pf=mc.mcpf.Hg(sample['g'])),
        mc.mclayer.Layer(d=np.inf, mua=0.0, mus=0.0, n=1.0,
                         pf=mc.mcpf.Hg(0.0))
    ])

    config['source'] = mc.mcsource.UniformFiber(
        fiber,
        direction=(0.0, 0.0, 1.0),
        position=(-config['fiber_sds']*0.5, 0.0, 0.0)
    )

    config['detectors'] = mc.mcdetector.Detectors(
        top = mc.mcdetector.Radial(
            mc.mcdetector.Axis(**config['radial_detector'])
        ),
        bottom = mc.mcdetector.Radial(
            mc.mcdetector.Axis(**config['radial_detector'])
        )
    )

    config['trace'] = mc.mctrace.Trace(
        maxlen=config['trace_maxlen'],
        options=mc.mctrace.Trace.TRACE_ALL,
        filter=mc.mctrace.Filter(
            z=(-fp_eps, fp_eps),
            pz=(-1.0, -fiber_cosmin),
            r=(0.0, fiber_cfg['dcore']*0.5, (config['fiber_sds']*0.5, 0.0))
        )
    )
    config['sv'] = mc.mcsv.SamplingVolume(
        xaxis=mc.mcsv.Axis(**config['sv']['xaxis']), 
        yaxis=mc.mcsv.Axis(**config['sv']['yaxis']),
        zaxis=mc.mcsv.Axis(**config['sv']['zaxis'])
    )

    return config

def sv_fiber_incidence(filename: str, config: dict = None,
                       cl_device: clinfo.cl.Device = None,
                       overwrite: bool = False,
                       verbose: bool = False):
    '''
    Generates data for visualization of the sampling volume in reflectance
    configuration as a function of the source fiber incidence.

    Parameters
    ----------
    filename: str
        Output file for the simulation results or None.
    config: dict
        Monte Carlo configuration as returned by
        :py:func:`default_configuration` and optionally customized.
    cl_device:
        Opencl device on which to run the simulations or None.
    overwrite: bool
        If False, existing output files will not be overwritten.
    verbose:
        If nonzero, progress reports are periodically printed to the console.

    Returns
    -------
    results: dict
        Simulation results with and simulation configuration data.
    '''
    if filename and os.path.isfile(filename) and not overwrite:
        raise FileExistsError('The output file "{}" already exists!')

    if cl_device is None:
        cl_device = clinfo.gpu()

    if config is None:
        config = default_configuration()

    config = process_config(config)

    mc_obj = mc.Mc(layers=config['layers'], source=config['source'],
                   detectors=config['detectors'], trace=config['trace'],
                   cl_devices=cl_device)
    mc_obj.rmax = config['rmax']
    mc_obj.source.position = (-config['fiber_sds']*0.5, 0.0, 0.0)
    ref_filter = mc_obj.trace.filter
    mc_obj.trace.filter = mc.mctrace.Filter(
        z=ref_filter.z,
        pz=ref_filter.pz,
        r=(ref_filter.r[0][0], ref_filter.r[0][1], (config['fiber_sds']*0.5, 0))
    )

    N = len(config['fiber_incidence'])
    sv_data = []
    sv = config['sv']
    for index, incidence in enumerate(config['fiber_incidence']):
        if verbose:
            print('Processing incidence angle {:.1f}°; {}/{}'.format(
                np.rad2deg(incidence), index + 1, N))

        mc_obj.source.direction = [np.sin(incidence), 0.0, np.cos(incidence)]
        num_packets = 0
        sv.clear()
        while num_packets < config['sv_num_packets']:
            trace_res, _, _ = mc_obj.run(config['batch_packets'])
            sv = mc_obj.sampling_volume(trace_res, sv)
            num_packets += len(trace_res)
            if verbose:
                print('    Collected {}/{} packets'.format(
                    num_packets, config['sv_num_packets']), end='\r')
        sv_data.append(np.array(sv.data))

    result = {
        'trace': mc_obj.trace.todict(),
        'layers': mc_obj.layers.todict(),
        'source': mc_obj.source.todict(),
        'sv': sv.todict(),
        'rmax': mc_obj.rmax,
        'batch_packets': config['batch_packets'],
        'num_packets': num_packets,
        'fiber_incidence': config['fiber_incidence'],
        'sv_data': sv_data
    }

    if filename:
        np.savez_compressed(filename, **result)

    return result

if __name__ == '__main__':
    from . import common

    common.main(sv_fiber_incidence, common.create_sv_animation_xz,
                default_configuration())