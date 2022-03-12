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
from typing import List

import jinja2

from xopto.dataset import DATASET_PATH, SV_TEMPLATE_PATH
from xopto.dataset.render import common

DEFAULT_WAVELENGTH = common.DEFAULT_WAVELENGTH

# Data common to all datasets.
RI_DIGITS = common.RI_DIGITS
RI_AIR = round(common.RI_AIR(DEFAULT_WAVELENGTH), RI_DIGITS)
RI_WATER = round(common.RI_WATER(DEFAULT_WAVELENGTH), RI_DIGITS) # 1.33
RI_FUSEDSILICA = round(common.RI_FUSED_SILICA(DEFAULT_WAVELENGTH), RI_DIGITS) # 1.452
PROBE_REFLECTIVITY = common.PROBE_REFLECTIVITY(DEFAULT_WAVELENGTH)
PROBE_DIAMETER = common.PROBE_DIAMETER

# Default dataset configuration.
CONFIG = {
    'rmax': 25e-3,
    'sv_num_packets': 1000000,
    'batch_packets': 100000,
    'root_storage_dir': os.path.join(DATASET_PATH, 'data'),

    'sample': {
        'semiinfinite': {
            # 'rmax': 25e-3, # use to overload the top level rmax attribute
            'layers': [
                {'d': 'np.inf', 'n': RI_AIR,   'mua': 0.0,   'mus': 0.0, 'g': 0.0},
                {'d': 'np.inf', 'n': RI_WATER, 'mua': 2.0e2, 'mus': 500.0e2, 'g': 0.95},
                {'d': 'np.inf', 'n': RI_AIR,   'mua': 0.0,   'mus': 0.0, 'g': 0.0}
            ],
            'dir': '1-layer-semiinfinite'   # storage directory name
        },
    },

    'source': {
        'fiber-200um-sds-500um': {
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'UniformFiber',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 200e-6, 'dcladding': 220e-6,
                           'ncore': RI_FUSEDSILICA, 'na': 0.22},
            },
            'args': [],
            'kwargs': {
                'spacing': 500e-6,  # source-detector-separation
                'n': 2,             # have 2 fibers
                'reflectivity': PROBE_REFLECTIVITY,
                'diameter': PROBE_DIAMETER
            },
            'dir': 'fiber-200um-0_22na', # storage directory name
            'n_above': RI_AIR, 'n_bellow': RI_AIR # refractive index for the two outer layers
        },
    },
    'sv': {
        'xaxis':{
            'args': [],
            'kwargs': {'start': -0.75e-3, 'stop': 0.75e-3, 'n': 300},
        },
        'yaxis':{
            'args': [],
            'kwargs': {'start': -0.75e-3, 'stop': 0.75e-3, 'n': 300},
        },
        'zaxis':{
            'args': [],
            'kwargs': {'start': 0.0, 'stop': 1.0e-3, 'n': 200},
        }
    },
    'trace': {
        'kwargs': {'maxlen': 1000}
    },
}

def render_sv_reflectance(target_dir: str = None, config: dict = None,
                          method: str = None,
                          cl_device: str = None, cl_index: int = 0,
                          cl_build_options: List[str] = None,
                          test: bool = False, verbose: bool = False):
    '''
    Render templates with the given configuration.

    Parameters
    ----------
    target_dir: str
        Root directory for the dataset. The scripts will be rendered into
        "run" subdirectory and the dataset data will be saved into the data
        subdirectory. If None, the parent directory of this file will serve
        as the root directory. 
    config: dict
        Configuration / context to use when rendering the run scripts. If None,
        the default configuration will be used.
    method: str
        Monte Carlo stepping method. One of "ar" "aw" or "mbl".
    cl_device: str
        Default OpenCL device name or None. The value can be also set through
        the CL_DEVICE environment variable.
    cl_index: int
        OpenCL device index (if multiple OpenCL devices of the same kind
        are installed). The value can be also set through the CL_INDEX
        environment variable.
    cl_build_options: List[str]
        A list of  OpenCL build options.
        See :py:class:`~xopto.cl.cloptions.ClBuildOption` for more details.
    test: bool
        Do a test run. The run scripts will be rendered but not saved. This
        option will automatically enable the verbose mode.
    verbose: bool
        Enables verbose reporting.
    '''
    if verbose:
        print('Rendering run scripts for Sampling Volume.')

    if config is None:
        config = CONFIG

    if method is None:
        method = 'aw'

    if target_dir is None:
        target_dir = os.getcwd()

    if test:
        verbose = True

    if cl_build_options is None:
        cl_build_options = []
    else:
        cl_build_options = [str(item) for item in cl_build_options]

    root_dataset_dir =  os.path.join(target_dir, 'data')
    run_script_dir = os.path.join(target_dir, 'run', 'sv', 'reflectance')

    if verbose:
        print('Root dataset directory set to "{}".'.format(target_dir))
        print('Rendering run scripts into "{}".'.format(run_script_dir))

    with open(os.path.join(SV_TEMPLATE_PATH, 'fiber',
                           'reflectance.template.py'), 'r') as fid:
        sv_template = fid.read()

    T_sv = jinja2.Template(sv_template)

    for sample_name, sample_data in config['sample'].items():
        for src_name, src_data in CONFIG['source'].items():
            if verbose:
                print('Rendering:', sample_name, src_name)
            rendered_template = T_sv.render(**{
                'sample': sample_data,
                'method': method,
                'cl_device': cl_device, 'cl_index': cl_index,
                'cl_build_options': cl_build_options,
                'rmax': src_data.get('rmax', sample_data.get('rmax', config['rmax'])),
                'sv_num_packets': src_data.get('sv_num_packets', config['sv_num_packets']),
                'batch_packets': config['batch_packets'],
                'sample': sample_data,
                'source': src_data,
                'sv': config['sv'],
                'trace': config['trace'],
                'root_dataset_dir': root_dataset_dir
            })
            filename = os.path.join(
                run_script_dir,
                '{}-{}.py'.format(sample_name, src_name)
            )
            if verbose:
                print('Creating output directory "{}".'.format(
                    run_script_dir))
                print('Saving run script to "{}"'.format(filename))
            if not test:
                os.makedirs(run_script_dir, exist_ok=True)
                with open(filename, 'w') as fid:
                    fid.write(rendered_template)
    if verbose:
        print('The run scripts will save data into "{}".'.format(
            root_dataset_dir))


if __name__ == '__main__':
    parser = common.prepare_cli('Render run scripts for optical '
                                'fiber reflectance Sampling Volume datasets')
    # no additional command line arguments are required
    kwargs = common.process_cli(parser)
    render_sv_reflectance(**kwargs)
