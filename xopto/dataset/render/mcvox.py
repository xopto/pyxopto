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

import numpy as np
import jinja2

from xopto.dataset import DATASET_PATH, MCVOX_TEMPLATE_PATH
from xopto.dataset.render import common
from xopto.dataset.render.mcml import DEFAULT_WAVELENGTH, RI_DIGITS


DEFAULT_WAVELENGTH = common.DEFAULT_WAVELENGTH
RI_DIGITS = common.RI_DIGITS
RI_WATER = round(common.RI_WATER(DEFAULT_WAVELENGTH), RI_DIGITS)

CONFIG = {
    'materials': {
        'surrounding': {'mua':   0.0001e2, 'mus':   1.0000e2, 'g': 1.0, 'n': RI_WATER},
        'epidermis':   {'mua':  16.5724e2, 'mus': 375.9398e2, 'g': 0.9, 'n': RI_WATER},
        'dermis':      {'mua':   0.4585e2, 'mus': 356.5406e2, 'g': 0.9, 'n': RI_WATER},
        'blood':       {'mua': 230.5427e2, 'mus':  93.9850e2, 'g': 0.9, 'n': RI_WATER}
    },
    'voxelization': {'nx': 201, 'ny': 201, 'nz': 200, 'voxel_size': 5.0e-6},
    'epidermis_thickness': 100e-6,
    'vessel_depth': 'np.round(np.linspace(200e-6, 800e-6, 25), 6)',
    'vessel_diameter': 200e-6,
    'num_packets': 1000e6,
    'batch_packets': 10e6,
    'rmax': 25e-3,
    'root_storage_dir': os.path.join(DATASET_PATH, 'data'),
}

def render_mcvox_fluence(target_dir: str = None, config: dict = None,
                         method: str = None, cache: bool = False,
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
    cache: bool
        Enables fluence accumulator cache.
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
        print('Rendering run scripts for MCVOX.')

    if config is None:
        config = CONFIG

    if target_dir is None:
        target_dir = os.getcwd()

    if test:
        verbose = True

    if cl_build_options is None:
        cl_build_options = []
    else:
        cl_build_options = [str(item) for item in cl_build_options]

    root_dataset_dir =  os.path.join(target_dir, 'data')
    run_script_dir = os.path.join(target_dir, 'run', 'mcvox')

    if verbose:
        print('Root dataset directory set to "{}".'.format(target_dir))
        print('Rendering run scripts into "{}".'.format(run_script_dir))

    with open(os.path.join(MCVOX_TEMPLATE_PATH, 
                           'fluence_skin_vessel.template.py'), 'r') as fid:
        sv_template = fid.read()

    T_sv = jinja2.Template(sv_template)

    if verbose:
        print('Rendering: fluence-skin-vessel')
    config = dict(config)
    config['cl_device'] = cl_device
    config['cl_index'] = cl_index
    config['method'] = method
    config['cache'] = cache
    config['root_dataset_dir'] =  root_dataset_dir
    rendered_template = T_sv.render(**config)
    filename = os.path.join(run_script_dir, 'fluence-skin-vessel.py')
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
    render_mcvox_fluence(**kwargs)
