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

from xopto.dataset import MCML_TEMPLATE_PATH, MCVOX_TEMPLATE_PATH
from xopto.dataset.render import common
from xopto.dataset.render.mcml import DEFAULT_WAVELENGTH


ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
SCRIPT_DIR = os.path.join(ROOT_DIR, 'run')

PF_DIR_FORMAT = lambda g: 'g-{:.2f}'.format(g).replace('.', '_')
FILENAME_FORMAT = lambda mua, musr: 'mua-{:.2f}-musr-{:.2f}-invcm'.format(
    mua*1e-2, musr*1e-2).replace('.', '_') + '.npz'

DEFAULT_WAVELENGTH = common.DEFAULT_WAVELENGTH
RI_DIGITS = common.RI_DIGITS
RI_WATER = round(common.RI_WATER(DEFAULT_WAVELENGTH), RI_DIGITS)
RI_FUSED_SILICA = round(common.RI_FUSED_SILICA(DEFAULT_WAVELENGTH), RI_DIGITS)
RI_AIR = round(common.RI_AIR(DEFAULT_WAVELENGTH), RI_DIGITS)

CONFIG = {
    'double_precision': False,
    'rmax': 1000.0e-3,
    'num_packets': 100e6,
    'fluence_k': 0x7FFFFF,

    'sample':{
        '1-layer-100mm': {
            'dir': '1-layer-100mm',
            'layers': [
                {'d': 'np.inf', 'n': RI_AIR, 'mua': 0.0, 'mus': 0.0, 'g': 0.0},
                {'d': 100.0e-3, 'n': RI_WATER, 'mua': 0.0, 'mus': 0.0, 'g': 0.0},
                {'d': 'np.inf', 'n': RI_AIR, 'mua': 0.0, 'mus': 0.0, 'g': 0.0}
            ],
            'mua': 'np.linspace(0.0, 5.0e2, 3)',
            'musr': 'np.linspace(5.0e2, 35.0e2, 3)',
            'g': (0.1, 0.5, 0.9),
        },
        '1-layer-1mm': {
            'dir': '1-layer-1mm',
            'layers': [
                {'d': 'np.inf', 'n': RI_AIR, 'mua': 0.0, 'mus': 0.0, 'g': 0.0},
                {'d': 1.0e-3, 'n': RI_WATER, 'mua': 0.0, 'mus': 0.0, 'g': 0.0},
                {'d': 'np.inf', 'n': RI_AIR, 'mua': 0.0, 'mus': 0.0, 'g': 0.0}
            ],
            'mua': 'np.linspace(0.0, 5.0e2, 3)',
            'musr': 'np.linspace(5.0e2, 35.0e2, 3)',
            'g': (0.1, 0.5, 0.9),
        },
        '2-layer-100um-1mm': {
            'dir': '2-layer-100um-1mm',
            'layers': [
                {'d': 'np.inf', 'n': RI_AIR,
                 'mua': 0.0, 'mus': 0.0, 'g': 0.0},
                {'d': 0.1e-3, 'n': RI_FUSED_SILICA,
                 'mua': 0.0, 'mus': 0.0, 'g': 0.0},
                {'d': 1.0e-3, 'n': RI_WATER,
                 'mua': 0.5e2, 'mus': 20.0e2/(1.0 - 0.8), 'g': 0.8},
                {'d': 'np.inf', 'n': RI_AIR,
                 'mua': 0.0, 'mus': 0.0, 'g': 0.0}
            ],
            'mua': 'np.linspace(0.0, 5.0e2, 3)',
            'musr': 'np.linspace(5.0e2, 35.0e2, 3)',
            'g': (0.1, 0.5, 0.9),
        },
    }
}
'''
MCML comparison dataset configuration.
'''

def render_mcml_comparison(target_dir: str = None, config: dict = None,
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
        print('Rendering run scripts for MCML.')

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
    mcml_run_script_dir = os.path.join(target_dir, 'run', 'mcml_comparison')
    mcvox_run_script_dir = os.path.join(target_dir, 'run', 'mcml_mcvox_comparison')

    if verbose:
        print('Root dataset directory set to "{}".'.format(target_dir))
        print('Rendering mcml run scripts into "{}".'.format(mcml_run_script_dir))
        print('Rendering mcvox run scripts into "{}".'.format(mcvox_run_script_dir))

    template_filename = os.path.join(
        MCML_TEMPLATE_PATH, 'mcml_comparison.template.py')
    with open(template_filename, 'r') as fid:
        mcml_template = fid.read()
    T_mcml = jinja2.Template(mcml_template)

    template_filename = os.path.join(
        MCVOX_TEMPLATE_PATH, 'mcml_mcvox_comparison.template.py')
    with open(template_filename, 'r') as fid:
        mcvox_template = fid.read()
    T_mcvox = jinja2.Template(mcvox_template)

    for sample_name, sample_data in config['sample'].items():
        context = {
            'double_precision': bool(config['double_precision']),
            'fluence_k': config['fluence_k'],
            'root_dataset_dir': root_dataset_dir,
            'rmax': float(config['rmax']),
            'num_packets': int(sample_data.get('num_packets',
                                               config['num_packets'])),
            'sample': sample_data,
            'cl_device': cl_device,
            'cl_index': cl_index,
            'cl_build_options': cl_build_options,
        }

        out_mcml = T_mcml.render(context)
        out_mcvox = T_mcvox.render(context)

        for run_script_dir, out in ((mcml_run_script_dir, out_mcml),
                                    (mcvox_run_script_dir, out_mcvox)):
            rendered_filename = os.path.join(run_script_dir, sample_name + '.py')
            if verbose:
                print('Writing rendered script to "{}".'.format(rendered_filename))
            if not test:
                os.makedirs(run_script_dir, exist_ok=True)
                with open(rendered_filename, 'w') as fid:
                    fid.write(out)

def path(sample: str, g: float, mua: float, musr: float,
         mcvox: bool = False) -> str:
    '''
    Returns path to a dataset file relative to the root directory of the
    dataset.

    Parameters
    ----------
    sample: str
        One of the sample configurations defined in CONFIG.
    g: float
        Scattering phase function anisotropy.
    mua: float
        Absorption coefficient (1/m).
    musr: float
        Reduced scattering coefficient (1/m).
    mcvox: bool
        Indicates if the dataset produced by the voxelized MC should be loaded.

    Returns
    -------
    pats: str
        Path to the dataset file relative to the toot directory of the dataset.
    '''
    if sample not in CONFIG['sample']:
        raise ValueError(
            'The sample argument must be one of {} but got "{}"!'.format(
                CONFIG['sample'].keys(), sample
            )
        )
    if g not in CONFIG['sample'][sample]['g']:
        raise ValueError(
            'The g parameter must be one of {} but got "{}"!'.format(
                CONFIG['sample'][sample]['g'], g
            )
    )
    dir = 'mcml_comparison' if not mcvox else 'mcml_mcvox_comparison'
    data_dir = os.path.join('data', dir, CONFIG['sample'][sample]['dir'],
                            'line', 'radial', 'hg', PF_DIR_FORMAT(g))

    filename = FILENAME_FORMAT(mua, musr)

    return os.path.join(data_dir, filename)

def load(sample: str, g: float, mua: float, musr: float,
         mcvox: bool = False, dataset_dir: str = None)-> dict:
    '''
    Load data for the given sample configuration and layer parameters.

    Parameters
    ----------
    sample: str
        One of the sample configurations defined in CONFIG.
    g: float
        Scattering phase function anisotropy.
    mua: float
        Absorption coefficient (1/m).
    musr: float
        Reduced scattering coefficient (1/m).
    mcvox: bool
        Indicates if the dataset produced by the voxelized MC should be loaded.
    dataset_dir: str
        Root dataset directory (one that contains "data" directory).
        If None or empty, the current working directory is used as
        the dataset directory.

    Returns
    -------
    data: dict
        Sample data.
    '''
    rel_filename = path(sample, g, mua, musr, mcvox)

    full_filename = os.path.join(dataset_dir, rel_filename)

    return dict(np.load(full_filename, allow_pickle=True))

def visualize(sample: str, g: float, mua: float, musr: float,
              mcvox: bool = False, dataset_dir: str = None) -> dict:
    '''
    Load data for the given sample configuration and layer parameters and
    plot the results.

    Parameters
    ----------
    sample: str
        One of the sample configurations defined in SAMPLE_CONTEXT.
    g: float
        Scattering phase function anisotropy.
    mua: float
        Absorption coefficient (1/m).
    musr: float
        Reduced scattering coefficient (1/m).
    mcvox: bool
        Indicates if the dataset produced by the voxelized MC should be loaded.
    dataset_dir: str
        Root dataset directory (one that contains "data" directory).
        The current working directory is used if None.

    Returns
    -------
    data: dict
        Sample data.
    '''
    import matplotlib.pyplot as pp

    if dataset_dir is None:
        dataset_dir =os.getcwd()

    data = load(dataset_dir, sample, g, mua, musr, mcvox)
    fig, [ax1, ax2, ax3] = pp.subplots(1, 3)
    e = np.linspace(0.0, 5.0e-3, 501)
    r = (e[:-1] + e[1:])*0.5

    pp.suptitle('Produced by mcml' if not mcvox else 'Produced by mcvox')

    ax1.semilogy(r*1e3, data['reflectance'])
    ax1.set_xlabel('r (mm)')
    ax1.set_title('Reflectance (total={:.4f})'.format(data['total_reflectance']))

    ax2.semilogy(r*1e3, data['transmittance'])
    ax2.set_xlabel('r (mm)')
    ax2.set_title('Transmittance (total={:.4f})'.format(data['total_transmittance']))

    extent = [0.0, 5.0, 5.0, 0.0]
    im = ax3.imshow(np.log10(data['fluence_data']), extent=extent)
    ax3.set_title('Deposition - log10')
    ax3.set_xlabel('r (mm)')
    ax3.set_ylabel('z (mm)')
    pp.colorbar(im)

    pp.show()

    return data

if __name__ == '__main__':
    parser = common.prepare_cli('Render run scripts for MCML comparison')
    # no additional command line arguments are required
    kwargs = common.process_cli(parser)
    render_mcml_comparison(**kwargs)
