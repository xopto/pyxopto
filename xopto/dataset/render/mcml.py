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
import numpy as np

from xopto.mcml import mc
from xopto.dataset import DATASET_PATH, MCML_TEMPLATE_PATH
from xopto.dataset.render import common

# Helper function that converts FWHM to sigma.
fwhm2sigma = mc.mcsource.GaussianBeam.fwhm2sigma

# Number of samples along the absorption and reduced scattering coefficient axes.
N_MUA = 21
N_MUSR = 21

# Spatial frequency domain - frequencies (1/m)
SFDI_FREQUENCY_RANGE = [0.00e3, 0.80e3]
SFDI_FREQUENCY_STEP = 0.01e3

# Data common to all datasets.
RI_DIGITS = common.RI_DIGITS
DEFAULT_WAVELENGTH = common.DEFAULT_WAVELENGTH

# General data for all datasets.
RI_AIR = round(common.RI_AIR(DEFAULT_WAVELENGTH), RI_DIGITS)
RI_WATER = round(common.RI_WATER(DEFAULT_WAVELENGTH), RI_DIGITS)
RI_FUSEDSILICA = round(common.RI_FUSED_SILICA(DEFAULT_WAVELENGTH), RI_DIGITS)
PROBE_REFLECTIVITY = common.PROBE_REFLECTIVITY(DEFAULT_WAVELENGTH)
PROBE_DIAMETER = common.PROBE_DIAMETER

# Data for the MIE datasets.
MIE_RI_DIGITS          = common.RI_DIGITS
MIE_WAVELENGTHS        = [DEFAULT_WAVELENGTH]
MIE_RI_POLYSTYRENE     = [round(common.RI_POLYSTYRENE(w), MIE_RI_DIGITS) for w in MIE_WAVELENGTHS]
MIE_RI_WATER           = [round(common.RI_WATER(w), MIE_RI_DIGITS) for w in MIE_WAVELENGTHS]
MIE_RI_FUSED_SILICA    = [round(common.RI_FUSED_SILICA(w), MIE_RI_DIGITS) for w in MIE_WAVELENGTHS]
MIE_RI_FIBER_CORE      = [round(common.RI_FUSED_SILICA(w), MIE_RI_DIGITS) for w in MIE_WAVELENGTHS]
MIE_DIAMETERS          = [0.25e-6, 0.5e-6, 1.0e-6, 2.0e-6, 4.0e-6]
MIE_PROBE_REFLECTIVITY = [PROBE_REFLECTIVITY]*len(MIE_WAVELENGTHS)

# Data for the SFDI dataset.
SFDI_ACCEPTANCE = np.deg2rad(10.0)
SFDI_INCIDENCE = np.deg2rad(20.0)

# Default dataset configuration.
CONFIG = {
    'rmax': 25e-3,
    'num_packets': 100e6,
    'root_storage_dir': os.path.join(DATASET_PATH, 'data'),

    'spatial_frequencies': 'np.arange({}, {}, {})'.format(
        SFDI_FREQUENCY_RANGE[0],
        SFDI_FREQUENCY_RANGE[1] + 0.5*SFDI_FREQUENCY_STEP,
        SFDI_FREQUENCY_STEP),

    'sample': {
        'semiinfinite': {
            # 'rmax': 25e-3, # use to overload the top level rmax attribute
            'layers': [
                {'d': 'np.inf', 'n': RI_AIR,   'mua': 0.0, 'mus': 0.0},
                {'d': 'np.inf', 'n': RI_WATER, 'mua': 0.0, 'mus': 0.0},
                {'d': 'np.inf', 'n': RI_AIR,   'mua': 0.0, 'mus': 0.0}
            ],
            'n_mua': N_MUA, 'n_musr': N_MUSR,   # number of samples along the mua and musr axis
            'mua': 'np.linspace(0.0, 5.0e2, {:d})'.format(N_MUA),       # for 1st layer only
            'musr': 'np.linspace(5.0e2, 35.0e2, {:d})'.format(N_MUSR),  # for 1st layer only
            'dir': '1-layer-semiinfinite'   # storage directory name
        },
    },

    'source': {
        'line': { # radial detectors on top and bottom, NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'Line',
            'args': [], 'kwargs': {},
            'dir': 'line', # storage directory name
            'n_above': RI_AIR, 'n_bellow': RI_AIR # refractive index for the two outer layers
        },
        'collimated-200um': { # radial detectors on top and bottom, NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'UniformBeam',
            'args': [], 'kwargs': {'diameter': 200.0e-6},
            'dir': 'collimated-200um', # storage directory name
            'n_above': RI_AIR, 'n_bellow': RI_AIR # refractive index for the two outer layers
        },
        'gaussian-100um': { # radial detectors on top and bottom, NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'GaussianBeam',
            'args': [], 'kwargs': {'sigma': fwhm2sigma(100e-6)},
            'dir': 'gaussian-fwhm-100um', # storage directory name
            'n_above': RI_AIR, 'n_bellow': RI_AIR # refractive index for the two outer layers
        },
        'fiber-200um': { # radial detector on top with with fiber NA, bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'UniformFiber',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 200e-6, 'dcladding': 220e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [], 'kwargs': {},
            'n_above': RI_FUSEDSILICA, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'fiber-200um-0.22na' # storage directory name
        },
        'six-around-one-200um': { # surface layout + SixAroundOne detector on top, radial on bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'SixAroundOne',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 200e-6, 'dcladding': 220e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [],
            'kwargs': {'diameter': PROBE_DIAMETER, 'spacing': 220e-6,
                       'reflectivity': PROBE_REFLECTIVITY},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'six-around-one-200um-0_22na' # storage directory name
        },
        'six-around-one-400um': { # surface layout + SixAroundOne detector on top, radial on bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'SixAroundOne',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 400e-6, 'dcladding': 420e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [],
            'kwargs': {'diameter': PROBE_DIAMETER, 'spacing': 420e-6,
                       'reflectivity': PROBE_REFLECTIVITY},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'six-around-one-400um-0_22na' # storage directory name
        },
        'six-linear-array-200um': { # surface layout + LinearArray detector on top, radial on bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'LinearArray',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 200e-6, 'dcladding': 220e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [],
            'kwargs': {'n': 6, 'diameter': PROBE_DIAMETER, 'spacing': 220e-6,
                       'reflectivity': PROBE_REFLECTIVITY},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'six-linear-200um-0_22na' # storage directory name
        },
        'single-fiber-100um': { # surface layout + LinearArray detector (1 fiber) on top, radial on bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'LinearArray',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 100e-6, 'dcladding': 120e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [],
            'kwargs': {'n': 1, 'diameter': PROBE_DIAMETER, 'spacing': 120e-6,
                       'reflectivity': PROBE_REFLECTIVITY},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'single-fiber-100um-0_22na' # storage directory name
        },
        'single-fiber-200um': { # surface layout + LinearArray detector (1 fiber) on top, radial on bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'LinearArray',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 200e-6, 'dcladding': 220e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [],
            'kwargs': {'n': 1, 'diameter': PROBE_DIAMETER, 'spacing': 220e-6,
                       'reflectivity': PROBE_REFLECTIVITY},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'single-fiber-200um-0_22na' # storage directory name
        },
        'single-fiber-400um': { # surface layout + LinearArray detector (1 fiber) on top, radial on bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'LinearArray',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 400e-6, 'dcladding': 420e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [],
            'kwargs': {'n': 1, 'diameter': PROBE_DIAMETER, 'spacing': 420e-6,
                       'reflectivity': PROBE_REFLECTIVITY},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'single-fiber-400um-0_22na' # storage directory name
        },
        'single-fiber-800um': { # surface layout + LinearArray detector (1 fiber) on top, radial on bottom with NA=1
            # 'num_packets': 100e6, # use to overload the top level num_packets attribute
            'type': 'LinearArray',
            'fiber': {
                'args': [],
                'kwargs': {'dcore': 800e-6, 'dcladding': 820e-6, 'na': 0.22,
                           'ncore': RI_FUSEDSILICA}
            },
            'args': [],
            'kwargs': {'n': 1, 'diameter': PROBE_DIAMETER, 'spacing': 820e-6,
                       'reflectivity': PROBE_REFLECTIVITY},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'single-fiber-800um-0_22na' # storage directory name
        },
        'sfdi-normal':{ # normal projection, normal detector
            'num_packets': 1000e6, # use to overload the top level num_packets attribute
            'rmax': 150e-3,     # overloads the top level and sample rmax attribute
            'is_sfdi': True,    # indicate that this is a SFDI source/detector configuration
            'detector': 'radial_sfdi',
            'type': 'Line',
            'args': [], 'kwargs': {'direction': [0.0, 0.0, 1.0]},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'sfdi-0deg-incidence' # storage directory name
        },
        'sfdi-normal-tilted':{ # normal projection, tilted detector
            'num_packets': 1000e6, # use to overload the top level num_packets attribute
            'rmax': 150e-3,     # overloads the top level and sample rmax attribute
            'is_sfdi': True,    # indicate that this is a SFDI source/detector configuration
            'type': 'Line',
            'detector': 'symmetricx_sfdi',
            'args': [],
            'args': [], 'kwargs': {'direction': [0.0, 0.0, 1.0]},
            'n_above': RI_AIR, 'n_bellow': RI_AIR, # refractive index for the two outer layers
            'dir': 'sfdi-0deg-incidence' # storage directory name
        }
    },
    'pf': {
        'hg': {
            'type': 'Hg',
            'g': [0.1, 0.3, 0.5, 0.7, 0.9], 'default': [0.0],
            'dir': 'hg', # pf dir
            'param_dir': ['g-{:.2f}'] # templates for parameter dirs
        },
        'gk': {
            'type': 'Gk',
            'g': [ 0.1, 0.3, 0.5, 0.7, 0.9],
            'a': [-0.5, 0.0, 0.5, 1.0, 1.5], 'default': [0.0, 0.0],
            'dir': 'gk', # pf dir
            'param_dir':['g-{:.2f}', 'a-{:.2f}'] # templates for parameter dirs
        },
        'mhg': {
            'type': 'MHg',
            'g': [0.1, 0.3, 0.5, 0.7, 0.9],
            'b': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], 'default': [0.0, 0.0],
            'dir': 'mhg', # pf dir
            'param_dir':[ 'g-{:.2f}', 'b-{:.2f}'] # templates for parameter dirs
        },
        'mpc': {
            'type': 'MPc',
            'n': [0.0, 1.0, 2.0, 4.0, 8.0],
            'b': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], 'default': [0.0, 0.0],
            'dir': 'mpc',   # pf dir
            'param_dir': ['n-{:.2f}', 'b-{:.2f}'] # templates for parameter dirs
        },
        'mie-sphere-polystyrene':{
            'lutsize': 4000,    # 2000 is the default value
            'diameters': MIE_DIAMETERS,
            'wavelengths': MIE_WAVELENGTHS,
            'nmedium': MIE_RI_WATER,
            'nparticle': MIE_RI_POLYSTYRENE,
            'nfibercore': MIE_RI_FIBER_CORE,
            'probe_reflectivity': MIE_PROBE_REFLECTIVITY,
            'dir': 'mie-polystyrene',
            'param_dir': ['diameter-{:.2f}um', 'wavelength-{:.0f}nm'],
            'default': [1.0, 1.0, 1e-6, 500e-9] # default parameter for initialization
        },
        'mie-sphere-fusedsilica':{
            'lutsize': 4000,    # 2000 is the default value
            'diameters': MIE_DIAMETERS,
            'wavelengths': MIE_WAVELENGTHS,
            'nmedium': MIE_RI_WATER,
            'nparticle': MIE_RI_FUSED_SILICA,
            'nfibercore': MIE_RI_FIBER_CORE,
            'probe_reflectivity': MIE_PROBE_REFLECTIVITY,
            'dir': 'mie-fusedsilica',
            'param_dir': ['diameter-{:.2f}um', 'wavelength-{:.0f}nm'],
            'default': [1.0, 1.0, 1e-6, 500e-9] # default parameter for initialization
        }
    },
    'detector':{
        'radial': {
            'top': {
                'raxis': {
                    'args': [],
                    'kwargs': {'start': 0.0, 'stop': 5.0e-3, 'n': 500}
                },
                'args': [], 'kwargs': {}
            },
            'bottom': {
                'raxis': {
                    'args': [],
                    'kwargs': {'start': 0.0, 'stop': 5.0e-3, 'n': 500}
                },
                'args': [], 'kwargs': {}
            },
            'dir': 'radial'
        },
        'radial_sfdi': {
            'top': {
                'raxis': {
                    'args': [],
                    'kwargs': {'start': 0.0, 'stop': 150.0e-3, 'n': 4000,
                               'logscale': True}
                },
                'args': [],
                'kwargs': {
                    'cosmin': np.cos(SFDI_ACCEPTANCE)
                }
            },
            'bottom': {
                'raxis': {
                    'args': [],
                    'kwargs': {'start': 0.0, 'stop': 150.0e-3, 'n': 4000,
                               'logscale': True}
                },
                'args': [],
                'kwargs': {
                    'cosmin': np.cos(SFDI_ACCEPTANCE)
                }
            },
            'dir': 'radial'
        },
        'symmetricx_sfdi': {
            'top': {
                'xaxis': { # SymmetricAxis
                    'args': [],
                    'kwargs': {'center': 0.0, 'range': 150.0e-3, 'n_half': 4000,
                               'logscale': True}
                },
                'args': [],
                'kwargs': {
                    'cosmin': np.cos(SFDI_ACCEPTANCE),
                    'direction': [
                        np.sin(SFDI_INCIDENCE),
                        0.0,
                        np.cos(SFDI_INCIDENCE)
                    ]
                }
            },
            'bottom': { # SymmetricAxis
                'raxis': {
                    'args': [],
                    'kwargs': {'center': 0.0, 'range': 150.0e-3,
                               'n_half': 4000, 'logscale': True}
                },
                'args': [],
                'kwargs': {
                    'cosmin': np.cos(SFDI_ACCEPTANCE)
                }
            },
            'dir': 'symmetricx-{:.0f}deg'.format(np.rad2deg(SFDI_INCIDENCE))
        }
    }
}

def render_mcml(target_dir: str = None, config: dict = None,
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
        print('Rendering run scripts for MCML.')

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
    run_script_dir = os.path.join(target_dir, 'run', 'mcml')

    if verbose:
        print('Root dataset directory set to "{}".'.format(target_dir))
        print('Rendering run scripts into "{}".'.format(run_script_dir))

    with open(os.path.join(MCML_TEMPLATE_PATH,
                           'analytical_pf.template.py'), 'r') as fid:
        analytical_pf_template = fid.read()
    with open(os.path.join(MCML_TEMPLATE_PATH,
                           'mie_pf.template.py'), 'r') as fid:
        mie_pf_template = fid.read()

    T_analytical_pf = jinja2.Template(analytical_pf_template)
    T_mie_pf = jinja2.Template(mie_pf_template)

    for sample_name, sample_data in config['sample'].items():
        for src_name, src_data in CONFIG['source'].items():
            for spf_name, spf_data in CONFIG['pf'].items():
                if verbose:
                    print('Rendering:', sample_name, src_name, spf_name)
                T = T_mie_pf if spf_name.startswith('mie-') else T_analytical_pf
                rendered_template = T.render(**{
                    'sample': sample_data,
                    'method': method,
                    'cache': cache,
                    'cl_device': cl_device, 'cl_index': cl_index,
                    'cl_build_options': cl_build_options,
                    'rmax': src_data.get('rmax', sample_data.get('rmax', config['rmax'])),
                    'num_packets': src_data.get('num_packets', config['num_packets']),
                    'sample': sample_data,
                    'source': src_data,
                    'detector': config['detector'],
                    'pf': spf_data,
                    'root_dataset_dir': root_dataset_dir,
                    'spatial_frequencies': config['spatial_frequencies']
                })
                filename = os.path.join(
                    run_script_dir,
                    '{}-{}-{}.py'.format(sample_name, src_name, spf_name)
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
    parser = common.prepare_cli('Render run scripts for the MCML datasets')
    # no additional command line arguments are required
    kwargs = common.process_cli(parser)
    render_mcml(**kwargs)
