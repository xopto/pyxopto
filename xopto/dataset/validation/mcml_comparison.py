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
import sys

import numpy as np
import matplotlib.pyplot as pp

from xopto.dataset.render.mcml_comparison import load, path, CONFIG
from xopto.mcml import mc

def _prepare(item):
    if isinstance(item, str):
        item = eval(item)
    return item

def _init_data():
    return {'mre': np.zeros((3, 3, 3)), 'rrmse': np.zeros((3, 3, 3))}

def _relative_error(test, reference):
    delta = test - reference
    if delta.size > 1:
        nonzero_mask = reference != 0.0
        relative = np.tile(np.nan, delta.shape)
        if nonzero_mask.size:
            relative[nonzero_mask] = delta[nonzero_mask]/reference[nonzero_mask]
        relative[delta == 0.0] = 0.0
    else:
        relative = float(np.nan)
        if reference != 0.0:
            relative = float(delta/reference)
        if delta == 0.0:
            relative = 0.0
    return relative

def _fluence_proj_re(test, reference, fluencerz):
    e = fluencerz.raxis.edges
    k = np.pi*(e[1:]**2 - e[:-1]**2)*fluencerz.dz
    k.shape = (1, k.size)
    test_raw = test*k
    reference_raw = reference*k
    rel_r = _relative_error(test_raw.sum(axis=0), reference_raw.sum(axis=0))
    rel_z = _relative_error(test_raw.sum(axis=1), reference_raw.sum(axis=1))
    return rel_r, rel_z

def compare(test_dir: str, reference_dir: str, out_dir: str = None,
            verbose: bool = False):
    '''
    Compares the results of two "MCML comparison" datasets that include data
    on several configurations computed with layered (mcml) and voxelized
    (mcvox) Monte Carlo simulator cores.

    Parameters
    ----------
    test_dir: str
        Path to the root directory of the test dataset.
    reference_dir: str
        Path to the root directory of the reference dataset.
    out_dir: str
        Optional directory for producing error plots. The error plots follow
        the directory structure of the dataset.
    verbose: bool
        Enable verbose reporting to stdout.

    Returns
    -------
    results: dict
        A collection of Relative Root-Mean-Square error (RRMSE) and
        Mean Relative Error (MRE) values for all the dataset configurations
        defined in the :py:attr:`~xopto.dataset.render.mcml_comparison.CONFIG`
        dict. The MRE and RRMSE values are given in multidimensional numpy
        arrays. The following code illustrates the organization of the
        returned data. note that the fields are fully expanded only for
        configuration :code:`1-layer-100mm` simulated with the Monte
        Carlo core for layered media (:code:`mcml`).

        .. code-block:: python

            {
                'mcml': {
                    '1-layer-100mm': {
                        'transmittance': {'mre': values, 'rrmse': values},
                        'reflectance': {'mre': values, 'rrmse': values},
                        'total_transmittance': {'mre': values, 'rrmse': values},
                        'total_reflectance': {'mre': values, 'rrmse': values},
                        'fluence_data': {'mre': values, 'rrmse': values},
                        'indexing': ('mua', 'musr', 'g')
                        'mua': mua_values,
                        'musr': msr_values,
                        'g': g_values
                    },
                    '1-layer-1mm': {...},
                    '2-layer-100um-1mm': {...}
            },
                'mcvox': {
                    '1-layer-100mm': {...},
                    '1-layer-1mm': {...},
                    '2-layer-100um-1mm': {...}
                }
            }

    '''
    raxis = mc.mcdetector.Axis(0.0, 5.0e-3, 500)
    r_mm = raxis.centers*1e3

    if verbose:
        print()

    plot = out_dir is not None
    results = {'mcvox':{}, 'mcml':{}}
    for mcvox, name in ((False, 'mcml'), (True, 'mcvox')):
        for key, sample in CONFIG['sample'].items():
            results[name][key] = {
                'reflectance': _init_data(),
                'transmittance': _init_data(),
                'total_reflectance': _init_data(),
                'total_transmittance': _init_data(),
                'specular': _init_data(),
                'indexing': ('mua', 'musr', 'g'),
                'fluence_data': _init_data()
            }
            rrmse = results[name][key]

            mua_values = np.asarray(_prepare(sample['mua']))
            musr_values = np.asarray(_prepare(sample['musr']))
            g_values = np.asarray(_prepare(sample['g']))
            rrmse.update(
                {'mua': mua_values, 'musr': musr_values, 'g': g_values})

            for mua_ind, mua in enumerate(mua_values):
                for musr_ind, musr in enumerate(musr_values):
                    for g_ind, g in enumerate(g_values):
                        test = load(key, g=g, mua=mua, musr=musr,
                                    mcvox=mcvox, dataset_dir=test_dir)
                        reference = load(key, g=g, mua=mua, musr=musr,
                                        mcvox=mcvox, dataset_dir=reference_dir)

                        rel_path = path(key, g=g, mua=mua, musr=musr,
                                        mcvox=mcvox)
                        out_file = os.path.splitext(rel_path)[0] + '.pdf'
                        
                        if verbose:
                            print('Processing: "{}"'.format(rel_path), end='\r')

                        if plot:
                            fig, ax = pp.subplots(2, 2, num=key)
                            ax[0, 0].set_title('Relative error (%)')
                            ax[0, 0].set_ylim(-5.0, 5.0)
                            ax[0, 0].set_xlabel('$r$ (mm)')
                            ax[0, 1].set_title('Deposition error (%)')
                            ax[0, 1].set_xlabel('$r$ (mm)')
                            ax[0, 1].set_ylabel('$z$ (mm)')
                            ax[1, 1].set_xlabel('$r, z$ (mm)')
                            ax[1, 1].set_ylim(-2.0, 2.0)

                        for item in ('fluence_data', 'total_reflectance',
                                    'reflectance', 'total_transmittance',
                                    'transmittance', 'specular'):

                            if item == 'fluence_data' and mua == 0.0:
                                continue

                            if item in ('reflectance', 'transmittance',
                                        'fluence_data', 'total_reflectance',
                                        'total_transmittance', 'specular'):
                                deltar = _relative_error(
                                    test[item], reference[item])
                                if np.any(np.isfinite(deltar)):
                                    rrmse_value = np.sqrt(np.nanmean(deltar**2))
                                    mre_value = np.nanmean(deltar)
                                else:
                                    rrmse_value = np.nan
                                    mre_value = np.nan

                            rrmse[item]['rrmse'][mua_ind, musr_ind, g_ind] = \
                                rrmse_value
                            rrmse[item]['mre'][mua_ind, musr_ind, g_ind] = \
                                mre_value

                            if item in ('reflectance', 'transmittance'):
                                if plot:
                                    ax[0, 0].plot(r_mm, deltar*1e2, label=item)
                            elif item in ('total_reflectance',
                                          'total_transmittance', 'specular'):
                                if plot:
                                    s = item.replace('_', ' ').capitalize()
                                    posy = {'total_reflectance': 0.75,
                                            'specular': 0.5,
                                            'total_transmittance': 0.25}[item]
                                    ax[1, 0].text(
                                        0.5, posy,
                                        '{} error: {:3f}%'.format(s, deltar*1e2),
                                        fontsize=8,
                                         verticalalignment='center',
                                         horizontalalignment='center',)
                            elif item == 'fluence_data':
                                if plot:
                                    fluencerz = mc.mcfluence.FluenceRz.fromdict(
                                        reference['fluence'].item())
                                    r_rel, z_rel = _fluence_proj_re(
                                        test[item], reference[item], fluencerz)
                                    im = ax[0, 1].imshow(
                                        deltar*1e2, vmin=-5.0, vmax=5.0,
                                        extent=[fluencerz.raxis.span[0]*1e3,
                                                fluencerz.raxis.span[1]*1e3,
                                                fluencerz.zaxis.span[1]*1e3,
                                                fluencerz.zaxis.span[0]*1e3])
                                    pp.colorbar(im, ax=ax[0, 1])

                                    ax[1, 1].plot(
                                        fluencerz.r*1e3, 1e2*r_rel,
                                        label='$r$-projection')
                                    ax[1, 1].plot(
                                        fluencerz.z*1e3, 1e2*z_rel,
                                        label='$z$-projection')

                                    ax[1, 1].legend()
                                    ax[1, 1].grid()

                        if plot:
                            ax[0, 0].grid()
                            ax[0, 0].legend()
                            pp.tight_layout()
                            full_filename = os.path.join(out_dir, out_file)
                            os.makedirs(os.path.dirname(full_filename),
                                        exist_ok=True)
                            fig.savefig(full_filename)
                            pp.close(fig)

    if verbose:
        print()

    return results

if __name__ == '__main__':
    import argparse

    cwd = os.getcwd()

    parser = argparse.ArgumentParser(
        description='Validation of MCML comparison dataset.')

    parser.add_argument(
        '-t', '--test', dest='test_dir', default=cwd, type=str,
        help='Root directory of the test dataset.')
    parser.add_argument(
        '-r', '--reference', dest='reference_dir', default=cwd, type=str,
        help='Root directory of the reference dataset.')
    parser.add_argument(
        '-o', '--output', dest='output_dir', default=cwd, type=str,
        help='Root directory for the comparison output.')

    parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='Enables verbose report mode.')


    args = parser.parse_args()

    test_dir = args.test_dir
    reference_dir = args.reference_dir
    output_dir = args.output_dir
    verbose = bool(args.verbose)

    if test_dir == reference_dir:
        raise ValueError('Directories of the reference and test dataset '
                         'must not be the same!')

    result = compare(test_dir, reference_dir, output_dir, verbose=verbose)

    np.savez_compressed(os.path.join(output_dir, 'validation.npz'), **result)
