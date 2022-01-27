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
import argparse
from xopto.materials import ri


DEFAULT_WAVELENGTH = 500e-9

RI_DIGITS = 3
RI_FUSED_SILICA = ri.glass.fusedsilica.malitson
RI_WATER = ri.water.daimon
RI_POLYSTYRENE = ri.polystyrene.nikolov
RI_AIR = lambda wavelength: 1.0

PROBE_REFLECTIVITY = lambda wavelength: 0.6
PROBE_DIAMETER = 6.0e-3

def prepare_cli(description: str) -> argparse.ArgumentParser:
    '''
    Prepare argument parser for the rendering calls.

    Parameters
    ----------
    description: str
        Description for the command line interface.

    Returns
    -------
    parser: argparse.ArgumentParser

    Notes
    -----
    Command line argument parser instance implements the following arguments:

        verbose: bool
        Turns on verbose mode.

        test: bool
            Run in test mode.

        output_dir: str
            Root directory of the dataset.

        cl_device: str
            OpenCL device name to use.

        cl_index: int
            OpenCL device index (zero-based).
    '''

    cwd = os.getcwd()

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        '-o', '--target-dir', dest='target_dir', default=cwd, type=str,
        help='Root directory of the dataset. The run scripts will be rendered '
             'into the "run/mcml/" subdirectory.')
    parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='Enables verbose report mode.')
    parser.add_argument(
        '-t', '--test', dest='test', action='store_true',
        help='Do a test run with verbose output. No directories will '
             'be created and no files will be saved.')
    parser.add_argument(
        '-d', '--cl-device', dest='cl_device', default='', type=str,
        help='OpenCL device to use for the simulations.')
    parser.add_argument(
        '-i', '--cl-index', dest='cl_index', default=0, type=int,
        help='OpenCL device index to use for the simulations.')

    parser.add_argument(
        '--cl-fast-relaxed-math', dest='cl_fast_relaxed_math',
        action='store_true',
        help='Adds OpenCL build option -cl-fast-relaxed-math.')
    parser.add_argument(
        '--cl-mad-enable', dest='cl_mad_enable',
        action='store_true',
        help='Adds OpenCL build option -cl-mad-enable.')
    parser.add_argument(
        '--cl-unsafe-math-optimizations', dest='cl_unsafe_math_optimizations',
        action='store_true',
        help='Adds OpenCL build option -cl-unsafe-math-optimizations.')
    parser.add_argument(
        '--cl-finite-math-only', dest='cl_finite_math_only',
        action='store_true',
        help='Adds OpenCL build option -cl-finite-math-only.')
    parser.add_argument(
        '--cl-no-signed-zeros', dest='cl_no_signed_zeros',
        action='store_true',
        help='Adds OpenCL build option -cl-no-signed-zeros.')

    return parser

def process_cli(parser: argparse.ArgumentParser) -> dict:
    '''
    Process the ArgumentParser instance created with :py:func:`prepare_cli`.

    Parameters
    ----------
    parser: argparse.ArgumentParser
        Argument parser as returned by the :py:func:`prepare_ci` function.

    Returns
    -------
    kwargs: dict
        Keyword arguments from the command line interface as a dict.

    Notes
    -----
    The command line implements the following keyword arguments:

        verbose: bool
            Turns on verbose mode.

        test: bool
            Run in test mode.

        output_dir: str
            Root directory of the dataset.

        cl_device: str
            OpenCL device name to use.

        cl_index: int
            OpenCL device index (zero-based).

        cl_build_options: List[str]
            A list of  OpenCL build options.
            See :py:class:`~xopto.cl.cloptions.ClBuildOption` for more
            details.
    '''
    args = parser.parse_args()

    cl_build_options = []
    if args.cl_fast_relaxed_math:
        cl_build_options.append('-cl-fast-relaxed-math')
    if args.cl_mad_enable:
        cl_build_options.append('-cl-mad-enable')
    if args.cl_unsafe_math_optimizations:
        cl_build_options.append('-cl-unsafe-math-optimizations')
    if args.cl_finite_math_only:
        cl_build_options.append('-cl-finite-math-ony')
    if args.cl_no_signed_zeros:
        cl_build_options.append('-cl-no-signed-zeros')

    return {'verbose': args.verbose, 'test': args.test,
            'target_dir': args.target_dir,
            'cl_device': args.cl_device,
            'cl_index': args.cl_index,
            'cl_build_options': cl_build_options}
