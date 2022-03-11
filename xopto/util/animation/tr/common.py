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

import os.path
import argparse
import matplotlib as mpl
from matplotlib.animation import FuncAnimation

import numpy as np

from xopto.util.animation.common import create_frame_animation

from xopto.mcml import mc
from xopto.cl import clinfo


def create_deposition_animation_xz(
        data_source: str or dict, filename: str = None,
        overwrite=False, logscale=False,
        fps: float = None, duration: float = None,
        axis_off: bool = False, title: str = None,
        autoscale: bool = True, writer: str or mpl.animation.MovieWriter = None,
        dpi: int = None, verbose: bool = False) -> FuncAnimation:
    '''
    Creates a continuous animation of the sampling volume projected
    onto the x-z plane.

    Parameters
    ----------
    data_source: str or dict
        File with the simulation results or simulation results as returned by
        :py:func:`sv_fiber_incidence`.
    filename: str
        Output file for the animation.
    overwrite: bool
        If False, existing animation files will not be overwritten.
    logscale: bool
        If True, logarithmically scale the sampling volume data.
    fps: float
        Frame rate of the animation (1/s). Another way to control the frame
        rate is specify the duration parameter.
    duration: float
        Duration of one full animation (s). The same effect can be achieved
        with the fps parameter.
    axis_off: bool 
        Turns off the axis and related labels.
    title: str
        Plot title.
    autoscale: bool
        Independently autoscale the intensity of each frame.
    writer: str or mpl.animation.MovieWriter
        Movie writer.
    dpi: int
        Resolution of the exported images.
    verbose: bool
        Turns on verbose progress report.

    Returns
    -------
    animation: FuncAnimation
        Animation instance
    '''
    if isinstance(data_source, str):
        data = np.load(data_source, allow_pickle=True)
    elif isinstance(data_source, dict):
        data = data_source
    else:
        raise TypeError('Data source must be a dict or a filename!')

    if writer is None:
        writer = 'imagemagick'

    if filename is None and isinstance(data_source, str):
        filename = os.path.splitext(data_source)[0] + '.gif'

    if filename is not None and os.path.isfile(filename) and not overwrite:
        if verbose:
            print('Animation file "{}" already exists!'.format(filename))
            print('Animation will not be saved!')
        filename = None

    fluence_dict = data['fluence']
    if not isinstance(fluence_dict, dict):
        fluence_dict = fluence_dict.item()

    fluence = mc.mcfluence.Fluencet.fromdict(fluence_dict)
    x_range = fluence.xaxis.edges[0]*1e3, fluence.xaxis.edges[-1]*1e3
    z_range = fluence.zaxis.edges[-1]*1e3, fluence.zaxis.edges[0]*1e3

    if title == 'time':
        title = ['{:.1f} ps'.format(item*1e12) for item in fluence.t]

    extent = [x_range[0], x_range[1], z_range[0], z_range[1]]

    data = np.array(data['fluence_data'])
    if data.ndim >= 4:
        data = data.mean(1)
    data = np.transpose(data, [2, 0, 1])

    frames = data

    return create_frame_animation(
        frames, filename, overwrite=overwrite, logscale=logscale,
        extent=extent, xrange=x_range, yrange=z_range, title=title,
        xlabel='$x$ (mm)', ylabel='$z$ (mm)', axis_off=axis_off,
        autoscale=autoscale, fps=fps, duration=duration,
        writer=writer, verbose=verbose,
        cbar=True, cbar_tick_format='{:.0f}',
        cbar_label='($\\mathrm{{m}}^{{3}}\\mathrm{{s}})^{-1}$',
        dpi=dpi)


def mc_tr_deposition_cli(description: str = None,
                         parser: argparse.ArgumentParser = None):
    '''
    Collects command line arguments and returns a dict of values.

    Parameters
    ----------
    description: str
        Command line description.
    parser: argparse.ArgumentParser
        Externally created argument parser can be used to add custom command
        line arguments.

    Returns
    -------
    args: argparse.Namespace
        Parsed command line arguments. If the parser was passed
        as an input argument, use the returned namespace to extract
        any custom parameters.
    arguments: dict
        Returns the collected arguments as a dict with the following keys:
    
        - num_packets: int
            Minimum number of packets that should reach the detector.
        - batch_packets:
            int Number of packets to launch in a single simulation run.
        - filename: str
            File for the simulation results.
        - overwrite: bool
            Overwrite existing file.
        - verbose: bool
            Enables verbose reports.
        - animation: bool
            Creates animation after completing the sampling volume simulations.
        - fps: float
            Frame rate (1/s) to use in the sampling volume animation.
        - duration: float
            Animation duration (s).
        - autoscale: bool
            Independently autoscale the intensity of each frame.
        - logscale: bool
            Apply logarithmic scale to the deposition data.
        - axis_off: bool
            Turns off the plot axis in animation.
        - dpi: int
            Resolution of the exported images.
        - opencl_device: str
            OpenCL device to use for the simulations.
        - opencl_index: int
            Index (zero-based) of the OpenCL device to use for the simulations.
    '''
    if description is None:
        description = 'Sampling volume simulation'

    parser = argparse.ArgumentParser(description='description')
    parser.add_argument(
        '-n', '--num-packets', dest='num_packets', default=1000, type=int,
        help='Minimum number of packets reaching the detector.')
    parser.add_argument(
        '-p', '--projection', dest='sv_projection', action='store_true',
        help='Reduce the number of sampling volume accumulators '
             'along the y axis to 1.')
    parser.add_argument(
        '-b', '--batch-packets', dest='batch_packets', default=10000000,
        type=int,
        help='Number of photon packets to launch in a single Monte Carlo '
             'simulation.')
    parser.add_argument(
        '-f', '--file', dest='filename', default='', type=str,
        help='File for the simulation results. Note that the extension .npz '
             'is addedd automatically.')
    parser.add_argument(
        '-o', '--overwrite', action='store_true',
        help='Set this flag to force overwrite of existing files.')
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Enables verbose progress report.')
    parser.add_argument(
        '-a', '--animation', action='store_true',
        help='Create animation after completing the Monte Carlo simulations.')
    parser.add_argument(
        '--title', dest='title', type=str, default='',
        help='Animation plot title.')
    parser.add_argument(
        '--axis-off', dest='axis_off', action='store_true',
        help='Turn off animation plot axis and related labels.')
    parser.add_argument(
        '--fps', dest='fps', type=float, default=0.0,
        help='Frame rate of the animation. Alternatively, use the duration '
             'parameter to set the FPS.')
    parser.add_argument(
        '--duration', dest='duration', type=float, default=0.0,
        help='Duration (s) of one full animation. Alternatively, use the fps '
             'parameter to set the duration of the animation.')
    parser.add_argument(
        '-l', '--logscale', dest='logscale', action='store_true',
        help='Log scale intensities in animations.')
    parser.add_argument(
        '--no-autoscale', dest='no_autoscale', action='store_true',
        help='Disable independent automatic intensity scaling of the '
             'individual frames in the animation.')
    parser.add_argument(
        '--dpi', dest='dpi', type=int, default=0,
        help='Resolution of the exported images. Set to 0 for default.')
    parser.add_argument(
        '-d', '--device', dest='opencl_device', default='', type=str,
        help='OpenCL device to use for the simulations.')
    parser.add_argument(
        '-i', '--index', dest='opencl_index', default=0, type=int,
        help='OpenCL device index to use for the simulations.')

    parser.add_argument(
        '--inclusion-n', dest='inclusion_n', default=-1.0, type=float,
        help='Refractive index of the cylindrical inclusion.')

    args = parser.parse_args()
    filename = args.filename
    if not filename:
        filename = None
    num_packets = max(1, args.num_packets)
    batch_packets = max(1, args.batch_packets)
    fps = args.fps
    if fps <= 0.0:
        fps = None
    duration = args.duration
    if duration <= 0.0:
        duration = None
    title = args.title
    axis_off = args.axis_off
    inclusion_n = args.inclusion_n
    if inclusion_n < 1.0:
        inclusion_n = None
    dpi = args.dpi
    if dpi <= 0:
        dpi = None

    return {'filename': filename,
            'num_packets': num_packets,
            'batch_packets': batch_packets,
            'overwrite': args.overwrite,
            'verbose': args.verbose,
            'animation': args.animation,
            'fps': fps,
            'duration': duration,
            'autoscale': not args.no_autoscale,
            'title': title,
            'axis_off': axis_off,
            'opencl_device': args.opencl_device,
            'opencl_index': args.opencl_index,
            'logscale': args.logscale,
            'dpi': dpi,
            'inclusion_n': inclusion_n,}, args

def main(sim_task: callable, ani_task: callable, config: dict):
    '''
    Code for the __main__ section of the sampling volume simulation modules.

    Parameters
    ----------
    sim_task: callable
        A function that runs the sampling volume simulations.
    ani_task: callable
        A function that runs/prepares the animations from the simulation
        results.
    config: dict
        Sampling volume simulation configuration as a dict. Note that the
        entries of this configuration dict can change during the call.
    '''
    kwargs, parsed_args = mc_tr_deposition_cli(
        'Time-resolved deposition for a medium with a cylindrical inclusion.')

    if kwargs['verbose']:
        print('Arguments returned by the command line interface:')
        print(kwargs)

    if kwargs['batch_packets'] is not None:
        config['batch_packets'] = kwargs['batch_packets']
    if kwargs['num_packets'] is not None:
        config['num_packets'] = kwargs['num_packets']

    if kwargs['inclusion_n'] is not None:
        config['inclusion']['n'] = kwargs['inclusion_n']

    if kwargs['opencl_device']:
        cl_device = clinfo.device(
            kwargs['opencl_device'], index=kwargs['opencl_index'])
    else:
        cl_device = clinfo.gpus()[kwargs['opencl_index']]

    try:
        result = sim_task(
            kwargs['filename'], config, cl_device=cl_device,
            overwrite=kwargs['overwrite'], verbose=kwargs['verbose'])
    except FileExistsError:
        if kwargs['verbose']:
            print('Existing file with simulation results found!')
            print('Skipping simulations!')
        pass

    kwargs['animation'] = True
    if kwargs['animation']:
        filename = kwargs['filename']
        if kwargs['filename'] is None:
            filename = result
        ani_task(
            filename, overwrite=kwargs['overwrite'],
            fps=kwargs['fps'], duration=kwargs['duration'],
            title=kwargs['title'], axis_off=kwargs['axis_off'],
            logscale=kwargs['logscale'], autoscale=kwargs['autoscale'],
            verbose=kwargs['verbose'])
