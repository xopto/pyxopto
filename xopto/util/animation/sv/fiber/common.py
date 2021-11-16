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
from matplotlib.animation import FuncAnimation

import numpy as np

from xopto.util.animation.common import create_2d_animation
from xopto.mcml import mc
from xopto.cl import clinfo


def create_sv_animation_xz(data_source: str or dict, filename: str = None,
                           overwrite=False, logscale=False,
                           fps: float = None, duration: float = None,
                           axis_off: bool = False, title: str = None,
                           autoscale: bool = True, verbose: bool = False) \
                            -> FuncAnimation:
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
    verbose: bool
        Turns on verbose progress report.
    '''
    if isinstance(data_source, str):
        data = np.load(data_source, allow_pickle=True)
    elif isinstance(data_source, dict):
        data = data_source
    else:
        raise TypeError('Data source must be a dict or a filename!')

    if filename is None and isinstance(data_source, str):
        filename = os.path.splitext(data_source)[0] + '.gif'

    if filename is not None and os.path.isfile(filename) and not overwrite:
        if verbose:
            print('Animation file "{}" already exists!'.format(filename))
            print('Animation will not be saved!')
        filename = None

    sv_dict = data['sv']
    if not isinstance(sv_dict, dict):
        sv_dict = sv_dict.item()

    sv = mc.mcsv.SamplingVolume.fromdict(sv_dict)
    x_range = sv.xaxis.edges[0]*1e3, sv.xaxis.edges[-1]*1e3
    z_range = sv.zaxis.edges[-1]*1e3, sv.zaxis.edges[0]*1e3 

    extent = [x_range[0], x_range[1], z_range[0], z_range[1]]

    data = np.array(data['sv_data'])
    if data.ndim > 3:
        data = data.mean(2)

    frames = np.vstack([data, data[-2::-1]])

    return create_2d_animation(
        frames, filename, overwrite=overwrite, logscale=logscale,
        extent=extent, xrange=x_range, yrange=z_range, title=title,
        xlabel='$x$ (mm)', ylabel='$y$ (mm)', axis_off=axis_off,
        autoscale=autoscale, fps=fps, duration=duration)


def create_path_animation_xz(data_source: str or dict, filename: str = None,
                             overwrite=False, fps: float = None,
                             duration: float = None, logscale: bool = False,
                             axis_off: bool = False, title: str = None,
                             autoscale: bool = True, verbose: bool = False):
    '''
    Creates a continuous animation of the packet paths that are projected
    onto the x-z plane.

    Parameters
    ----------
    data_source: str or dict
        File with the simulation results or simulation results as returned by
        :py:func:`sv_fiber_incidence`.
    filename: str
        Output file for the animation. Use an empty string for no title.
    overwrite: bool
        If False, existing animation files will not be overwritten.
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
    logscale: bool
        Not used by this animation.
    autoscale: bool
        Not used by this animation.
    verbose: bool
        Turns on verbose progress report.
    '''
    if isinstance(data_source, str):
        data = np.load(data_source, allow_pickle=True)
    elif isinstance(data_source, dict):
        data = data_source
    else:
        raise TypeError('Data source must be a dict or a filename!')

    if filename is None and isinstance(data_source, str):
        filename = os.path.splitext(data_source)[0] + '.gif'

    if filename is not None and os.path.isfile(filename) and not overwrite:
        if verbose:
            print('Animation file "{}" already exists!'.format(filename))
            print('Animation will not be saved!')
        filename = None

    sv_dict = data['sv']
    if not isinstance(sv_dict, dict):
        sv_dict = sv_dict.item()
    sv = mc.mcsv.SamplingVolume.fromdict(sv_dict)
    x_range = sv.xaxis.edges[0], sv.xaxis.edges[-1]
    z_range = sv.zaxis.edges[-1], sv.zaxis.edges[0] 

    # extent = [x_range[0]*1e3, x_range[1]*1e3, z_range[0]*1e3, z_range[1]*1e3]

    trace_data = data['trace_data']
    trace_n= data['trace_n']
    num_packets, trace_len_max = trace_data.shape

    fig, ax = pp.subplots()
    plots = []

    num_frames = trace_n.max()

    # fill the remaining space in the trace with the terminal location
    for packet_index in range(num_packets):
        entries = trace_n[packet_index]
        for coordinate in ('x', 'y', 'z'):
            trace_data[packet_index][coordinate][entries:] = \
                trace_data[packet_index][coordinate][entries - 1]
        plots.append(ax.plot([], [])[0])

    if fps is None and duration is not None:
        fps = num_frames/duration

    if fps is None:
        fps = 0.125*num_frames

    def ani_init():
        ax.set_xlim(x_range[0]*1e3, x_range[1]*1e3)
        ax.set_ylim(z_range[0]*1e3, z_range[1]*1e3)
        if not axis_off:
            ax.set_xlabel('x (mm)')
            ax.set_ylabel('z (mm)')
            if title:
                ax.set_title(title)
        else:
            ax.set_axis_off()
        return plots

    def ani_update(path_segment):
        for packet_index in range(num_packets):
            x = trace_data[packet_index, :path_segment + 1]['x']
            z = trace_data[packet_index, :path_segment + 1]['z']
            plots[packet_index].set_data(x*1e3, z*1e3)
        return plots

    ani = FuncAnimation(fig, ani_update, range(num_frames - 1),
                        init_func=ani_init, blit=False)
    if filename is not None:
        ani.save(filename, writer='imagemagick', fps=fps)
    pp.show()

def mc_sv_fiber_cli(description: str = None,
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
        - autoscale: bool
            Independently autoscale the intensity of each frame.
        - opencl_device: str
            OpenCL device to use for the simulations.
        - opencl_index: int
            Index (zero-based) of the OpenCL device to use for the simulations.
    '''
    if description is None:
        description = 'Sampling volume simulation'

    parser = argparse.ArgumentParser(description='description')
    parser.add_argument(
        '-n', '--num-packets', dest='sv_num_packets', default=1000, type=int,
        help='Minimum number of packets reaching the detector.')
    parser.add_argument(
        '-p', '--projection', dest='sv_projection', action='store_true',
        help='Reduce the number of sampling volume accumulators '
             'along the y axis to 1.')
    parser.add_argument(
        '-b', '--batch-packets', dest='batch_packets', default=1000000, type=int,
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
        '-d', '--device', dest='opencl_device', default='', type=str,
        help='OpenCL device to use for the simulations.')
    parser.add_argument(
        '-i', '--index', dest='opencl_index', default=0, type=int,
        help='OpenCL device index to use for the simulations.')

    parser.add_argument(
        '--sample-g', dest='sample_g', default=-1.0, type=float,
        help='Anisotropy of the sample. Use a negative value of less than '
             '-1 for default.')
    parser.add_argument(
        '--sample-mua', dest='sample_mua', default=-1.0, type=float,
        help='Absorption coefficient of the sample (1/m). '
             'Use a negative value for default.')
    parser.add_argument(
        '--sample-mus', dest='sample_mus', default=-1.0, type=float,
        help='Scattering coefficient of the sample (1/m). '
             'Use a negative value for default.')
    parser.add_argument(
        '--sample-n', dest='sample_n', default=-1.0, type=float,
        help='Refractive index of the sample. '
             'Use a negative value for default.')
    parser.add_argument(
        '--fiber-sds', dest='fiber_sds', default=-1.0, type=float,
        help='Source detector separation (m).'
             'Use a negative value for default.')

    args = parser.parse_args()
    filename = args.filename
    if not filename:
        filename = None
    sv_num_packets = max(1, args.sv_num_packets)
    batch_packets = max(1, args.batch_packets)
    fps = args.fps
    if fps <= 0.0:
        fps = None
    duration = args.duration
    if duration <= 0.0:
        duration = None
    title = args.title
    axis_off = args.axis_off
    sample_mua = args.sample_mua
    if sample_mua < 0.0:
        sample_mua = None
    sample_mus = args.sample_mus
    if sample_mus < 0.0:
        sample_mus = None
    sample_g = args.sample_g
    if sample_g <= -1.0:
        sample_g = None
    sample_n = args.sample_n
    if sample_n <= 0.0:
        sample_n = None
    fiber_sds = args.fiber_sds
    if fiber_sds < 0.0:
        fiber_sds = None


    return {'filename': filename, 'batch_packets': batch_packets,
            'sv_num_packets': sv_num_packets, 'overwrite': args.overwrite,
            'sv_projection': args.sv_projection,
            'verbose': args.verbose, 'animation': args.animation,
            'fps': fps, 'duration': duration,
            'autoscale': not args.no_autoscale,
            'title': title, 'axis_off': axis_off,
            'opencl_device': args.opencl_device,
            'opencl_index': args.opencl_index, 'logscale': args.logscale,
            'sample_mua': sample_mua, 'sample_mus': sample_mus,
            'sample_g': sample_g, 'sample_n': sample_n,
            'fiber_sds': fiber_sds}, args

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
    kwargs, parsed_args = mc_sv_fiber_cli(
        'Sampling volume for optical fiber transmittance')

    if kwargs['verbose']:
        print('Arguments returned by the command line interface:')
        print(kwargs)

    if kwargs['batch_packets'] is not None:
        config['batch_packets'] = kwargs['batch_packets']
    if kwargs['sv_num_packets'] is not None:
        config['sv_num_packets'] = kwargs['sv_num_packets']
    if kwargs['sv_projection'] and 'sv' in config:
        config['sv']['yaxis']['n'] = 1
    if kwargs['sample_g'] is not None:
        config['sample']['g'] = kwargs['sample_g']
    if kwargs['sample_mua'] is not None:
        config['sample']['mua'] = kwargs['sample_mua']
    if kwargs['sample_mus'] is not None:
        config['sample']['mus'] = kwargs['sample_mus']
    if kwargs['sample_n'] is not None:
        config['sample']['n'] = kwargs['sample_n']
    if kwargs['fiber_sds'] is not None and 'fiber_sds' in config:
        if isinstance(config['fiber_sds'], (float, int)):
            config['fiber_sds'] = kwargs['fibers_sds']
        elif isinstance(config['fiber_sds'], (np.ndarray, tuple, list)):
            config['fiber_sds'] = np.linspace(
                config['fiber_sds'][0],
                kwargs['fiber_sds'],
                len(config['fiber_sds'])
            )

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
