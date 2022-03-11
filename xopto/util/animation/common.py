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

from typing import Tuple
import os.path
import tempfile

import numpy as np
import matplotlib.pyplot as pp
import matplotlib as mpl
from matplotlib.animation import FuncAnimation


def create_frame_animation(frames: np.ndarray, filename: str = None,
                           overwrite=False, logscale=False,
                           fps: float = None, duration: float = None,
                           xlabel=None, ylabel=None, title: str = None,
                           axis_off: bool = False, tight_layout: bool = False,
                           cbar: bool = False, cbar_label : str = None,
                           cbar_tick_format: str = None, cmap: str = None,
                           extent: Tuple[float, float, float, float] = None,
                           xrange: Tuple[float, float] = None,
                           yrange: Tuple[float, float] = None,
                           autoscale: bool = True,
                           imshow_kwargs: dict = None,
                           writer: str or mpl.animation.MovieWriter = None,
                           dpi: int = None,
                           verbose: bool = False) -> FuncAnimation:
    '''
    Creates a continuous animation of frames.

    Parameters
    ----------
    data: np.ndarray
        A 3D array of frames organized as [frame, y, x].
    filename: str
        Output file for the animation.
    overwrite: bool
        If False, existing animation files will not be overwritten.
    logscale: bool
        If True, logarithmically scale the frame data.
    fps: float
        Frame rate of the animation (1/s). Another way to control the frame
        rate is specify the duration parameter.
    duration: float
        Duration of one full animation (s). The same effect can be achieved
        with the fps parameter.
    title: str
        Plot title.
    xlabel: str
        Plot x axis label.
    ylabel: str
        Plot y axis label.
    axis_off: bool 
        Turns off the axis and related labels.
    tight_layout: bool
        Apply tight layout to the plot.
    cbar: bool
        Show a colorbar.
    cbar_label: str
        Text displayed by the colorbar label.
    cbar_tick_format: str
        Colorbar tick format string.
    cmap: str
        Colormap of the image.
    extent: Tuple[float, float, float]
        Size of the frame as [xmin, xmax, ymin, ymax].
    xrange: Tuple[float, float]
        Range of the x axis as (xmin, xmax).
    yrange: Tuple[float, float]
        Range of the y axis as (ymin, ymax).
    autoscale: bool
        Independently autoscale the intensity of each frame.
    writer: str or mpl.animation.MovieWriter
        Movie writer.
    imshow_kwargs: dict
        Additional keyword arguments for pyplot.imshow.
    dpi: int
        Resolution of the exported images.
    verbose: bool
        Turns on verbose progress report.

    Returns
    -------
    animation: FuncAnimation
        Animation object.
    '''
    if cbar_tick_format is None:
        cbar_tick_format = '{}'

    if imshow_kwargs is None:
        imshow_kwargs = {}

    frames = np.asarray(frames)
    if frames.ndim != 3:
        raise ValueError('The frames data array must have exactly 3 dims!')

    if filename is not None and os.path.isfile(filename) and not overwrite:
        if verbose:
            print('Animation file "{}" already exists!'.format(filename))
            print('Animation will not be saved!')
        filename = None

    if logscale:
        l = frames[frames > 0.0].min()
        frames[frames <= 0.0] = l
        frames = np.log10(frames)

    num_frames = frames.shape[0]
    imshow_kwargs.update({'cmap': cmap, 'extent': extent})
    vmin, vmax = frames.min(), frames.max()
    if not autoscale:
        imshow_kwargs.update({'vmin': vmin, 'vmax': vmax})

    fig, ax = pp.subplots(dpi=dpi)
    img = ax.imshow(frames[0], **imshow_kwargs)
    if cbar:
        colorbar = pp.colorbar(img)
        if logscale:
            pp.savefig(tempfile.TemporaryFile(), format='png')
            ticks = colorbar.ax.get_yticklabels()
            new_ticks = []
            for tick in ticks:
                if isinstance(tick, mpl.text.Text):
                    value = tick.get_text()
                    if value:
                        value = float(value)
                elif isinstance(tick, (int, float)):
                    value = float(tick)
                if isinstance(value, (int, float)):
                    new_ticks.append('$10^{' + cbar_tick_format.format(value) + '}$')
            colorbar.ax.set_yticklabels(new_ticks)
        if cbar_label:
            colorbar.set_label(cbar_label)

    if fps is None and duration is not None:
        fps = num_frames/duration

    if fps is None:
        fps = 25

    fargs = [{'frame_index': 0}]

    def ani_init():
        if xrange is not None:
            ax.set_xlim(xrange[0], xrange[1])
        if yrange is not None:
            ax.set_ylim(yrange[0], yrange[1])
        if not axis_off:
            if xlabel is not None:
                ax.set_xlabel(xlabel)
            if ylabel is not None:
                ax.set_ylabel(ylabel)
            if title is not None:
                if isinstance(title, str):
                    ax.set_title(title)
                else:
                    ax.set_title(title[0])
        else:
            ax.set_axis_off()

        if tight_layout:
            pp.tight_layout()

        fargs[0]['frame_index'] = 0

        return img,

    def ani_update(frame, *fargs):
        fargs[0]['frame_index'] += 1
        if verbose:
            print('Processing frame {}/{}'.format(
                fargs[0]['frame_index'], num_frames))
        img.set_array(frame)
        if autoscale:
            vmin, vmax = frame.min(), frame.max()
            img.set_clim(vmin, vmax)
        if isinstance(title, (tuple, list)):
            ax.set_title(title[fargs[0]['frame_index'] % len(title)])
        return img,

    ani = FuncAnimation(fig, ani_update, frames, fargs=fargs,
                        init_func=ani_init, blit=False)

    if filename is not None:
        ani.save(filename, writer=writer, fps=fps)

    pp.show()

    return ani


def create_path_animation(x: np.ndarray, y: np.ndarray, filename: str = None,
                          overwrite=False, logscale=False,
                          fps: float = None, duration: float = None,
                          xlabel=None, ylabel=None, title: str = None,
                          axis_off: bool = False, tight_layout: bool = False,
                          xrange: Tuple[float, float] = None,
                          yrange: Tuple[float, float] = None,
                          plot_kwargs: dict = None,
                          writer: str or mpl.animation.MovieWriter = None,
                          verbose: bool = False) -> FuncAnimation:
    '''
    Creates a continuous animation of the paths specified in the x and y
    arrays.

    Parameters
    ----------
    x: np.ndarray
        X coordinates of the path as a vector or a 2D array of shape
        (num_paths, path_length) or (1, path_length) or (path_length,).
    y: np.ndarray
        Y coordinates of the path as a vector or a 2D array of shape
        (num_paths, path_length) or (path_length,).
    filename: str
        Output file for the animation.
    overwrite: bool
        If False, existing animation files will not be overwritten.
    logscale: bool
        If True, logarithmically scale the data.
    fps: float
        Frame rate of the animation (1/s). Another way to control the frame
        rate is specify the duration parameter.
    duration: float
        Duration of one full animation (s). The same effect can be achieved
        with the fps parameter.
    title: str
        Plot title.
    xlabel: str
        Plot x axis label.
    ylabel: str
        Plot y axis label.
    axis_off: bool 
        Turns off the axis and related labels.
    tight_layout: bool
        Apply tight layout to the plot.
    xrange: Tuple[float, float]
        Range of the x axis as (xmin, xmax).
    yrange: Tuple[float, float]
        Range of the y axis as (ymin, ymax).
    writer: str or mpl.animation.MovieWriter
        Movie writer.
    plot_kwargs: dict
        Additional keyword arguments for pyplot.plot.
    verbose: bool
        Turns on verbose progress report.

    Returns
    -------
    animation: FuncAnimation
        Animation object.
    '''

    if filename is not None and os.path.isfile(filename) and not overwrite:
        if verbose:
            print('Animation file "{}" already exists!'.format(filename))
            print('Animation will not be saved!')
        filename = None

    if plot_kwargs is None:
        plot_kwargs = {}

    x = np.asarray(x)
    y = np.asarray(y)

    if logscale:
        l = x[x > 0.0].min()
        x[x <= 0.0] = l
        x = np.log10(x)

        l = y[y > 0.0].min()
        y[y <= 0.0] = l
        y = np.log10(y)

    if y.ndim == 1:
        y.shape = (1, y.size)
    if x.ndim == 1:
        x.shape = (1, x.size)

    if x.shape[1] != y.shape[1]:
        raise ValueError('The number of points in the x and y '
                         'arrays must be equal!')

    if x.shape[0] > 1 and x.shape[0] != y.shape[0]:
        raise ValueError('The number of paths in the x and y arrays '
                         'must be equal!')

    fig, ax = pp.subplots()
    plots = []

    num_paths = y.shape[0]
    num_steps = y.shape[1]

    if fps is None and duration is not None:
        fps = num_frames/duration

    if fps is None:
        fps = 25

    def ani_init():
        if xrange:
            ax.set_xlim(xrange[0], xrange[1])
        if yrange:
            ax.set_ylim(yrange[0], yrange[1])

        if not axis_off:
            if xlabel is not None:
                ax.set_xlabel(xlabel)
            if ylabel is not None:
                ax.set_ylabel(ylabel)
            if title is not None:
                if isinstance(title, str):
                    ax.set_title(title)
                else:
                    ax.set_title(title[0])
        else:
            ax.set_axis_off()

        if tight_layout:
            pp.tight_layout()

        plots.clear()
        for _ in range(num_paths):
            plots.append(ax.plot([], [], **plot_kwargs)[0])

        return plots

    def ani_update(step):
        for path_index in range(num_paths):
            if x.shape[0] > 1:
                x_data = x[path_index, :step + 1]
            else:
                x_data = x[0, :step + 1]
            y_data = y[path_index, :step + 1]
            
            plots[path_index].set_data(x_data, y_data)

            if isinstance(title, (tuple, list)):
                ax.set_title(title[step % len(title)])
        return plots

    ani = FuncAnimation(fig, ani_update, range(num_steps - 1),
                        init_func=ani_init, blit=False)
    if filename is not None:
        ani.save(filename, writer=writer, fps=fps)
    pp.show()

    return ani


if __name__ == '__main__':
    num_frames = 25
    H = W = 100
    R = 5

    x = np.arange(-W/2.0, W/2.0, dtype=np.float)
    y = np.arange(-H/2.0, H/2.0, dtype=np.float)
    Y, X = np.meshgrid(y, x, indexing='ij')
    frames = np.ones([num_frames, H, W])
    c_x, c_y = -W*0.5, -H*0.5
    for frame_index in range(num_frames):
        mask = (X - c_x)**2 + (Y - c_y)**2 < R**2
        frames[frame_index][mask] = 100
        c_x += W/num_frames
        c_y += H/num_frames

    create_frame_animation(
        frames, xlabel='$x$ (mm)', ylabel='$y$ (mm)',
        xrange=(x[0], x[-1]), yrange=(y[0], y[-1]),
        extent=(x[0], x[-1], y[0], y[-1]), logscale=True,
        cbar=True, cbar_label='Intensity (a.u.)',
        tight_layout=True, writer='imagemagick'
    )
