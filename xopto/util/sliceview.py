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

from typing import Callable
import matplotlib.pyplot as pp
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
import numpy as np


class SliceView:
    def __init__(self, data: np.ndarray, axis: int = 0,
                 slices: np.ndarray = None,
                 logscale: bool or Callable[[np.ndarray], np.ndarray]=False,
                 autoscale: bool = False,
                 xlabel: str = None, ylabel: str = None, title: str = None,
                 **kwargs):
        '''
        Creates a matplotlib-based slice viewer of 3D data cubes along
        one axis.

        Parameters
        ----------
        data: np.ndarray
            A 3D Data array.
        along: int
            The axis along which to slice the data.
            Assuming that the dimensions of the data are (z, y, x), the
            produced slices are as follows:

            - along = 0; slicing along the z axis, slice dimensions (y, x)
            - along = 1; slicing along the y axis, slice dimensions (z, x)
            - along = 1; slicing along the x axis, slice dimensions (z, y)

            Note that the first dimension of the slice will be shown along the
            vertical/y axis of the slice image and the second dimension of the
            slice along the horizontal/x axis of the slice image.
        slices: np.ndarray
            A vector that defines the positions of slices. If None, a default
            sequence of positions from 0 to num_slices - 1 is created.
        logscale: bool or Callable[[np.ndarray], np.ndarray]
            A bool value enables or disables the logarithmic scaling of the
            data that is computed as:

            - :code:`np.log10(slice_data + 1)` if data are strictly positive
            - :code:`np.log10(slice_dat - data.min() + 1)` if data include negative values

            If a callable, all data are scaled with this function.
        autoscale: bool
            Automatically adjust the intensity range of the image to
            the current slice minimum and maximum intensity.
        title: str
            Slice title. This can be a format string with placeholders for
            the slice number "{slice}", total number of slices "{num_slices}"
            and position of the slice "{pos}".
        xlabel: str
            Label of the x axis.
        ylabel: str
            Label of the y axis.
        kwargs: dict
            Additional keyword arguments that are passed to pyplot.imshow.
        '''
        self._axis = int(axis)
        data = np.asarray(data)
        if title is None:
            title = 'Slice {slice}/{num_slices} @ {pos:.4f}'

        self._title = None
        self._title_fmtstr = title

        if axis != 0:
            order = {0:(0, 1, 2), 1:(1, 0, 2), 2:(2, 0, 1)}.get(axis)
            data = np.transpose(data, order)
        self._data = data
        self._data_range = (data.min(), data.max())

        self._logmin = np.finfo(float).eps

        if slices is None:
            slices = np.arange(data.shape[axis])
        self._slices = slices

        if callable(logscale):
            self._scale_fun = logscale
            logscale = True
        else:
            self._scale_fun = self._apply_logscale

        self._logscale = bool(logscale)

        self._autoscale = bool(autoscale)

        self._fig, self._ax_image = pp.subplots()
        self._fig.canvas.manager.set_window_title('Slice View')

        if 'vmin' in kwargs:
            kwargs.pop('vmin')
        if 'vmax' in kwargs:
            kwargs.pop('vmax')

        self._image = self._ax_image.imshow(
            data[0], vmin=self._data_range[0], vmax=self._data_range[1],
            **kwargs
        )
        pp.colorbar(self._image, orientation='vertical')

        pp.subplots_adjust(bottom=0.30, left=0)
        self._slice_ax = pp.axes([0.30, 0.16, 0.58, 0.03])
        self._vmin_ax =  pp.axes([0.30, 0.08, 0.58, 0.03])
        self._vmax_ax =  pp.axes([0.30, 0.03, 0.58, 0.03])
        self._scale_ax = pp.axes([0.04, 0.01, 0.15, 0.10])
        self._auto_ax =  pp.axes([0.04, 0.12, 0.15, 0.08])

        self._vmin_slider = Slider(
            self._vmin_ax, 'vmin',
            valinit=0.0, valmin=0.0, valmax=100.0, valfmt='%4.1f %%')
        self._vmin_slider.on_changed(self._vmin_changed)

        self._vmax_slider = Slider(
            self._vmax_ax, 'vmax',
            valinit=100.0, valmin=0.0, valmax=100.0, valfmt='%4.1f %%')
        self._vmax_slider.on_changed(self._vmax_changed)

        self._slice_slider = Slider(
            self._slice_ax, 'Slice', valinit=slices[0],
            valmin=slices[0], valmax=slices[-1], valfmt='%.4f')
        self._slice_slider.on_changed(self._slice_changed)

        self._scale_radio = RadioButtons(
            self._scale_ax, ('lin', 'log'), active=self._scale_active())
        self._scale_radio.on_clicked(self._set_scale)

        self._auto_button = CheckButtons(
            self._auto_ax, ['autoscale'], [self._autoscale])
        self._auto_button.on_clicked(self._autoscale_changed)

        if self._title_fmtstr is not None:
            self._title = self._ax_image.set_title(self._make_title())

        if xlabel is not None:
            self._ax_image.set_xlabel(xlabel)

        if ylabel is not None:
            self._ax_image.set_ylabel(ylabel)

        self._slice_changed(0)

    def _make_title(self, slice_index=0):
        return self._title_fmtstr.format(
            slice=slice_index + 1, num_slices=self._slices.size,
            pos=self._slices[slice_index]
        )

    def _apply_logscale(self, data: np.ndarray) -> np.ndarray:
        if self._data_range[0] < 0.0:
            return np.log10((1.0 - self._data_range[0]) + data)
        else:
            return np.log10(1.0 + data)

    def _slice_changed(self, pos):
        slice_index = np.argmin(np.abs(self._slices - pos))
        img_slice = self._data[slice_index]
        if self._logscale:
            img_slice = self._scale_fun(img_slice)
        self._image.set_array(img_slice)
        if self._autoscale:
            self._image.autoscale()
        if self._title is not None:
            self._title.set_text(self._make_title(slice_index))

        self._fig.canvas.draw_idle()

    def _set_scale(self, val):
        if val == 'log':
            self._logscale = True
        else:
            self._logscale = False
        if not self._autoscale:
            self._vmin_changed(self._vmin_slider.val)
            self._vmax_changed(self._vmax_slider.val)
        self._slice_changed(self._slice_slider.val)

    def _scale_active(self):
        return {False: 0, True:1}.get(self._logscale, 0)

    def _autoscale_changed(self, value):
        self._autoscale = not self._autoscale
        if not self._autoscale:
            self._vmin_changed(self._vmin_slider.val)
            self._vmax_changed(self._vmax_slider.val)
        self._slice_changed(self._slice_slider.val)

    def _vmin_changed(self, value):
        vmin = self._data_range[0] + \
            (self._data_range[1] - self._data_range[0])*(value*0.01)
        if self._logscale:
            vmin = self._scale_fun(vmin)
        self._image.set_clim(vmin=vmin)

    def _vmax_changed(self, value):
        vmax = self._data_range[0] + \
            (self._data_range[1] - self._data_range[0])*(value*0.01)
        if self._logscale:
            vmax = self._scale_fun(vmax)
        self._image.set_clim(vmax=vmax)

    def _get_ax(self):
        return self._ax_image

    axis = property(_get_ax, None, None, 'Image axis.')

    def _get_fig(self):
        return self._fig

    fig = property(_get_fig, None, None, 'Figure instance.')

    def show(self):
        '''
        Show the SliceView window.
        '''
        pp.show()

def show():
    pp.show()


class SliceViewCyl:
    def __init__(self, data: np.ndarray, axis: int = 0,
                 slices: np.ndarray = None,
                 logscale: bool or Callable[[np.ndarray], np.ndarray]=False,
                 autoscale: bool = False,
                 xlabel: str = None, ylabel: str = None, title: str = None,
                 **kwargs):
        '''
        Creates a matplotlib-based slice viewer of 3D data cubes along
        one axis.

        Parameters
        ----------
        data: np.ndarray
            A 3D Data array.
        along: int
            The axis along which to slice the data.
            Assuming that the dimensions of the data are (z, fi, r), the
            produced slices are as follows:

            - along = 0; slicing along the z axis, slice dimensions (fi, r)
            - along = 1; slicing along the fi axis, slice dimensions (z, r)
            - along = 1; slicing along the z axis, slice dimensions (z, fi)

        slices: np.ndarray
            A vector that defines the positions of slices. If None, a default
            sequence of positions from 0 to num_slices - 1 is created.
        logscale: bool or Callable[[np.ndarray], np.ndarray]
            A bool value enables or disables the logarithmic scaling of the
            data that is computed as:

            - :code:`np.log10(slice_data + 1)` if data are strictly positive
            - :code:`np.log10(slice_dat - data.min() + 1)` if data include negative values

            If a callable, all data are scaled with this function.
        autoscale: bool
            Automatically adjust the intensity range of the image to
            the current slice minimum and maximum intensity.
        title: str
            Slice title. This can be a format string with placeholders for
            the slice number "{slice}", total number of slices "{num_slices}"
            and position of the slice "{pos}".
        xlabel: str
            Label of the x axis.
        ylabel: str
            Label of the y axis.
        kwargs: dict
            Additional keyword arguments that are passed to pyplot.imshow or
            pyplot.pcolormesh.
        '''
        self._axis = int(axis)
        data = np.asarray(data)
        if title is None:
            title = 'Slice {slice}/{num_slices} @ {pos:.4f}'

        self._title = None
        self._title_fmtstr = title

        if axis != 0:
            order = {0:(0, 1, 2), 1:(1, 0, 2), 2:(2, 0, 1)}.get(axis)
            data = np.transpose(data, order)
        self._data = data
        self._data_range = (data.min(), data.max())

        self._logmin = np.finfo(float).eps

        if slices is None:
            slices = np.arange(data.shape[axis])
        self._slices = slices

        if callable(logscale):
            self._scale_fun = logscale
            logscale = True
        else:
            self._scale_fun = self._apply_logscale

        self._logscale = bool(logscale)

        self._autoscale = bool(autoscale)

        self._polar_axis = False
        subplot_kw = {}
        if self._axis == 0:
            self._polar_axis = True
            subplot_kw = {'projection': 'polar'}

        self._fig, self._ax_image = pp.subplots(subplot_kw=subplot_kw)
        self._fig.canvas.manager.set_window_title('Cylindrical slice View')

        if 'vmin' in kwargs:
            kwargs.pop('vmin')
        if 'vmax' in kwargs:
            kwargs.pop('vmax')

        if self._polar_axis:
            self._data = data = np.transpose(data, [0, 2, 1])
            R = Fi = None
            if 'R' in kwargs:
                R = kwargs.pop('R')
            if 'Fi' in kwargs:
                Fi = kwargs.pop('Fi')

            if R is not None and Fi is not None:
                self._image = pp.pcolormesh(
                    Fi, R, data[0],
                    vmin=self._data_range[0], vmax=self._data_range[1])
            else:
                fi = np.linspace(0, 2*np.pi, data.shape[2])
                r = np.arange(data.shape[1], dtype=np.float)
                R, Fi = np.meshgrid(r, fi, indexing='ij')

                self._image = pp.pcolormesh(
                    Fi, R, data[0],
                    vmin=self._data_range[0], vmax=self._data_range[1])
        else:
            self._image = self._ax_image.imshow(
                data[0], vmin=self._data_range[0], vmax=self._data_range[1],
                **kwargs
            )

        pp.colorbar(self._image, orientation='vertical')

        pp.subplots_adjust(bottom=0.30, left=0)
        self._slice_ax = pp.axes([0.30, 0.16, 0.58, 0.03])
        self._vmin_ax =  pp.axes([0.30, 0.08, 0.58, 0.03])
        self._vmax_ax =  pp.axes([0.30, 0.03, 0.58, 0.03])
        self._scale_ax = pp.axes([0.04, 0.01, 0.15, 0.10])
        self._auto_ax =  pp.axes([0.04, 0.12, 0.15, 0.08])

        self._vmin_slider = Slider(
            self._vmin_ax, 'vmin',
            valinit=0.0, valmin=0.0, valmax=100.0, valfmt='%4.1f %%')
        self._vmin_slider.on_changed(self._vmin_changed)

        self._vmax_slider = Slider(
            self._vmax_ax, 'vmax',
            valinit=100.0, valmin=0.0, valmax=100.0, valfmt='%4.1f %%')
        self._vmax_slider.on_changed(self._vmax_changed)

        self._slice_slider = Slider(
            self._slice_ax, 'Slice', valinit=slices[0],
            valmin=slices[0], valmax=slices[-1], valfmt='%.4f')
        self._slice_slider.on_changed(self._slice_changed)

        self._scale_radio = RadioButtons(
            self._scale_ax, ('lin', 'log'), active=self._scale_active())
        self._scale_radio.on_clicked(self._set_scale)

        self._auto_button = CheckButtons(
            self._auto_ax, ['autoscale'], [self._autoscale])
        self._auto_button.on_clicked(self._autoscale_changed)

        if self._title_fmtstr is not None:
            self._title = self._ax_image.set_title(self._make_title())

        if xlabel is not None:
            self._ax_image.set_xlabel(xlabel)

        if ylabel is not None:
            self._ax_image.set_ylabel(ylabel)

        self._slice_changed(0)

    def _make_title(self, slice_index=0):
        return self._title_fmtstr.format(
            slice=slice_index + 1, num_slices=self._slices.size,
            pos=self._slices[slice_index]
        )

    def _apply_logscale(self, data: np.ndarray) -> np.ndarray:
        if self._data_range[0] < 0.0:
            return np.log10((1.0 - self._data_range[0]) + data)
        else:
            return np.log10(1.0 + data)

    def _slice_changed(self, pos):
        slice_index = np.argmin(np.abs(self._slices - pos))
        img_slice = self._data[slice_index]
        if self._logscale:
            img_slice = self._scale_fun(img_slice)
        if self._polar_axis:
            self._image.set_array(img_slice.flatten())
        else:
            self._image.set_array(img_slice)
        if self._autoscale:
            self._image.autoscale()
        if self._title is not None:
            self._title.set_text(self._make_title(slice_index))

        self._fig.canvas.draw_idle()

    def _set_scale(self, val):
        if val == 'log':
            self._logscale = True
        else:
            self._logscale = False
        if not self._autoscale:
            self._vmin_changed(self._vmin_slider.val)
            self._vmax_changed(self._vmax_slider.val)
        self._slice_changed(self._slice_slider.val)

    def _scale_active(self):
        return {False: 0, True:1}.get(self._logscale, 0)

    def _autoscale_changed(self, value):
        self._autoscale = not self._autoscale
        if not self._autoscale:
            self._vmin_changed(self._vmin_slider.val)
            self._vmax_changed(self._vmax_slider.val)
        self._slice_changed(self._slice_slider.val)

    def _vmin_changed(self, value):
        vmin = self._data_range[0] + \
            (self._data_range[1] - self._data_range[0])*(value*0.01)
        if self._logscale:
            vmin = self._scale_fun(vmin)
        self._image.set_clim(vmin=vmin)

    def _vmax_changed(self, value):
        vmax = self._data_range[0] + \
            (self._data_range[1] - self._data_range[0])*(value*0.01)
        if self._logscale:
            vmax = self._scale_fun(vmax)
        self._image.set_clim(vmax=vmax)

    def _get_ax(self):
        return self._ax_image

    axis = property(_get_ax, None, None, 'Image axis.')

    def _get_fig(self):
        return self._fig

    fig = property(_get_fig, None, None, 'Figure instance.')

    def show(self):
        '''
        Show the SliceView window.
        '''
        pp.show()

def show():
    pp.show()

if __name__ == '__main__':
    data = np.random.randint(0, 255, (20,256,256)).astype(np.float)
    for i in range(data.shape[0]):
        data[i] = data[i]*(1.0 - float(i)/data.shape[0])
    sv = SliceView(data, slices=np.linspace(0, 1, 20), cmap='gray')

    pp.show()
