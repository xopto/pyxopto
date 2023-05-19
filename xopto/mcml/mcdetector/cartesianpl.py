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

import numpy as np
import scipy.constants
from xopto.mcml.mcutil import interp

from xopto.mcml.mcdetector.base import Detector
from xopto.mcml.mcdetector.radialpl import RadialPl
from xopto.mcml import cltypes, mcobject, mcoptions
from xopto.mcml.mcutil import axis


class CartesianPl(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClCartesian(cltypes.Structure):
            '''
            Structure that that represents a Cartesian detector
            in the Monte Carlo simulator core.

            Fields
            ------
            direction: mc_point3f_t
                Reference direction of the detector.
            x_min: mc_fp_t
                The leftmost edge of the first accumulator along the x axis.
            inv_dx: mc_fp_t
                Inverse of the spacing between the accumulators along the
                x axis.
            y_min: mc_fp_t
                The leftmost edge of the first accumulator along the
                y axis.
            inv_dy: mc_fp_t
                Inverse of the spacing between the accumulators along the
                y axis.
            pl_min: mc_fp_t
                The leftmost edge of the first optical path length accumulator.
            inv_dpl: mc_fp_t
                Inverse of the width of the optical path length accumulators.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle relative to the
                reference direction.
            n_x: mc_size_t
                The number of accumulators along the x axis.
            n_y: mc_size_t
                The number of accumulators along the y axis.
            n_pl: mc_size_t
                The number of path length accumulators.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            pl_log_scale: mc_int_t
                A flag indicating logarithmic scale of the path length axis.
            '''
            _pack_ = 1
            _fields_ = [
                ('direction', T.mc_point3f_t),
                ('x_min', T.mc_fp_t),
                ('inv_dx', T.mc_fp_t),
                ('y_min', T.mc_fp_t),
                ('inv_dy', T.mc_fp_t),
                ('pl_min', T.mc_fp_t),
                ('inv_dpl', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n_x', T.mc_size_t),
                ('n_y', T.mc_size_t),
                ('n_pl', T.mc_size_t),
                ('offset', T.mc_size_t),
                ('pl_log_scale', T.mc_int_t),
            ]
        return ClCartesian

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the accumulator in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_point3f_t direction;',
            '	mc_fp_t x_min;',
            '	mc_fp_t inv_dx;',
            '	mc_fp_t y_min;',
            '	mc_fp_t inv_dy;',
            '	mc_fp_t pl_min;',
            '	mc_fp_t inv_dpl;',
            '	mc_fp_t cos_min;',
            '	mc_size_t n_x;',
            '	mc_size_t n_y;',
            '	mc_size_t n_pl;',
            '	mc_size_t offset;',
            '	mc_int_t pl_log_scale;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the accumulator in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_detector('.format(loc),
            '		__mc_detector_mem const Mc{}Detector *detector){{'.format(Loc),
            '	dbg_print("Mc{}Detector - CartesianPl detector:");'.format(Loc),
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_float(INDENT "x_min (mm):", detector->x_min*1e3f);',
            '	dbg_print_float(INDENT "inv_dx (1/mm):", detector->inv_dx*1e-3f);',
            '	dbg_print_float(INDENT "y_min (mm):", detector->y_min*1e3f);',
            '	dbg_print_float(INDENT "inv_dy (1/mm):", detector->inv_dy*1e-3f);',
            '	dbg_print_float(INDENT "pl_min (um):", detector->pl_min*1e6f);',
            '	dbg_print_float(INDENT "inv_dpl (1/um)", detector->inv_dpl*1e-6f);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_float(INDENT "n_x:", detector->n_x);',
            '	dbg_print_float(INDENT "n_y:", detector->n_y);',
            '	dbg_print_float(INDENT "n_pl:", detector->n_pl);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '	dbg_print_int(INDENT "pl_log_scale:", detector->pl_log_scale);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim,',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir,',
            '		mc_fp_t weight){',
            '',
            '	__mc_detector_mem const struct Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '	__global mc_accu_t *address;',
            '	mc_int_t index_x, index_y, index_pl;',
            '	mc_size_t index;'
            '',
            '	dbg_print_status(mcsim, "{} Cartesian detector hit");'.format(Loc),
            '',
            '	index_x = mc_int((pos->x - detector->x_min)*detector->inv_dx);',
            '	index_x = mc_clip(index_x, 0, detector->n_x - 1);',
            '',
            '	index_y = mc_int((pos->y - detector->y_min)*detector->inv_dy);',
            '	index_y = mc_clip(index_y, 0, detector->n_y - 1);',
            '',
            '	mc_fp_t pl = mcsim_optical_pathlength(mcsim);',
            '	if (detector->pl_log_scale)',
            '	    pl = mc_log(mc_fmax(pl, FP_PLMIN));',
            '	index_pl = mc_int((pl - detector->pl_min)*detector->inv_dpl);',
            '	index_pl = mc_clip(index_pl, 0, detector->n_pl - 1);',
            '',
            '	index = (index_pl*detector->n_y + index_y)*detector->n_x + index_x;',
            '',
            '	address = mcsim_accumulator_buffer_ex(',
            '		mcsim, index + detector->offset);',
            '',
            '	mc_point3f_t detector_direction = detector->direction;',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(mc_dot_point3f(dir, &detector_direction)));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} Cartesian detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def cl_options(self, mc, target=None) -> mcoptions.RawOptions:
        '''
        OpenCL kernel options defined by this object.
        '''
        return [('MC_TRACK_OPTICAL_PATHLENGTH', True)]

    def __init__(self, xaxis: axis.Axis or 'CartesianPl', yaxis: axis.Axis = None,
                 plaxis: axis.Axis = None, cosmin=0.0,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        2D Cartesian reflectance/transmittance accumulator in the x-y plane.

        The grid of the Cartesian accumulators corresponds to a 3D numpy array
        with the first dimension representing the optical path length axis,
        the second dimension representing the y axis and the third dimension
        representing the x axis (reflectance[pl, y, x] or
        transmittance[pl, y, x]).

        Parameters
        ----------
        xaxis: axis.Axis or CartesianPl
            Object that defines accumulators along the x axis.
        yaxis: axis.Axis
            Object that defines accumulators along the y axis. If None, the
            y axis will equal x axis.
        plaxis: axis.Axis
            Object that defines the accumulators along the optical path length
            axis (this axis supports log-scale).
        cosmin: float
            Cosine of the maximum acceptance angle (relative to the direction)
            of the detector.
         direction: (float, float, float)
            Reference direction/orientation of the source.

        Note
        ----
        The first dimension of the accumulator represents the optical
        path length axis, the second dimension represents the y axis and the
        third dimension represents the x axis. 
        '''
        if isinstance(xaxis, CartesianPl):
            detector = xaxis
            xaxis = type(detector.xaxis)(detector.xaxis)
            yaxis = type(detector.yaxis)(detector.yaxis)
            plaxis = type(detector.plaxis)(detector.plaxis)
            cosmin = detector.cosmin
            direction = detector.direction
            raw_data = np.copy(detector.raw)
            nphotons = detector.nphotons
        else:
            if yaxis is None:
                yaxis = axis.Axis(xaxis)
            if plaxis is None:
                plaxis = axis.Axis(0.0, 1.0, 1)

            raw_data = np.zeros((plaxis.n, yaxis.n, xaxis.n))
            nphotons = 0

        super().__init__(raw_data, nphotons)
        self._cosmin = 0.0
        self._direction = np.zeros((3,))

        self._x_axis = xaxis
        self._y_axis = yaxis
        self._pl_axis = plaxis

        self._set_cosmin(cosmin)
        self._set_direction(direction)

        self._accumulators_area = self._x_axis.step*self._y_axis.step

    def _get_x_axis(self) -> axis.Axis:
        return self._x_axis
    xaxis = property(_get_x_axis, None, None, 'Axis object of the x axis.')

    def _get_y_axis(self) -> axis.Axis:
        return self._y_axis
    yaxis = property(_get_y_axis, None, None, 'Axis object of the y axis.')

    def _get_plaxis(self) -> axis.Axis:
        return self._pl_axis
    plaxis = property(_get_plaxis, None, None, 'Path length axis object.')

    def _get_taxis(self) -> axis.Axis:
        return axis.Axis(
            self._pl_axis.start/scipy.constants.c,
            self._pl_axis.stop/scipy.constants.c,
            self._pl_axis.n, logscale=self._pl_axis.logscale)
    taxis = property(_get_taxis, None, None, 'Time axis (s) derived from '
                                             'the path length axis object.')

    def _get_cosmin(self) -> Tuple[float, float]:
        return self._cosmin
    def _set_cosmin(self, value: float or Tuple[float, float]):
        self._cosmin = min(max(float(value), 0.0), 1.0)
    cosmin = property(_get_cosmin, _set_cosmin, None,
                      'Cosine of the maximum acceptance angle.')

    def _get_direction(self) -> Tuple[float, float, float]:
        return self._direction
    def _set_direction(self, direction: Tuple[float, float, float]):
        self._direction[:] = direction
        norm = np.linalg.norm(self._direction)
        if norm == 0.0:
            raise ValueError('Direction vector norm/length must not be 0!')
        self._direction *= 1.0/norm
    direction = property(_get_direction, _set_direction, None,
                        'Detector reference direction.')

    def _get_x(self) -> np.ndarray:
        return self._x_axis.centers
    x = property(_get_x, None, None,
                 'Centers of the accumulators along the x axis.')

    def _get_y(self) -> np.ndarray:
        return self._y_axis.centers
    y = property(_get_y, None, None,
                 'Centers of the accumulators along the y axis.')

    def _get_x_edges(self) -> np.ndarray:
        return self._x_axis.edges
    xedges = property(_get_x_edges, None, None,
                      'Edges of the accumulators along the x axis.')

    def _get_y_edges(self) -> np.ndarray:
        return self._y_axis.edges
    yedges = property(_get_y_edges, None, None,
                      'Edges of the accumulators along the y axis.')

    def _get_nx(self) -> int:
        return self._x_axis.n
    nx = property(_get_nx, None, None,
                  'Number of accumulators along the x axis.')

    def _get_ny(self) -> int:
        return self._y_axis.n
    ny = property(_get_ny, None, None,
                  'Number of accumulators along the y axis.')

    def _get_pl(self):
        return self._pl_axis.centers
    pl = property(_get_pl, None, None,
                  'Centers of the optical path length axis accumulators.')

    def _get_pledges(self):
        return self._pl_axis.edges
    pledges = property(_get_pledges, None, None,
                       'Edges of the optical path length axis accumulators.')

    def _get_npl(self):
        return self._pl_axis.n
    npl = property(_get_npl, None, None,
                   'Number of accumulators in the optical path length axis.')

    def _get_t(self):
        return self._pl_axis.centers*(1.0/scipy.constants.c)
    t = property(_get_t, None, None,
                  'Centers of the optical path length axis accumulators '
                  'expressed in propagation time (s).')

    nt = property(_get_npl, None, None,
                  'Number of accumulators in the time axis derived from the '
                  'optical path length axis.')

    def _get_tedges(self):
        return self._pl_axis.edges*(1.0/scipy.constants.c)
    tedges = property(_get_tedges, None, None,
                      'Edges (s) of the time axis accumulators derived from '
                      'the optical path length axis.')

    def meshgrid(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Returns 3D arrays of x, y and path length coordinates of the centers
        of accumulators that match the size of the reflectance / transmittance
        arrays.
        The grid of the Cartesian accumulators corresponds to a 2D numpy array
        with the first dimension representing the y axis, the second dimension
        representing the x axis (reflectance[y, x] or transmittance[y, x])
        and the third dimension the optical path length axis.

        Returns
        -------
        pl: np.ndarray
            A 3D array of optical path length coordinates.
        x: np.ndarray
            A 3D array of x coordinates.
        y: np.ndarray
            A 3D array of y coordinates.
        '''
        PL, Y, X = np.meshgrid(self.pl, self.y, self.x, indexing='ij')
        return X, Y, PL

    def _get_normalized(self) -> np.ndarray:
        return self.raw*(1.0/(max(self.nphotons, 1.0)*self._accumulators_area))
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def radial(self, raxis: axis.Axis, order=2, nfi=360, center=(0.0, 0.0)) \
            -> Detector:
        '''
        Converts this x-y accumulator to a radial (Radial) accumulator.

        Parameters
        ----------
        raxis: axis.Axis
            Radial axis of the accumulator.
        order: int
            Order of interpolation used in the conversion process.
        nfi: int
            Angular discretization used during conversion.
        center: (float, float)
            Center of transformation to radial coordinates as a tuple
            (x_center, y_center).

        Returns
        -------
        radial: xopto.mcml.mcdetector.base.Detector
            RadialPl representation of the Cartesian detector.
        '''
        plaxis = axis.Axis(self.plaxis)
        nfi = int(nfi)
        dfi = 2.0*np.pi/nfi
        fi = np.linspace(0.5*dfi, 2.0*np.pi - 0.5*dfi, nfi)

        a = np.pi*(raxis.edges[1:]**2 - raxis.edges[:-1]**2)
        w = a/(nfi*self._accumulators_area)
        w.shape = (1, raxis.n)

        R, Fi = np.meshgrid(raxis.centers, fi)
        Xi, Yi = R*np.cos(Fi) + center[0], R*np.sin(Fi) + center[1]

        raw = np.zeros(self.npl, raxis.n)
        for pl_index, pl in enumerate(self.pl):
            raw_slice = interp.interp2(Xi, Yi, self.x,
                                       self.y, self.raw[pl_index], order=order)
            raw_slice = (raw_slice*w).sum(0)
            raw[pl_index] = raw_slice
    
        radialpl = RadialPl(raxis, plaxis, self.cosmin)
        radialpl.set_raw_data(raw, self.nphotons)

        return radialpl

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`Cartesian.cl_type` method for a detailed
        list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
            Filled ctypes structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        allocation = mc.cl_allocate_rw_accumulator_buffer(self, self.shape)
        target.offset = allocation.offset

        target.direction.fromarray(self._direction)

        target.x_min = self._x_axis.start
        if self._x_axis.n > 1:
            target.inv_dx = 1.0/self._x_axis.step
        else:
            target.inv_dx = 0.0
        target.y_min = self._y_axis.start
        if self._y_axis.n > 1:
            target.inv_dy = 1.0/self._y_axis.step
        else:
            target.inv_dy = 0.0
        target.cos_min = self.cosmin
        target.n_x = self._x_axis.n
        target.n_y = self._y_axis.n

        target.pl_min = self._pl_axis.scaled_start
        if self._pl_axis.step != 0.0:
            target.inv_dpl = 1.0/self._pl_axis.step
        else:
            target.inv_dpl = 0.0
        target.pl_log_scale = self._pl_axis.logscale
        target.n_pl = self._pl_axis.n


        return target

    def todict(self) -> dict:
        '''
        Save the accumulator configuration without the accumulator data to
        a dictionary. Use the :meth:`Cartesian.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'CartesianPl',
            'pl_axis': self._pl_axis.todict(),
            'x_axis': self._x_axis.todict(),
            'y_axis': self._y_axis.todict(),
            'cosmin': self._cosmin,
            'direction': self._direction.tolist()
        }

    @staticmethod
    def fromdict(data: dict) -> 'CartesianPl':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the py:meth:`CartesianPl.todict` method.
        '''
        data_ = dict(data)
        detector_type = data_.pop('type')
        if detector_type != 'CartesianPl':
            raise TypeError('Expected "CartesianPl" type but got "{}"!'.format(
                detector_type))

        xaxis_data = data_.pop('x_axis')
        xaxis_type = xaxis_data.pop('type')

        yaxis_data = data_.pop('y_axis')
        yaxis_type = yaxis_data.pop('type')

        plaxis_data = data_.pop('pl_axis')
        plaxis_type = plaxis_data.pop('type')

        return CartesianPl(
            getattr(axis, xaxis_type)(**xaxis_data),
            getattr(axis, yaxis_type)(**yaxis_data),
            getattr(axis, plaxis_type)(**plaxis_data),
            **data_
        )

    def plot(self, scale: str = 'log', axis: str ='z', autoscale: bool = True,
             show: bool = True, **kwargs):
        '''
        Show detector slices or integral projections.

        Parameters
        ----------
        scale: str
            Data scaling can be "log" for logarithmic or "lin" for linear.
        axis: str
            The axis of slicing ("pl", "t", "y" or "x") or ome of the
            projection planes ("xy", "xpl", "xt", "ypl" or "yt"). Note that
            the time axis (s) is derived from the path length axis.
        autoscale: bool
            Scale the color coding of individual slices to the corresponding
            range of weights. If True, the color coding changes from slice
            to slice.
        show: bool
            Show the plot window if True.
        **kwargs: dict
            Additional keyword arguments passed to pyplot.imshow.
        '''
        from xopto.util import sliceview

        aspect = 'auto'
        if 'aspect' in kwargs:
            aspect = kwargs.pop('aspect')

        data = self.reflectance

        if axis == 'plproj': axis = 'xy'
        if axis == 'yproj': axis = 'xpl'
        if axis == 'xproj': axis = 'ypl'

        # slicing or projection axis
        ax = {'pl': 0, 't': 0, 'y': 1, 'x': 2,
              'xy': 0, 'xpl': 1, 'xt': 1, 'ypl': 2, 'yt': 2}.get(axis)

        # time or path length axis
        if 'pl' in axis:
            pl_t_axis, pl_t_label = self._pl_axis, 'pl'
        else:
            pl_t_axis, pl_t_label = self.taxis, 't'

        if axis == 't':
            title_fmt_str = 'Slice {{slice}}/{} @ {} = {{pos:.6g}}'
        else:
            title_fmt_str = 'Slice {{slice}}/{} @ {} = {{pos:.6f}}'
        title = title_fmt_str.format(data.shape[ax], axis)

        logscale = scale == 'log'

        fig = None

        if ax == 0:
            extent = [self._x_axis.start, self._x_axis.stop,
                      self._y_axis.start, self._y_axis.stop]
            slices = pl_t_axis.centers
            xlabel, ylabel = 'x', 'y'
        elif ax == 1:
            extent = [self._x_axis.start, self._x_axis.stop,
                      pl_t_axis.start, pl_t_axis.stop]
            slices = self._y_axis.centers
            xlabel, ylabel = 'x', pl_t_label
        elif ax == 2:
            extent = [self._y_axis.start, self._y_axis.stop,
                      pl_t_axis.start, pl_t_axis.stop]
            slices = self._x_axis.centers
            xlabel, ylabel = 'y', pl_t_label

        window_title = 'CartesianPl SliceView'

        if axis in ('xy', 'xpl', 'xt', 'ypl' or 'yt'):
            import matplotlib.pyplot as pp

            fig = pp.figure()
            data_slice = data.sum(axis=ax)
            low = data_slice.min()
            if scale == 'log':
                if low < 0:
                    data_slice = np.log(data_slice + (1.0 - low))
                else:
                    data_slice = np.log(data_slice + 1.0)

            pp.imshow(data_slice, extent=extent, origin='lower',
                      aspect=aspect, **kwargs)
            pp.xlabel(xlabel)
            pp.ylabel(ylabel)
            pp.title('Integral projection along the {:s} axis'.format(
                {'xy': 'pl', 'xpl': 'y',
                 'xt': 'y', 'ypl': 'x', 'yt': 'x'}.get(axis)))
            fig.canvas.manager.set_window_title(window_title)
            pp.tight_layout()
            if show:
                pp.show()

        else:
            sv = sliceview.SliceView(
                data, axis=ax, slices=slices, title=title, logscale=logscale,
                extent=extent, xlabel=xlabel, ylabel=ylabel,
                autoscale=autoscale, aspect=aspect, **kwargs)
            sv.fig.canvas.manager.set_window_title(window_title)
            if show:
                sv.show()

    def __str__(self):
        return 'CartesianPl(xaxis={}, yaxis={}, plaxis={}, cosmin={}, '\
               'direction=({}, {}, {}))'.format(
                   self._x_axis, self._y_axis, self._pl_axis,
                   self._cosmin, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
