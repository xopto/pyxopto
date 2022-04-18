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
from xopto.mcml.mcutil import interp

from xopto.mcml.mcdetector.base import Detector
from xopto.mcml.mcdetector.radial import Radial
from xopto.mcml import cltypes, mcobject
from xopto.mcml.mcutil import axis


class Cartesian(Detector):
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
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle relative to the
                reference direction.
            n_x: mc_size_t
                The number of accumulators along the x axis.
            n_y: mc_size_t
                The number of accumulators along the y axis.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            '''
            _pack_ = 1
            _fields_ = [
                ('direction', T.mc_point3f_t),
                ('x_min', T.mc_fp_t),
                ('inv_dx', T.mc_fp_t),
                ('y_min', T.mc_fp_t),
                ('inv_dy', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n_x', T.mc_size_t),
                ('n_y', T.mc_size_t),
                ('offset', T.mc_size_t)
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
            '	mc_fp_t cos_min;',
            '	mc_size_t n_x;',
            '	mc_size_t n_y;',
            '	mc_size_t offset;',
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
            '	dbg_print("Mc{}Detector - Cartesian detector:");'.format(Loc),
            '	dbg_print_point3f(INDENT "direction:", &detector->direction);',
            '	dbg_print_float(INDENT "x_min (mm):", detector->x_min*1e3f);',
            '	dbg_print_float(INDENT "inv_dx (1/mm):", detector->inv_dx*1e-3f);',
            '	dbg_print_float(INDENT "y_min (mm):", detector->y_min*1e3f);',
            '	dbg_print_float(INDENT "inv_dy (1/mm):", detector->inv_dy*1e-3f);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_float(INDENT "n_x:", detector->n_x);',
            '	dbg_print_float(INDENT "n_y:", detector->n_y);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
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
            '	mc_int_t index_x, index_y;',
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
            '	index = index_y*detector->n_x + index_x;',
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

    def __init__(self, xaxis, yaxis=None, cosmin=0.0,
                 direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)):
        '''
        2D Cartesian reflectance/transmittance accumulator in the x-y plane.

        The grid of the Cartesian accumulators corresponds to a 2D numpy array
        with the first dimension representing the y axis and second dimension
        representing the x axis (reflectance[y, x] or transmittance[y, x]).

        Parameters
        ----------
        xaxis: axis.Axis
            Object that defines accumulators along the x axis.
        yaxis: axis.Axis
            Object that defines accumulators along the y axis. If None, the
            y axis will equal x axis.
        cosmin: float
            Cosine of the maximum acceptance angle (relative to the direction)
            of the detector.
         direction: (float, float, float)
            Reference direction/orientation of the source.
        '''
        if isinstance(xaxis, Cartesian):
            detector = xaxis
            xaxis = type(detector.xaxis)(detector.xaxis)
            yaxis = type(detector.yaxis)(detector.yaxis)
            cosmin = detector.cosmin
            direction = detector.direction
            raw_data = np.copy(detector.raw)
            nphotons = detector.nphotons
        else:
            if yaxis is None:
                yaxis = axis.Axis(xaxis)

            raw_data = np.zeros((yaxis.n, xaxis.n))
            nphotons = 0

        super().__init__(raw_data, nphotons)
        self._cosmin = 0.0
        self._direction = np.zeros((3,))

        self._x_axis = xaxis
        self._y_axis = yaxis

        self._set_cosmin(cosmin)
        self._set_direction(direction)

        self._accumulators_area = self._x_axis.step*self._y_axis.step

    def _get_x_axis(self) -> axis.Axis:
        return self._x_axis
    xaxis = property(_get_x_axis, None, None, 'Axis object of the x axis.')

    def _get_y_axis(self) -> axis.Axis:
        return self._y_axis
    yaxis = property(_get_y_axis, None, None, 'Axis object of the y axis.')

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

    def meshgrid(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns 2D arrays of x and y coordinates of the centers of accumulators
        that match the size of the reflectance / transmittance arrays.
        The grid of the Cartesian accumulators corresponds to a 2D numpy array
        with the first dimension representing the y axis and second dimension
        representing the x axis (reflectance[y, x] or transmittance[y, x]).

        Returns
        -------
        x: np.ndarray
            A 2D array of x coordinates.
        y: np.ndarray
            A 2D array of y coordinates.
        '''
        Y, X = np.meshgrid(self.y, self.x, indexing='ij')
        return X, Y

    def _get_normalized(self) -> np.ndarray:
        return self.raw*(1.0/(max(self.nphotons, 1.0)*self._accumulators_area))
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def radial(self, raxis, order=2, nfi=360, center=(0.0, 0.0)) -> Detector:
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
            Radial representation of the Cartesian detector.
        '''
        nfi = int(nfi)
        dfi = 2.0*np.pi/nfi
        fi = np.linspace(0.5*dfi, 2.0*np.pi - 0.5*dfi, nfi)

        a = np.pi*(raxis.edges[1:]**2 - raxis.edges[:-1]**2)
        w = a/(nfi*self._accumulators_area)
        w.shape = (1, raxis.n)

        R, Fi = np.meshgrid(raxis.centers, fi)
        Xi, Yi = R*np.cos(Fi) + center[0], R*np.sin(Fi) + center[1]

        raw = interp.interp2(Xi, Yi, self.x, self.y, self.raw, order=order)
        raw = (raw*w).sum(0)

        radial = Radial(raxis, self.cosmin)
        radial.set_raw_data(raw, self.nphotons)

        return radial

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
            'type':'Cartesian',
            'x_axis': self._x_axis.todict(),
            'y_axis': self._y_axis.todict(),
            'cosmin': self._cosmin,
            'direction': self._direction.tolist()
        }

    @staticmethod
    def fromdict(data: dict) -> 'Cartesian':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the py:meth:`Cartesian.todict` method.
        '''
        data_ = dict(data)
        detector_type = data_.pop('type')
        if detector_type != 'Cartesian':
            raise TypeError('Expected "Cartesian" type but got "{}"!'.format(
                detector_type))
        xaxis_data = data_.pop('x_axis')
        xaxis_type = xaxis_data.pop('type')

        yaxis_data = data_.pop('y_axis')
        yaxis_type = yaxis_data.pop('type')

        return Cartesian(
            getattr(axis, xaxis_type)(**xaxis_data),
            getattr(axis, yaxis_type)(**yaxis_data),
            **data_
        )

    def plot(self, scale: str = 'log', raw: bool = False, show: bool = True):
        '''
        Show the detector contet as a 2D image.

        Parameters
        ----------
        scale: str
            Data scaling can be "log" for logarithmic or "lin" for linear.
        raw: bool
            Set to True to show the raw data. Default is False and shows the
            normalized (reflectance) content.
        show: bool 
        '''
        import matplotlib.pyplot as pp

        extent = [self._x_axis.start, self._x_axis.stop,
                  self._y_axis.start, self._y_axis.stop]

        data = self.raw if raw else self.reflectance
        which = 'raw' if raw else 'reflectance'
        
        if scale == 'log':
            mask = data > 0.0
            if mask.size > 0:
                log_data = np.tile(np.log10(data[mask].min()), data.shape)
                log_data[mask] = np.log10(data[mask])
                data = log_data

        fig, ax = pp.subplots()
        img = ax.imshow(data, extent=extent, origin='lower')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        pp.colorbar(img)

        fig.canvas.manager.set_window_title(
            'Cartesian detector - {} - {}'.format(scale, which))

        if show:
            pp.show()

    def __str__(self):
        return 'Cartesian(xaxis={}, yaxis={}, cosmin={}, '\
               'direction=({}, {}, {}))'.format(
                   self._x_axis, self._y_axis, self._cosmin, *self._direction)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
