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

import time
from typing import Tuple, List, Dict

import numpy as np

from xopto.mcbase import cltypes
from xopto.mcbase import mctypes
from xopto.mcbase import mcoptions
from xopto.mcbase import mcobject
from xopto.mcbase.mcutil.axis import Axis

class SamplingVolume(mcobject.McObject):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClSamplingVolume(cltypes.Structure):
            '''
            OpenCL structure that is used by the Monte carlo simulator.

            Parameters
            ----------
            mc: McObject
                A Monte Carlo simulator instance.

            Returns
            -------
            struct: cltypes.Structure
                A structure type that represents the sampling volume in
                the Monte Carlo kernel.

            Fields
            ------
            top_left: McTypes.mc_point3f_t
                Coordinates of the top-left corner of the sampling volume accumulators.
            voxel_size: McTypes.mc_point3f_t
                Size of a voxel.
            shape: McTypes.mc_point3_t
                Shape/size of the accumulator along the x, y and z axis.
            multiplier: mc_fp_t
                Multiplier that is applied to the product of the photon
                packet weight and path length traversed through the voxel
                before converting to integer weight. This value should be
                approximately the inverse of the voxel size (the smallest dim).
                This should eliminate degradation of precission due to integer
                accumulators.  
            offset: mc_size_t
                Offset of the first element of the accumulator buffer.
            k: McTypes.mc_int_t
                Integer factor that converts floating point photon packet
                weight to integer value compatible with the sampling volume
                accumulators.
            '''
            _fields_ = [
                ('top_left', T.mc_point3f_t),
                ('voxel_size', T.mc_point3f_t),
                ('shape', T.mc_point3_t),
                ('multiplier', T.mc_fp_t),
                ('offset', T.mc_size_t),
                ('k', T.mc_int_t),
            ]

        return ClSamplingVolume

    def __init__(self, xaxis: Axis or 'SamplingVolume', yaxis: Axis = None,
                 zaxis: Axis = None):
        '''
        Sampling volume object constructor. Default constructor disables the
        functionality by creating a zero-size accumulator array.

        Parameters
        ----------
        xaxis: Axis or SamplingVolume
            Axis that defines accumulators along the x axis.
            If a SamplingVolume instance, a new copy is created.
        yaxis: Axis
            Axis that defines accumulators along the y axis.
        zaxis: Axis
            Axis that defines accumulators along the z axis.

        Note
        ----
            Sampling volume accumulator buffer data type is inherited from the
            Monte Carlo simulator mc_accu_t type.
        '''
        data = None
        weight = 0.0
        if isinstance(xaxis, SamplingVolume):
            sv = xaxis
            xaxis = Axis(sv.xaxis)
            yaxis = Axis(sv.yaxis)
            zaxis = Axis(sv.zaxis)

            if sv.data is not None:
                data = np.copy(sv.data)

            weight = sv.weight

        self._x_axis = xaxis
        self._y_axis = yaxis
        self._z_axis = zaxis

        self._data = data
        self._weight = weight

        self._k = mctypes.McFloat32.mc_fp_maxint

        if self._x_axis.n*self._y_axis.n*self._z_axis.n <= 0:
            raise ValueError('Sampling volume accumulator array has one or '\
                             'dimensions equal to zero!')

    def todict(self) -> dict:
        '''
        Save the sampling volume configuration without the accumulator data to
        a dictionary.

        Returns
        -------
        data: dict
            Sampling volume configuration as a dictionary.
        '''
        return {
            'type':'SamplingVolume',
            'xaxis':self._x_axis.todict(),
            'yaxis':self._y_axis.todict(),
            'zaxis':self._z_axis.todict()
        }
    
    @classmethod
    def fromdict(cls, data) -> 'SamplingVolume':
        '''
        Create a new instance of :py:class:`SamplingVolume` from dict data
        exported by the :py:meth:`todict` method.

        Parameters
        ----------
        data: dict
            Instance data exported to a dict.

        Returns
        -------
        sv: SamplingVolume
            A new :py:class:`SamplingVolume` instance initialized with the data.
        '''
        data_ = dict(data)

        type_name = data_.pop('type')
        if type_name != cls.__name__:
            raise TypeError('Expected data for type "{}" but got "{}"!'.format(
                cls.__name__, type_name))

        return SamplingVolume(
            Axis.fromdict(data_.pop('xaxis')),
            Axis.fromdict(data_.pop('yaxis')),
            Axis.fromdict(data_.pop('zaxis')),
            **data_
        )

    def _get_shape(self) -> Tuple[int, int, int]:
        return (self._z_axis.n, self._y_axis.n, self._x_axis.n,)
    shape = property(_get_shape, None, None, 'Sampling volume array shape.')

    def _get_x(self) -> np.ndarray:
        return self._x_axis.centers
    x = property(_get_x, None, None, 'Accumulator centers along the x axis.')

    def _get_y(self) -> np.ndarray:
        return self._y_axis.centers
    y = property(_get_y, None, None, 'Accumulator centers along the y axis.')

    def _get_z(self) -> np.ndarray:
        return self._z_axis.centers
    z = property(_get_z, None, None, 'Accumulator centers along the z axis.')

    def _get_x_axis(self) -> Axis:
        return self._x_axis
    xaxis = property(_get_x_axis, None, None,
                     'Accumulator axis object along the x axis.')

    def _get_y_axis(self) -> Axis:
        return self._y_axis
    yaxis = property(_get_y_axis, None, None,
                     'Accumulator axis object along the y axis.')

    def _get_z_axis(self) -> Axis:
        return self._z_axis
    zaxis = property(_get_z_axis, None, None,
                     'Accumulator axis object along the z axis.')

    def _get_k(self) -> int:
        return self._k
    k = property(_get_k, None, None,
                 'Samping volume floating point to '
                 'accumulator integer conversion coefficient.')

    def _get_data(self) -> np.ndarray:
        return self._data
    def _set_data(self, data: np.ndarray):
        self._data = data
    data = property(_get_data, _set_data, None,
                    'Raw sampling volume accumulator data if any.')

    def _get_weight(self) -> float:
        return self._weight
    def _set_weight(self, w: float):
        self._weight = w
    weight = property(_get_weight, _set_weight, None,
                    'Total weight of the accumulated photon packets.')

    def clear(self):
        '''
        Clear the data buffer (fill with 0) if one exists.
        '''
        if self._data is not None:
            self._data.fill(0)

    def update_data(self, mc: mcobject.McObject,
                    accumulators: List[np.ndarray], total_weight: int,
                    **kwargs):
        '''
        Update the sampling volume accumulator data with simulation results.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that produced the data.
        accumulators: List[np.ndarray]
            List of allocated accumulators (this implementation uses only one).
        total_weight: int
            Total detected weight of the processed photon packets in integer
            units. 
        kwargs: dict
            Additional keyword arguments not used by this implementation.
        '''
        multiplier = self._multiplier(mc)

        new_data = np.reshape(accumulators[0], self.shape)
        if self._data is not None:
            self._data += new_data*(1.0/(self.k*multiplier))
            self._weight += float(total_weight)/self._k
        else:
            self._data = new_data*(1.0/(self.k*multiplier))
            self._weight = float(total_weight)/self._k

    def _multiplier(self, mc: mcobject.McObject):
        '''
        Multiplier factor for the weight - voxel path length product.
        '''
        l = np.min([self._x_axis.step, self._y_axis.step, self._z_axis.step])
        if l <= 0.0:
            l = 1.0/self._k
        return 1.0/l

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`SamplingVolume.cl_type` for
        a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: ClSamplingVolume
            Structure that is filled with the source data.
        buffer: np.ndarray
            Accumulator buffer or None. Should be checked for proper size. Use
            :py:attr:`xopto.mcml.Mc.types.np_accu` method to determine the
            numpy type of the accumulator used in the Monte Carlo simulator.

        Returns
        -------
        target: ClSamplingVolume
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        allocation = mc.cl_allocate_rw_accumulator_buffer(self, self.shape)
        target.offset = allocation.offset

        target.voxel_size.x = self._x_axis.step
        target.voxel_size.y = self._y_axis.step
        target.voxel_size.z = self._z_axis.step

        target.top_left.x = self._x_axis.start
        target.top_left.y = self._y_axis.start
        target.top_left.z = self._z_axis.start

        target.shape.x = self._x_axis.n
        target.shape.y = self._y_axis.n
        target.shape.z = self._z_axis.n

        target.multiplier = self._multiplier(mc)

        target.k = self._k

        return target

    def plot(self, scale: str = 'log', axis: str ='z', autoscale: bool = True,
             show: bool = True):
        '''
        Show sampling volume slices or integral projections.

        Parameters
        ----------
        scale: str
            Data scaling can be "log" for logarithmic or "lin" for linear.
        axis: str
            The axis of slicing ("x", "y" or "z") or a projection along the
            selected coordinate axis ("xproj", "yproj", "zproj").
            Alternatively, specify the projection plane as one of
            ("xy", "xz", or "yz").
        autoscale: bool
            Scale the color coding of individual slices to the corresponding
            range of weights. If True, the color coding changes from slice
            to slice.
        show: bool 
        '''
        from xopto.util import sliceview

        if axis == 'xy': axis = 'zproj'
        if axis == 'xz': axis = 'yproj'
        if axis == 'yz': axis = 'xproj'

        ax = {'z':0, 'y':1, 'x':2}.get(axis[0], 0)
        title = 'Slice {{slice}}/{} @ {} = {{pos:.6f}}'.format(
            self.data.shape[ax], axis)
        logscale = scale == 'log'

        fig = None

        if ax == 0:
            extent = [self._x_axis.start, self._x_axis.stop,
                      self._y_axis.start, self._y_axis.stop]
            slices = self._z_axis.centers
            xlabel, ylabel = 'x', 'y'
        elif ax == 1:
            extent = [self._x_axis.start, self._x_axis.stop,
                      self._z_axis.start, self._z_axis.stop]
            slices = self._y_axis.centers
            xlabel, ylabel = 'x', 'z'
        elif ax == 2:
            extent = [self._y_axis.start, self._y_axis.stop,
                      self._z_axis.start, self._z_axis.stop]
            slices = self._x_axis.centers
            xlabel, ylabel = 'y', 'z'

        if axis in ('xproj', 'yproj', 'zproj'):
            import matplotlib.pyplot as pp

            fig = pp.figure()
            data_slice = self.data.sum(axis=ax)
            low = data_slice.min()
            if scale == 'log':
                if low < 0:
                    data_slice = np.log(data_slice + (1.0 - low))
                else:
                    data_slice = np.log(data_slice + 1.0)

            pp.imshow(data_slice, extent=extent, origin='lower')
            pp.xlabel(xlabel)
            pp.ylabel(ylabel)
            pp.title('Integral projection along the {:s} axis'.format(axis[0]))
            fig.canvas.manager.set_window_title('Sampling Volume SliceView')
            pp.tight_layout()
            if show:
                pp.show()

        else:
            sv = sliceview.SliceView(
                self.data, axis=ax, slices=slices, title=title, logscale=logscale,
                extent=extent, xlabel=xlabel, ylabel=ylabel, origin='lower',
                autoscale=autoscale)
            sv.fig.canvas.manager.set_window_title('Sampling Volume SliceView')
            if show:
                sv.show()

    def __str__(self):
        return "SamplingVolume(xaxis={}, yaxis={}, zaxis={})".format(
            self._x_axis, self._y_axis, self._z_axis)

    def __repr__(self):
        return self.__str__() + \
            ' # object at 0x{:>08X}.'.format(id(self))
