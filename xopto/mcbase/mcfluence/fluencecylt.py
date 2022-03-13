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

from typing import Dict, List, Tuple

import numpy as np
from xopto.mcbase import cltypes
from xopto.mcbase import mctypes
from xopto.mcbase import mcoptions
from xopto.mcbase import mcobject
from xopto.mcbase.mcutil.axis import Axis


class FluenceCylt(mcobject.McObject):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClFluenceCylt(cltypes.Structure):
            '''
            OpenCL structure that used by the Monte carlo simulator.

            Fields
            ------
            center: McTypes.cl_point2f
                Center of the cylindrical coordinate system in the x-y plane.
            r_min: McTypes.cl_float
                Minimum r coordinate.
            fi_min: McTypes.cl_float
                Minimum polar angle coordinate.
            z_min: McTypes.cl_float
                Minimum z coordinate.
            t_min: mc_fp_t
                Left edge of the first bin along the t axis.
            inv_dr: McTypes.cl_float
                Inverse spacings of the fluence accumulators along the radial axis.
            inv_dfi: McTypes.cl_float
                Inverse spacings of the fluence accumulators along the
                polar angle axis.
            inv_dz: McTypes.cl_float
                Inverse spacings of the fluence accumulators along the
                z axis.
            inv_dt: McTypes.cl_float
                Inverse spacings of the fluence accumulators in the t
                axis.
            n_r: McTypes.cl_int
                Number of accumulators along the r axis.
            n_fi: McTypes.cl_int
                Number of accumulators along the polar angle axis.
            n_z: McTypes.cl_int
                Number of accumulators along the z axis.
            n_t: McTypes.cl_int
                Number of accumulators along the t axis.
            offset: McTypes.mc_int_t
                Offset of the first element of the fluence accumulator buffer.
            k: McTypes.cl_float
                Integer factor that converts floating point photon packet
                weight to integer value compatible with the fluence
                accumulators.
            '''
            _fields_ = [
                ('center', T.mc_point2f_t),
                ('r_min', T.mc_fp_t),
                ('fi_min', T.mc_fp_t),
                ('z_min', T.mc_fp_t),
                ('t_min', T.mc_fp_t),
                ('inv_dr', T.mc_fp_t),
                ('inv_dfi', T.mc_fp_t),
                ('inv_dz', T.mc_fp_t),
                ('inv_dt', T.mc_fp_t),
                ('n_r', T.mc_int_t),
                ('n_fi', T.mc_int_t),
                ('n_z', T.mc_int_t),
                ('n_t', T.mc_int_t),
                ('offset', T.mc_size_t),
                ('k', T.mc_int_t),
            ]

        return ClFluenceCylt

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McFluence{',
            '	mc_point2f_t center;',
            '	mc_fp_t r_min;',
            '	mc_fp_t fi_min;',
            '	mc_fp_t z_min;',
            '	mc_fp_t t_min;',
            '	mc_fp_t inv_dr;',
            '	mc_fp_t inv_dfi;',
            '	mc_fp_t inv_dz;',
            '	mc_fp_t inv_dt;',
            '	mc_int_t n_r;',
            '	mc_int_t n_fi;',
            '	mc_int_t n_z;',
            '	mc_int_t n_t;',
            '	mc_size_t offset;',
            '	mc_int_t k;',
            '};',
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject):
        return '\n'.join((
            'void dbg_print_fluence(__mc_fluence_mem const McFluence *fluence){',
            '	dbg_print("Cylindrical timer-resolved McFluence fluence:");',
            '	dbg_print_float(INDENT "center.x (mm):", fluence->center.x*1e3f);',
            '	dbg_print_float(INDENT "center.y (mm):", fluence->center.y*1e3f);',
            '	dbg_print_float(INDENT "r_min (mm):", fluence->r_min*1e3f);',
            '	dbg_print_float(INDENT "fi_min (rad):", fluence->fi_min);',
            '	dbg_print_float(INDENT "z_min (mm):", fluence->z_min*1e-3f);',
            '	dbg_print_float(INDENT "t_min (ns):", fluence->t_min*1e9f);',
            '	dbg_print_float(INDENT "inv_dr (1/mm):", fluence->inv_dr*1e-3f);',
            '	dbg_print_float(INDENT "inv_dfi (1/rad):", fluence->inv_dfi);',
            '	dbg_print_float(INDENT "inv_dz (1/mm):", fluence->inv_dz*1e-3f);',
            '	dbg_print_float(INDENT "inv_dt (1/ns):", fluence->inv_dz*1e-9f);',
            '	dbg_print_int(INDENT "n_r:", fluence->n_r);',
            '	dbg_print_int(INDENT "n_fi:", fluence->n_fi);',
            '	dbg_print_int(INDENT "n_z:", fluence->n_z);',
            '	dbg_print_int(INDENT "n_t:", fluence->n_t);',
            '	dbg_print_size_t(INDENT "offset:", fluence->offset);',
            '	dbg_print_int(INDENT "k:", fluence->k);',
            '',
            '	#if MC_FLUENCE_MODE_RATE',
            '		dbg_print(INDENT "Mode: fluence");',
            '	#else',
            '		dbg_print(INDENT "Mode: deposition");',
            '	#endif',
            '};',
            '',
            '#if MC_FLUENCE_MODE_RATE',
            'inline void mcsim_fluence_deposit_at(',
            '	McSim *mcsim, mc_point3f_t const *position, ',
            '	mc_fp_t weight, mc_fp_t mua){',
            '#else',
            'inline void mcsim_fluence_deposit_at(',
            '	McSim *mcsim, mc_point3f_t const *position, mc_fp_t weight){',
            '#endif',
            '	mc_fp_t indexf_r, indexf_fi, indexf_z, indexf_t;',
            '	__mc_fluence_mem McFluence const *fluence = mcsim_fluence(mcsim);',
            '',
            '	mc_fp_t dx = mcsim_position_x(mcsim) - fluence->center.x;',
            '	mc_fp_t dy = mcsim_position_y(mcsim) - fluence->center.y;',
            '',
            '	mc_fp_t r = mc_sqrt(dx*dx + dy*dy);',
            '	mc_fp_t fi = mc_atan2(dy, dx) + FP_PI;',
            '',
            '	mc_fp_t dt = mcsim_optical_pathlength(mcsim)*FP_INV_C - fluence->t_min;',
            '',
            '	indexf_r = (r - fluence->r_min)*fluence->inv_dr;',
            '	indexf_z = (mcsim_position_z(mcsim) - fluence->z_min)*fluence->inv_dz;',
            '	indexf_fi = (fi - fluence->fi_min)*fluence->inv_dfi;',
            '	indexf_t = dt*fluence->inv_dt;',
            '',
            '	if (indexf_r >= 0 && indexf_z >= 0 && ',
            '			indexf_fi >= 0 && indexf_t >= 0 &&',
            '			indexf_r < fluence->n_r && ',
            '			indexf_z < fluence->n_z && '
            '			indexf_fi < fluence->n_fi &&',
            '			indexf_t < fluence->n_t){',
            '		mc_size_t index, index_r, index_fi, index_z, index_t;',
            '		index_r = mc_uint(indexf_r);',
            '		index_z = mc_uint(indexf_z);',
            '		index_fi = mc_uint(indexf_fi);',
            '		index_t = mc_uint(indexf_t);',
            '',
            '		index = ((index_z*fluence->n_fi + index_fi)*fluence->n_r + index_r)*fluence->n_t + index_t;',
            '		#if MC_ENABLE_DEBUG',
            '		mc_point4_t index_xyzt = {index_r, index_fi, index_z, index_t};',
            '		dbg_print("Fluence depositing:");',
            '		dbg_print_float(INDENT     "weight                     :", weight);',
            '		dbg_print_point4(INDENT    "voxel address (r, fi, z, t):", &index_xyzt);',
            '		dbg_print_int(INDENT       "flat index                 :", index);',
            '		dbg_print_size_t(INDENT    "offset                     :", fluence->offset);',
            '		#endif',
            '',
            '		#if MC_FLUENCE_MODE_RATE',
            '		weight *= (mua != FP_0) ? mc_fdiv(FP_1, mua) : FP_0;',
            '		#endif'
            '',
            '		uint32_t ui32w = (uint32_t)(weight*fluence->k + FP_0p5);',
            '',
            '		mcsim_fluence_weight_deposit_ll(mcsim, fluence->offset + index, ui32w);',
            '	};',
            '};',
        ))

    '''
    Maximum integer '0x7FFFFF' (8388607) that can be represented by a floating
    point number is used by default to convert photon packet weight
    (floating point) to accumulator data type (unsigned integer).
    '''

    def cl_options(self, mc: mcobject.McObject) -> mcoptions.RawOptions:
        return [('MC_USE_FLUENCE', True),
                ('MC_FLUENCE_MODE_RATE', self.mode == 'fluence'),
                ('MC_TRACK_OPTICAL_PATHLENGTH', True)]

    def __init__(self, raxis: Axis or 'FluenceCylt', fiaxis: Axis = None,
                 zaxis: Axis = None, taxis: Axis = None,
                 center: Tuple[float, float] = (0.0, 0.0),
                 mode: str = 'deposition'):
        '''
        Cylindrical Fluence object constructor.
        Default constructor disables the fluence functionality by creating
        a zero-size fluence accumulator array.

        Parameters
        ----------
        raxis: Axis or FluenceCylt
            Axis that defines accumulators along the radial axis.
            If FluenceCylt instance, a new copy is created.
        fiaxis: Axis
            Axis that defines accumulators along the polar angle axis. The
            range of this axis typically span from 0 to 2*:math:`\pi`.
        zaxis: Axis
            Axis that defines accumulators along the z axis.
        taxis: Axis
            Axis that defines accumulators along the temporal axis.
        center: Tuple(float, float)
            Center of the cylindrical coordinate system in the x-y plane.
        mode: str from ('deposition', 'fluence')
            Mode that is used to accumulate the photon packet weight:

            - fluence - fluence rate ( 1/m :superscript:`2`)
            - deposition - absorbed energy (sum of photon packet weights
              absorbed in the voxel 1/m :superscript:`3`).

        Note
        ----
            The fluence accumulator buffer data type is inherited from the
            Monte Carlo simulator mc_accu_t type.
        '''
        data = None
        nphotons = 0
        k = mctypes.McFloat32.mc_fp_maxint

        if isinstance(raxis, FluenceCylt):
            fluence = raxis
            raxis = Axis(fluence.raxis)
            fiaxis = Axis(fluence.fiaxis)
            zaxis = Axis(fluence.zaxis)
            taxis = Axis(fluence.taxis)
            center = fluence.center
            nphotons = fluence.nphotons
            mode = fluence.mode
            k = fluence.k

            if fluence.raw is not None:
                data = np.copy(fluence.raw)

        if mode not in ('fluence', 'deposition'):
            raise ValueError(
                'The value of mode parameter must be '
                '"fluence" or "deposition" but got {}!'.format(mode))

        if raxis is None:
            raxis = Axis(0.0, 1.0, 1)
        if fiaxis is None:
            fiaxis = Axis(0.0, 2*np.pi, 1)
        if zaxis is None:
            zaxis = Axis(0.0, 1.0, 1)
        if taxis is None:
            taxis = Axis(0.0, 1.0, 1)

        if raxis.logscale:
            raise ValueError('FluenceCylt does not support logarithmic radial axis!')
        if fiaxis.logscale:
            raise ValueError('FluenceCylt does not support logarithmic polar angle axis!')
        if zaxis.logscale:
            raise ValueError('FluenceCylt does not support logarithmic z axis!')
        if taxis.logscale:
            raise ValueError('FluenceCylt does not support logarithmic t axis!')

        self._r_axis = raxis
        self._fi_axis = fiaxis
        self._z_axis = zaxis
        self._t_axis = taxis
        self._center = np.zeros((2,))
        self._set_center(center)

        self._mode = mode

        self._data = data
        self._nphotons = nphotons

        self._k = k

        if self._r_axis.n*self._fi_axis.n*self._z_axis.n*self._t_axis.n <= 0:
            raise ValueError('Fluence accumulator array has one or more array '\
                             'dimensions equal to zero!')

    def _get_nphotons(self) -> int:
        return self._nphotons
    nphotons = property(_get_nphotons, None, None,
                        'The number of photon packets that produced '\
                        'the raw data accumulator content.')

    def _get_mode(self) -> int:
        return self._mode
    mode = property(_get_mode, None, None, 'The accumulator mode.')

    def todict(self) -> dict:
        '''
        Save the fluence configuration without the accumulator data to
        a dictionary.

        Returns
        -------
        data: dict
            Fluence configuration as a dictionary.
        '''
        return {
            'type':'FluenceCylt',
            'mode':self._mode,
            'raxis':self._x_axis.todict(),
            'fiaxis':self._y_axis.todict(),
            'zaxis':self._z_axis.todict(),
            'taxis':self._t_axis.todict(),
            'center':self._center.tolist()
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'FluenceCylt':
        '''
        Create a FluenceCylt instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`FluenceCylt.todict` method.
        '''
        data_ = dict(data)
        fluence_type = data_.pop('type')
        if fluence_type != cls.__name__:
            raise TypeError('Expected "{}" type bot got "{}"!'.format(
                cls.__name__, fluence_type))
        r_axis = Axis.fromdict(data_.pop('raxis'))
        fi_axis = Axis.fromdict(data_.pop('fiaxis'))
        z_axis = Axis.fromdict(data_.pop('zaxis'))
        t_axis = Axis.fromdict(data_.pop('taxis'))

        return cls(r_axis, fi_axis, z_axis, t_axis, **data_)

    def _get_shape(self) -> Tuple[int, int, int, int]:
        return (self._z_axis.n, self._fi_axis.n, self._r_axis.n, self._t_axis.n)
    shape = property(_get_shape, None, None, 'Fluence array shape.')

    def _get_center(self) -> np.ndarray:
        return self._center
    def _set_center(self, center: np.ndarray or Tuple[float, float]):
        self._center[:] = center
    center = property(_get_center, _set_center, None,
                      'Center of the cylindrical coordinate system in the '
                      'x-y plane.')

    def _get_r(self) -> np.ndarray:
        return self._r_axis.centers
    r = property(_get_r, None, None, 'Accumulator centers along the radial axis.')

    def _get_dr(self) -> np.ndarray:
        return abs(self._r_axis.step)
    dr = property(_get_dr, None, None,
                  'The size of voxels along the radial axis.')

    def _get_fi(self) -> np.ndarray:
        return self._fi_axis.centers
    fi = property(_get_fi, None, None,
                  'Accumulator centers along the polar angle axis.')

    def _get_dfi(self) -> np.ndarray:
        return abs(self._fi_axis.step)
    dfi = property(_get_dfi, None, None,
                   'The size of voxels along the polar angle axis.')

    def _get_z(self) -> np.ndarray:
        return self._z_axis.centers
    z = property(_get_z, None, None, 'Accumulator centers along the z axis.')

    def _get_dz(self) -> np.ndarray:
        return abs(self._z_axis.step)
    dz = property(_get_dz, None, None, 'The size of voxels along the z axis.')

    def _get_t(self) -> np.ndarray:
        return self._t_axis.centers
    t = property(_get_t, None, None, 'Accumulator centers along the t axis.')

    def _get_dt(self) -> np.ndarray:
        return abs(self._t_axis.step)
    dt = property(_get_dt, None, None, 'The size of voxels along the t axis.')

    def _get_r_axis(self) -> Axis:
        return self._r_axis
    raxis = property(_get_r_axis, None, None,
                     'Accumulator axis object along the radial axis.')

    def _get_fi_axis(self) -> Axis:
        return self._fi_axis
    fiaxis = property(_get_fi_axis, None, None,
                     'Accumulator axis object along the polar angle axis.')

    def _get_z_axis(self) -> Axis:
        return self._z_axis
    zaxis = property(_get_z_axis, None, None,
                     'Accumulator axis object along the z axis.')

    def _get_t_axis(self) -> Axis:
        return self._t_axis
    taxis = property(_get_t_axis, None, None,
                     'Accumulator axis object along the t axis.')

    def _get_k(self) -> int:
        return self._k
    def _set_k(self, k: int):
        self._k = max(1, min(int(k), int(2**31 - 1)))
    k = property(_get_k, _set_k, None, 'Fluence floating point to accumulator'
                                       'integer conversion coefficient.')

    def _get_raw_data(self) -> np.ndarray:
        return self._data
    def _set_raw_data(self, data: np.ndarray):
        self._data = data
    raw = property(_get_raw_data, _set_raw_data, None,
                   'Raw fluence accumulator data if any.')

    def _get_data(self):
        r = self._r_axis.edges
        k = 1.0/(self.nphotons*(r[1:]**2 - r[:-1]**2)*self.dfi*self.dz*self.dt)
        k.shape = (1, 1, k.size, 1)
        return self._data*k
    data = property(_get_data, None, None,
                    'Fluence accumulator - deposition or fluence rate.')

    def update_data(self, mc: mcobject.McObject,
                    data: Dict[np.dtype, List[np.ndarray]],
                    nphotons: int, **kwargs):
        '''
        Update fluence accumulator data with simulation results.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that produced the data.
        data: Dict[np.dtype, List[np.ndarray]]
            List of allocated accumulators (this implementation uses only
            one accumulator buffer).
        nphotons: int
            The number of photon packets that produced the raw data
            accumulator content.
        kwargs: dict
            Additional keyword arguments not used by this implementation.
        '''
        accumulators = data[np.dtype(mc.types.np_accu)]

        if self._data is not None:
            self._data.flat += accumulators[0]*(1.0/self.k)
            self._nphotons += nphotons
        else:
            self._data = accumulators[0]*(1.0/self.k)
            self._data.shape = self.shape
            self._nphotons = nphotons

    def update(self, obj : 'FluenceCylt'):
        '''
        Update the fluence accumulator with data from the given fluence object.

        Parameters
        ----------
        obj: FluenceCylt
            Update the fluence accumulator of this instance with the data
            from fluence instance obj.
        '''
        if self._data is not None:
            if self.shape != obj.shape:
                raise TypeError(
                    'Cannot update with fluence data of incompatible shape!')

            self._data += obj.raw
            self._nphotons += obj.nphotons
        else:
            self._data = obj.raw
            self._nphotons = obj.nphotons

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator. See the :py:meth:`FluenceCylt.cl_type`
        for a detailed list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: ClFluenceCylt
            CStructure that is filled with the source data.
        buffer: np.ndarray
            Accumulator buffer or None. Should be checked for proper size. Use
            py:attr:`mc.types.np_accu` attribute to determine the
            numpy type of the accumulator used in the Monte Carlo simulator.

        Returns
        -------
        target: ClFluenceCylt
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        allocation = mc.cl_allocate_rw_accumulator_buffer(self, self.shape)
        target.offset = allocation.offset

        target.center.x = self._center[0]
        target.center.y = self._center[1]

        target.r_min = self._r_axis.start
        target.fi_min = self._fi_axis.start
        target.z_min = self._z_axis.start
        target.t_min = self._t_axis.start

        target.inv_dr = 1.0/self._r_axis.step
        target.inv_dfi = 1.0/self._fi_axis.step
        target.inv_dz = 1.0/self._z_axis.step
        target.inv_dt = 1.0/self._t_axis.step

        target.n_r = self._r_axis.n
        target.n_fi = self._fi_axis.n
        target.n_z = self._z_axis.n
        target.n_t = self._t_axis.n

        target.k = self._k

        return target

    def plot(self, scale: str = 'log', axis: str ='z', autoscale: bool = True,
             show: bool = True):
        '''
        Show fluence slices or integral projections.

        Parameters
        ----------
        scale: str
            Data scaling can be "log" for logarithmic or "lin" for linear.
        axis: str
            The axis of slicing ("z", "fi" or "r") or a projection along the
            selected coordinate axis ("rproj", "fiproj", "zproj").
            Alternatively, specify the projection plane as one of
            ("rfi", "rz", or "fiz").
        autoscale: bool
            Scale the color coding of individual slices to the corresponding
            range of weights. If True, the color coding changes from slice
            to slice.
        show: bool 
        '''
        from xopto.util import sliceview

        data = self.data

        if axis == 'rfi': axis = 'zproj'
        if axis == 'rz': axis = 'fiproj'
        if axis == 'fiz': axis = 'rproj'

        ax = {'z':0, 'f':1, 'r':2}.get(axis[0], 0)

        logscale = scale == 'log'

        fig = None

        if ax == 0:
            extent = [self._r_axis.start, self._r_axis.stop,
                      self._fi_axis.start, self._fi_axis.stop]
            slices = self._z_axis.centers
            xlabel, ylabel = 'r', 'fi'
        elif ax == 1:
            extent = [self._r_axis.start, self._r_axis.stop,
                      self._z_axis.start, self._z_axis.stop]
            slices = self._fi_axis.centers
            xlabel, ylabel = 'r', 'z'
        elif ax == 2:
            extent = [self._fi_axis.start, self._fi_axis.stop,
                      self._z_axis.start, self._z_axis.stop]
            slices = self._r_axis.centers
            xlabel, ylabel = 'fi', 'z'

        window_title = 'FluenceCylt SliceView - {} mode'.format(self.mode)

        if axis in ('rproj', 'fiproj', 'zproj'):
            data_slice = data.sum(axis=ax)
            polar = axis == 'zproj'
            if polar:
                # axis order must be (r, fi, t) but is currently (fi, r, t)
                data_slice = np.transpose(data_slice, (1, 0, 2))

            title = 'Slice {{slice}}/{} @ {} = {{pos:.1f}} ps'.format(
                data_slice.shape[2], axis)

            r, fi = self.r, self.fi
            R, Fi = np.meshgrid(r, fi, indexing='ij')
            sv = sliceview.SliceViewCyl(
                data_slice,
                polar=polar,
                axis=2, slices=self.t*1e12,
                title=title, logscale=logscale,
                extent=extent, xlabel=xlabel, ylabel=ylabel, origin='lower',
                autoscale=autoscale, R=R, Fi=Fi, aspect='auto')
            sv.fig.canvas.manager.set_window_title(window_title)
            if show:
                sv.show()

        else:
            title = 'Slice {{slice2}}/{} @ {} = {{pos2:.6f}}'.format(
                data.shape[ax], axis)

            polar = ax == 0
            if polar:
                # axis order must be (z, r, fi, t) but is currently (z, fi, r, t)
                data = np.transpose(data, (0, 2, 1, 3))

            r, fi = self.r, self.fi
            R, Fi = np.meshgrid(r, fi, indexing='ij')
            sv = sliceview.DualSliceViewCyl(
                data, polar=ax == 0,
                axis1=3, slices1=self.t*1e12,
                slice1_label='Time', slice1_valfmt='%.1f ps',
                axis2=ax, slices2=slices,
                title=title, logscale=logscale,
                extent=extent, xlabel=xlabel, ylabel=ylabel, origin='lower',
                autoscale=autoscale, R=R, Fi=Fi, aspect='auto')
            sv.fig.canvas.manager.set_window_title(window_title)
            if show:
                sv.show()

    def __str__(self):
        return "FluenceCylt(raxis={}, fiaxis={}, zaxis={}, taxis={}, " \
               "center={})".format(
                   self._r_axis, self._fi_axis, self._z_axis, self._t_axis,
                   self._center)

    def __repr__(self):
        return self.__str__() + \
            ' # object at 0x{:>08X}.'.format(id(self))
