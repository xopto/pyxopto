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

from xopto.mccyl.mcdetector.base import Detector
from xopto.mccyl import cltypes, mcobject
from xopto.mccyl.mcutil import axis


class FiZ(Detector):
    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        class ClCartesian(cltypes.Structure):
            '''
            Structure that that represents a polar φ-z detector
            in the Monte Carlo simulator core.

            Fields
            ------
            fi_min: mc_fp_t
                The leftmost edge of the first accumulator along the φ axis.
            inv_dfi: mc_fp_t
                Inverse of the spacing between the accumulators along the φ axis.
            z_min: mc_fp_t
                The leftmost edge of the first accumulator along the
                z axis.
            inv_zy: mc_fp_t
                Inverse of the spacing between the accumulators along the
                z axis.
            cos_min: mc_fp_t
                Cosine of the maximum acceptance angle relative to the
                radial direction.
            n_fi: mc_size_t
                The number of accumulators along the φ axis.
            n_z: mc_size_t
                The number of accumulators along the z axis.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            '''
            _pack_ = 1
            _fields_ = [
                ('fi_min', T.mc_fp_t),
                ('inv_dfi', T.mc_fp_t),
                ('z_min', T.mc_fp_t),
                ('inv_dz', T.mc_fp_t),
                ('cos_min', T.mc_fp_t),
                ('n_fi', T.mc_size_t),
                ('n_z', T.mc_size_t),
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
            '	mc_fp_t fi_min;',
            '	mc_fp_t inv_dfi;',
            '	mc_fp_t z_min;',
            '	mc_fp_t inv_dz;',
            '	mc_fp_t cos_min;',
            '	mc_size_t n_fi;',
            '	mc_size_t n_z;',
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
            '	dbg_print("Mc{}Detector - FiZ detector:");'.format(Loc),
            '	dbg_print_float(INDENT "fi_min (deg):", detector->fi_min*FP_RAD2DEG);',
            '	dbg_print_float(INDENT "inv_dfi (1/deg):", detector->inv_dfi*FP_DEG2RAD);',
            '	dbg_print_float(INDENT "z_min (mm):", detector->z_min*1e3f);',
            '	dbg_print_float(INDENT "inv_dz (1/mm):", detector->inv_dz*1e-3f);',
            '	dbg_print_float(INDENT "cos_min:", detector->cos_min);',
            '	dbg_print_float(INDENT "n_fi:", detector->n_fi);',
            '	dbg_print_float(INDENT "n_z:", detector->n_z);',
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
            '	mc_int_t index_fi, index_z;',
            '	mc_size_t index;'
            '',
            '	dbg_print_status(mcsim, "{} FiZ detector hit");'.format(Loc),
            '',
            '	mc_fp_t fi = mc_atan2(pos->y, pos->x);',
            '	index_fi = mc_int((fi - detector->fi_min)*detector->inv_dfi);',
            '	index_fi = mc_clip(index_fi, 0, detector->n_fi - 1);',
            '',
            '	index_z = mc_int((pos->z - detector->z_min)*detector->inv_dz);',
            '	index_z = mc_clip(index_z, 0, detector->n_z - 1);',
            '',
            '	index = index_z*detector->n_fi + index_fi;',
            '',
            '	address = mcsim_accumulator_buffer_ex(',
            '		mcsim, index + detector->offset);',
            '',
            '	mc_point3f_t normal;',
            '	radial_normal(pos, &normal);',
            '',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min <= mc_fabs(mc_dot_point3f(dir, &normal)));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print_uint("{} FiZ detector depositing int:", ui32w);'.format(Loc),
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};'
        ))

    def __init__(self, fiaxis: axis.Axis or 'FiZ', zaxis: axis.Axis = None,
                 cosmin: float = 0.0):
        '''
        2D Cylindrical reflectance/transmittance accumulator in the φ-z plane.

        The grid of the Cartesian accumulators corresponds to a 2D numpy array
        with the first dimension representing the φ axis and second dimension
        representing the z axis (reflectance[z, φ] or transmittance[z, φ]).

        Parameters
        ----------
        fiaxis: axis.Axis
            Object that defines accumulators along the φ axis.
        zaxis: axis.Axis
            Object that defines accumulators along the z axis.
        cosmin: float
            Cosine of the maximum acceptance angle (relative to the direction)
            of the detector computed relative to the detector normal.
        '''
        if isinstance(fiaxis, FiZ):
            detector = fiaxis
            fiaxis = type(detector.fiaxis)(detector.fiaxis)
            zaxis = type(detector.zaxis)(detector.zaxis)
            cosmin = detector.cosmin
            raw_data = np.copy(detector.raw)
            nphotons = detector.nphotons
        else:
            if zaxis is None:
                zaxis = axis.Axis(-1.0, 1.0, 1)

            raw_data = np.zeros((zaxis.n, fiaxis.n))
            nphotons = 0

        super().__init__(raw_data, nphotons)
        self._cosmin = 0.0

        self._fi_axis = fiaxis
        self._z_axis = zaxis

        self._set_cosmin(cosmin)

        self._r_sample = 1.0
        self._accumulators_area = self._fi_axis.step*self._z_axis.step

    def _get_fi_axis(self) -> axis.Axis:
        return self._fi_axis
    fiaxis = property(_get_fi_axis, None, None,
                      'Axis object of the φ axis.')

    def _get_z_axis(self) -> axis.Axis:
        return self._z_axis
    zaxis = property(_get_z_axis, None, None, 'Axis object of the z axis.')

    def _get_cosmin(self) -> Tuple[float, float]:
        return self._cosmin
    def _set_cosmin(self, value: float or Tuple[float, float]):
        self._cosmin = min(max(float(value), 0.0), 1.0)
    cosmin = property(_get_cosmin, _set_cosmin, None,
                      'Cosine of the maximum acceptance angle.')

    def _get_fi(self) -> np.ndarray:
        return self._fi_axis.centers
    fi = property(_get_fi, None, None,
                 'Centers of the accumulators along the φ axis.')

    def _get_z(self) -> np.ndarray:
        return self._z_axis.centers
    z = property(_get_z, None, None,
                 'Centers of the accumulators along the z axis.')

    def _get_fi_edges(self) -> np.ndarray:
        return self._fi_axis.edges
    fiedges = property(_get_fi_edges, None, None,
                       'Edges of the accumulators along the φ axis.')

    def _get_z_edges(self) -> np.ndarray:
        return self._z_axis.edges
    zedges = property(_get_z_edges, None, None,
                      'Edges of the accumulators along the z axis.')

    def _get_nfi(self) -> int:
        return self._fi_axis.n
    nfi = property(_get_nfi, None, None,
                   'Number of accumulators along the φ axis.')

    def _get_nz(self) -> int:
        return self._z_axis.n
    nz = property(_get_nz, None, None,
                  'Number of accumulators along the z axis.')

    def meshgrid(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Returns 2D arrays of z and φ coordinates of the centers of accumulators
        that match the size of the reflectance / transmittance arrays.
        The grid of the Cartesian accumulators corresponds to a 2D numpy array
        with the first dimension representing the z axis and second dimension
        representing the φ axis (reflectance[z, φ] or transmittance[z, φ]).

        Returns
        -------
        z: np.ndarray
            A 2D array of φ coordinates.
        φ: np.ndarray
            A 2D array of z coordinates.
        '''
        Z, Fi = np.meshgrid(self.z, self.fi, indexing='ij')
        return Z, Fi

    def update_data(self, mc: mcobject.McObject, *args, **kwargs):
        # save the sample radius
        self._r_sample = mc.layers[1].d*0.5
        return super().update_data(mc, *args, **kwargs)

    def _get_normalized(self) -> np.ndarray:
        accumulator_area = self._accumulators_area*self._r_sample
        return self.raw*(1.0/(max(self.nphotons, 1.0)*accumulator_area))
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

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

        target.fi_min = self._fi_axis.start
        if self._fi_axis.n > 1:
            target.inv_dfi = 1.0/self._fi_axis.step
        else:
            target.inv_dfi = 0.0
        target.z_min = self._z_axis.start
        if self._z_axis.n > 1:
            target.inv_dz = 1.0/self._z_axis.step
        else:
            target.inv_dz = 0.0
        target.cos_min = self.cosmin
        target.n_fi = self._fi_axis.n
        target.n_z = self._z_axis.n

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
            'type':'FiZ',
            'fi_axis': self._fi_axis.todict(),
            'z_axis': self._z_axis.todict(),
            'cosmin': self._cosmin
        }

    @staticmethod
    def fromdict(data: dict) -> 'FiZ':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the py:meth:`FiZ.todict` method.
        '''
        data_ = dict(data)
        detector_type = data_.pop('type')
        if detector_type != 'FiZ':
            raise TypeError('Expected "FiZ" type but got "{}"!'.format(
                detector_type))
        fiaxis_data = data_.pop('fi_axis')
        fiaxis_type = fiaxis_data.pop('type')

        zaxis_data = data_.pop('z_axis')
        zaxis_type = zaxis_data.pop('type')

        return FiZ(
            getattr(axis, fiaxis_type)(**fiaxis_data),
            getattr(axis, zaxis_type)(**zaxis_data),
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

        extent = [self._fi_axis.start, self._fi_axis.stop,
                  self._z_axis.start, self._z_axis.stop]

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
        ax.set_xlabel('φ')
        ax.set_ylabel('z')
        pp.colorbar(img)

        fig.canvas.manager.set_window_title(
            'FiZ detector - {} - {}'.format(scale, which))

        if show:
            pp.show()

    def __str__(self):
        return 'FiZ(fiaxis={}, zaxis={}, cosmin={})'.format(
                   self._fi_axis, self._z_axis, self._cosmin)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
