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

from typing import Tuple, List

import numpy as np

from xopto.mcml.mcdetector.base import Detector
from xopto.mcml import cltypes, mctypes, mcobject, mcoptions
from xopto.mcml.mcutil.fiber import MultimodeFiber, FiberLayout
from xopto.mcml.mcutil import geometry
from xopto.mcml.mcutil import axis


class FiberArrayPl(Detector):
    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        T = mc.types
        n = self.n
        class ClFiberArrayPl(cltypes.Structure):
            '''
            Structure that that represents a detector in the Monte Carlo
            simulator core.

            Fields
            ------
            transformation: mc_point3f_t*n
                Transforms coordinates from Monte Carlo to the detector / fiber
                coordinates for each fiber.
            core_position: mc_point2f_t*n
                Position of the individual fibers.
            core_r_squared: mc_fp_t*n
                Squared radius of the optical fibers.
            cos_min: mc_fp_t*n
                Cosine of the maximum acceptance angle (in the detector / fiber
                coordinate space).
            pl_min: mc_fp_t
                The leftmost edge of the first optical path length accumulator.
            inv_dpl: mc_fp_t
                Inverse of the width of the optical path length accumulators.
            n_pl: mc_size_t
                The number of path length accumulators.
            offset: mc_int_t
                The offset of the first accumulator in the Monte Carlo
                detector buffer.
            pl_log_scale: mc_int_t
                A flag indicating logarithmic scale of the path length axis.
            '''
            _fields_ = [
                ('transformation', T.mc_matrix3f_t*n),
                ('core_position', T.mc_point2f_t*n),
                ('core_r_squared', T.mc_fp_t*n),
                ('cos_min', T.mc_fp_t*n),
                ('pl_min', T.mc_fp_t),
                ('inv_dpl', T.mc_fp_t),
                ('n_pl', T.mc_size_t),
                ('offset', T.mc_size_t),
                ('pl_log_scale', T.mc_int_t),
            ]
        return ClFiberArrayPl

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        n = self.n
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES Mc{}Detector{{'.format(Loc),
            '	mc_matrix3f_t transformation[{}];'.format(n),
            '	mc_point2f_t core_position[{}];'.format(n),
            '	mc_fp_t core_r_squared[{}];'.format(n),
            '	mc_fp_t cos_min[{}];'.format(n),
            '	mc_fp_t pl_min;',
            '	mc_fp_t inv_dpl;',
            '	mc_size_t n_pl;',
            '	mc_size_t offset;',
            '	mc_int_t pl_log_scale;',
            '};'
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of the detector in the Monte Carlo simulator.
        '''
        loc = self.location
        Loc = loc.capitalize()
        n = self.n
        return '\n'.join((
            'void dbg_print_{}_detector(__mc_detector_mem const Mc{}Detector *detector){{'.format(loc, Loc),
            '	dbg_print("Mc{}Detector - FiberArrayPl fiber array detector:");'.format(Loc),
            '	for(mc_size_t index=0; index < {}; ++index){{'.format(n),
            '		dbg_print_size_t(INDENT "Fiber index:", index);',
            '		dbg_print_matrix3f(INDENT INDENT "transformation:", &detector->transformation[index]);',
            '		dbg_print_point2f(INDENT INDENT "core_position:", &detector->core_position[index]);',
            '		dbg_print_float(INDENT INDENT "core_r_squared (mm2):", detector->core_r_squared[index]*1e6f);',
            '		dbg_print_float(INDENT INDENT "cos_min:", detector->cos_min[index]);',
            '	};',
            '	dbg_print_float(INDENT "pl_min (um)", detector->pl_min*1e6f);',
            '	dbg_print_float(INDENT "inv_dpl (1/um)", detector->inv_dpl*1e-6f);',
            '	dbg_print_size_t(INDENT "n_pl:", detector->n_pl);',
            '	dbg_print_size_t(INDENT "offset:", detector->offset);',
            '	dbg_print_int(INDENT "pl_log_scale:", detector->pl_log_scale);',
            '};',
            '',
            'inline void mcsim_{}_detector_deposit('.format(loc),
            '		McSim *mcsim, ',
            '		mc_point3f_t const *pos, mc_point3f_t const *dir, ',
            '		mc_fp_t weight){',
            '',
            '	__global mc_accu_t *address;',
            '',
            '	dbg_print_status(mcsim, "{} FiberArrayPl fiber array detector hit");'.format(loc),
            '',
            '	__mc_detector_mem const Mc{}Detector *detector = '.format(Loc),
            '		mcsim_{}_detector(mcsim);'.format(loc),
            '',
            '',
            '	mc_size_t fiber_index = {}; /* invalid index ... no fiber hit */'.format(n),
            '',
            '	mc_fp_t dx, dy, r_squared;',
            '	mc_point3f_t mc_pos, detector_pos;',
            '',
            '	pragma_unroll_hint({})'.format(n),
            '	for(mc_size_t index=0; index < {}; ++index){{'.format(n),
            '		mc_pos.x = pos->x - detector->core_position[index].x;',
            '		mc_pos.y = pos->y - detector->core_position[index].y;',
            '		mc_pos.z = FP_0;',
            '',
            '		transform_point3f(',
            '			&detector->transformation[index], &mc_pos, &detector_pos);',
            '		dx = detector_pos.x;',
            '		dy = detector_pos.y;',
            '		r_squared = dx*dx + dy*dy;',
            '',
            '		if (r_squared <= detector->core_r_squared[index]){',
            '			/* hit this fiber */',
            '			fiber_index = index;',
            '			break;',
            '		};',
            '	};',
            '',
            '	if (fiber_index >= {})'.format(n),
            '		return;',
            '',
            '	mc_fp_t pl = mcsim_optical_pathlength(mcsim);',
            '	if (detector->pl_log_scale)',
            '	    pl = mc_log(mc_fmax(pl, FP_PLMIN));',
            '	mc_int_t pl_index = mc_int((pl - detector->pl_min)*detector->inv_dpl);',
            '	pl_index = mc_clip(pl_index, 0, detector->n_pl - 1);',
            '',
            '	mc_size_t index = pl_index*{} + fiber_index;'.format(n),
            '',
            '	address = mcsim_accumulator_buffer_ex(',
            '		mcsim, detector->offset + index);',
            '',
            '	/* Transfor direction vector component z into the detector space. */',
            '	mc_fp_t pz = transform_point3f_z(',
            '		&detector->transformation[fiber_index], dir);',
            '	dbg_print_float("Packet direction z:", pz);',
            '	uint32_t ui32w = weight_to_int(weight)*',
            '		(detector->cos_min[fiber_index] <= mc_fabs(pz));',
            '',
            '	if (ui32w > 0){',
            '		dbg_print("{} FiberArrayPl fiber array detector depositing:");'.format(loc),
            '		dbg_print_uint(INDENT "uint weight:", ui32w);',
            '		dbg_print_size_t(INDENT "to fiber index:", fiber_index);',
            '		accumulator_deposit(address, ui32w);',
            '	};',
            '};',
        ))

    def cl_options(self, mc, target=None) -> mcoptions.RawOptions:
        '''
        OpenCL kernel options defined by this object.
        '''
        return [('MC_TRACK_OPTICAL_PATHLENGTH', True)]

    def __init__(self, fibers: List[FiberLayout], plaxis: axis.Axis = None):
        '''
        Optical fiber probe detector with an array of optical fibers that
        are optionally tilted(direction parameter). The optical fibers are
        always polished in a way that forms a tight optical contact with
        the surface of the sample.

        Parameters
        ----------
        fibers: List[FiberLayout]
            List of optical fiber configuration that form the fiber array.
        plaxis: axis.Axis
            Object that defines the accumulators along the optical path length
            axis (this axis supports log-scale).
        '''
        if isinstance(fibers, FiberArrayPl):
            la = fibers
            fibers = la.fibers
            raw_data = np.copy(la.raw)
            nphotons = la.nphotons
            plaxis = la.plaxis
        else:
            if plaxis is None:
                plaxis = axis.Axis(0.0, 1.0, 1)

            nphotons = 0
            raw_data = np.zeros((plaxis.n, len(fibers)))

        super().__init__(raw_data, nphotons)

        self._fibers = fibers
        self._pl_axis = plaxis

    def update(self, other: 'FiberArrayPl' or dict):
        '''
        Update this detector configuration from the other detector. The
        other detector must be of the same type as this detector or a dict with
        appropriate fields.

        Parameters
        ----------
        other: FiberArrayPl or dict
            This source is updated with the configuration of the other source.
        '''
        if isinstance(other, FiberArrayPl):
            self.fibers = other.fibers
        elif isinstance(other, dict):
            self.fibers = other.get('fibers', self.fibers)

    def _get_fibers(self) -> List[FiberLayout]:
        return self._fibers
    def _set_fibers(self, fibers: List[FiberLayout]):
        if len(self._fibers) != len(fibers):
            raise ValueError('The number of optical fibers must not change!')
        if type(self._fibers[0]) != type(fibers[0]):
            raise TypeError('The type of optical fibers must not change!')
        self._fibers[:] = fibers
    fibers = property(_get_fibers, _set_fibers, None,
                     'List of optical fibers.')

    def _get_n_fiber(self) -> int:
        return len(self._fibers)
    n = property(_get_n_fiber, None, None,
                 'Number of optical fiber in the array.')

    def __len__(self):
        return len(self._fibers)

    def __iter__(self):
        return iter(self._fibers)

    def _get_plaxis(self) -> axis.Axis:
        return self._pl_axis
    plaxis = property(_get_plaxis, None, None, 'Path length axis object.')

    def _get_pl(self):
        return self._pl_axis.centers
    pl = property(_get_pl, None, None,
                  'Centers of the optical pathlength axis accumulators.')

    def _get_pledges(self):
        return self._pl_axis.edges
    pledges = property(_get_pledges, None, None,
                       'Edges of the optical pathlength axis accumulators.')

    def _get_npl(self):
        return self._pl_axis.n
    npl = property(_get_npl, None, None,
                   'Number of accumulators in the optical pathlength axis.')

    def check(self):
        '''
        Check if the configuration has errors and raise exceptions if so.
        '''
        for fiber_cfg in self._fibers:
            for other_cfg in self._fibers:
                if other_cfg != fiber_cfg:
                    d = np.linalg.norm(fiber_cfg.position - other_cfg.position)
                    if d < max(fiber_cfg.fiber.dcladding, other_cfg.fiber.dcladding):
                        raise ValueError('Some of the fibers in the '
                                         'detector array overlap!')

        return True

    def fiber_position(self, index: int) -> Tuple[float, float]:
        '''
        Returns the position of the fiber center as a tuple (x, y).

        Parameters
        ----------
        index: int
            Fiber index from 0 to n -1.

        Returns
        -------
        position: (float, float)
            The position of the fiber center as a tuple (x, y).
        '''
        n = len(self._fibers)
        if index >= n or index < -n:
            raise IndexError('The fiber index is out of valid range!')
        return tuple(self._fibers[index].position[:2])

    def __getitem__(self, what):
        return self._fibers[what]

    def __setitem__(self, what, value: FiberLayout or List[FiberLayout]):
        self._fibers[what] = value

    def _get_normalized(self) -> np.ndarray:
        return self.raw*(1.0/self.nphotons)
    normalized = property(_get_normalized, None, None, 'Normalized.')
    reflectance = property(_get_normalized, None, None, 'Reflectance.')
    transmittance = property(_get_normalized, None, None, 'Transmittance.')

    def cl_pack(self, mc: mcobject.McObject,
                target : cltypes.Structure = None) -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`FiberArrayPl.cl_type` method for a detailed list
        of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Ctypes structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        for index, fiber_cfg in enumerate(self._fibers):
            adir = fiber_cfg.direction[0], fiber_cfg.direction[1], \
                   abs(fiber_cfg.direction[2])
            T = geometry.transform_base(adir, (0.0, 0.0, 1.0))

            cos_min = (1.0 - (fiber_cfg.fiber.na/fiber_cfg.fiber.ncore)**2)**0.5
            core_r_squared = 0.25*fiber_cfg.fiber.dcore**2

            target.transformation[index].fromarray(T)
            target.core_position[index].fromarray(fiber_cfg.position)
            target.core_r_squared[index] = core_r_squared
            target.cos_min[index] = cos_min

        allocation = mc.cl_allocate_rw_accumulator_buffer(self, self.shape)
        target.offset = allocation.offset

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
        a dictionary. Use the :meth:`FiberArrayPl.fromdict` method to create a new
        accumulator instance from the returned data.

        Returns
        -------
        data: dict
            Accumulator configuration as a dictionary.
        '''
        return {
            'type':'FiberArrayPl',
            'fibers': [item.todict() for item in self._fibers],
            'pl_axis': self.plaxis.todict(),
        }

    @staticmethod
    def fromdict(data: dict) -> 'FiberArrayPl':
        '''
        Create an accumulator instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`FiberArrayPl.todict` method.
        '''
        detector_type = data.pop('type')
        if detector_type != 'FiberArrayPl':
            raise TypeError(
                'Expected a "FiberArrayPl" type bot got "{}"!'.format(
                    detector_type))
        fibers = data.pop('fibers')
        fibers = [FiberLayout.fromdict(item) for item in fibers]

        pl_axis_data = data.pop('pl_axis')
        pl_axis_type = pl_axis_data.pop('type')
        plaxis = getattr(axis, pl_axis_type)(**pl_axis_data)

        return FiberArrayPl(fibers, plaxis=plaxis, **data)

    def __str__(self):
        return 'FiberArrayPl(fibers={}, plaxis={})'.format(
            self._fibers, self._pl_axis)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
