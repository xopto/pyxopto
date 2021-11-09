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

from typing import Tuple, List, Dict

import numpy as np

from xopto.mcbase import cltypes
from xopto.mcbase import mcoptions
from xopto.mcbase import mcobject
from xopto.mcbase.mcutil.buffer import BufferAllocation

def _range(value: Tuple[float, float]) -> Tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise TypeError(
            'Range must be specified as a tuple of two float: (low, high)!')
    return float(value[0]), float(value[1])

def _origin2(value: Tuple[float, float]) -> Tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise TypeError(
            'Origin must be specified as a tuple of two float: (x, y)!')
    return float(value[0]), float(value[1])

def _dir3(value: Tuple[float, float, float], k: float = 1.0) \
        -> Tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise TypeError(
            'Direction must be specified as a tuple of three float: (x, y, z)!')
    k = float(k)
    return float(value[0])*k, float(value[1])*k, float(value[2])*k


class Filter(object):
    def __init__(self, x: tuple = None, y: tuple = None,
                 z: tuple = None, pz: tuple = None, r: tuple = None,
                 dir: tuple = None, pl: tuple = None):
        '''
        Trace filter constructor. Creates a filter of traced photon packets
        according to the final position of the photon packet. Set the unused
        filter parameters to None - default value. Filters are
        connected through a logical and.

        Parameters
        ----------
        x: (float, float) or ((float, float), (float, float), ...)
            Final x coordinate of the photon packet must be within the
            specified interval including the boundaries (xmin, xmax). Multiple
            intervals can be specified as a list or tuple of intervals.

        y: (float, float) or ((float, float), (float, float), ...)
            Final y coordinate of the photon packet must be within the
            specified interval including the boundaries (ymin, ymax). Multiple
            intervals can be specified as a list or tuple of intervals.

        z: (float, float) or ((float, float), (float, float), ...)
            Final z coordinate of the photon packet must be within the
            specified interval including the boundaries (zmin, zmax). Multiple
            intervals can be specified as a list or tuple of intervals.

        pz: (float, float)
            Final z component of the photon packet propagation direction must
            be within the specified interval. Note that the z component of
            the propagation direction of photon packets terminated at the
            top sample boundary is negative!. Multiple intervals can be
            specified as a list or tuple of intervals.

        dir: (float, float, (float, float, float))
            The cosine of the angle between the given direction and the final
            direction of the photon packet must be within the specified
            interval (cosmin, cosmax, (px, py, pz)).
            Multiple filters can be specified as a list or tuple.

        r: (float, float, (float, float)) or ((float, float, (float, float)), ...)
            Final r=sqrt((x-xo)**2 + (y - yo)**2) coordinate of the photon
            packet must be within the specified radius (rmin, rmax). The
            radius is calculated with respect to the given origin
            (x0, y0). Multiple filters can be specified as a list or tuple
            of (rmin, rmax, (x0, y0)).

        pl: (float, float) or ((float, float), (float, float), ...)
            Final optical path length traveled by the photon packet must be
            within the specified interval including the boundaries
            (plmin, plmax). Multiple intervals can be specified as a list or
            tuple of intervals.
        '''
        if isinstance(x, Filter):
            filt = x
            x = filt.x
            y = filt.y
            z = filt.z
            pz = filt.pz
            r = filt.r
            dir = filt.dir
            pl = filt.pl

        if x is not None:
            x = self._chek_range('x', x)
        if y is not None:
            y = self._chek_range('y', y)
        if z is not None:
            z = self._chek_range('z', z)
        if pz is not None:
            pz = self._chek_range('pz', pz)
        if r is not None:
            r = self._chek_r(r)
        if dir is not None:
            dir = self._chek_dir(dir)
        if pl is not None:
            pl = self._chek_range('pl', pl)

        self._x = x
        self._y = y
        self._z = z
        self._pz = pz
        self._r = r
        self._dir = dir
        self._pl = pl

    def todict(self) -> dict:
        '''
        Save the filter configuration to a dictionary.

        Returns
        -------
        data: dict
            Filter configuration as a dictionary.
        '''
        return {
            'type':self.__class__.__name__,
            'x':self._x,
            'y':self._y,
            'z':self._z,
            'pz':self._pz,
            'dir':self._dir,
            'r':self._r,
            'pl':self._pl
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Filter':
        '''
        Create a Filter instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`Filter.todict` method.
        '''
        filter_type = data.pop('type')
        if filter_type != 'Filter':
            raise TypeError('Expected "{}" type bot got "{}"!'.format(
                cls.__name__, filter_type))

        return cls(**data)


    def _get_x(self) -> tuple:
        return self._x
    x = property(_get_x, None, None,
                 'Filter of the x coordinate as ((x_min, x_max), ...).')

    def _get_y(self) -> tuple:
        return self._y
    y = property(_get_y, None, None,
                 'Filter of the y coordinate as ((y_min, y_max), ...).')

    def _get_z(self) -> tuple:
        return self._z
    z = property(_get_z, None, None,
                 'Filter of the z coordinate as ((z_min, z_max), ...).')

    def _get_pz(self) -> tuple:
        return self._pz
    pz = property(_get_pz, None, None,
                 'Filter of the z component of the propagation direction as '
                 '(pz_min, pz_max).')

    def _get_dir(self) -> tuple:
        return self._dir
    dir = property(_get_dir, None, None,
                   'Filter of the propagation direction as '
                   '(cosmin, cosmax, (px, py, pz)).')

    def _get_r(self) -> tuple:
        return self._r
    r = property(_get_r, None, None,
                 'Filter of r=sqrt(x*x + y*y) as '
                 '((r_min, r_max, (x_origin, y_origin)), ...).')

    def _get_pl(self) -> tuple:
        return self._pl
    pl = property(_get_pl, None, None,
                  'Filter of the optical pathlength as '
                  '((pl_min, pl_max), ...).')

    def __call__(self, trace_obj: 'Trace', update: bool = True) \
            -> Tuple['Trace', int]:
        '''
        Apply the filter to a trace object and create a new trace object for
        the filtered data only if the value of update parameter is False.
        Parameters
        ----------
        trace_obj: Trace
            The trace object to filter.
        update: bool
            If True (default), the existing trace object is updated. If False,
            a new trace object is created for the filtered data.

        Returns
        -------
        filtered: Trace
            Filtered trace.
        n_dropped: int
            Number of photon packets that matched the filter but exceeded
            the maximum length of the trace.
        '''
        if not isinstance(trace_obj, Trace):
            raise TypeError(
                'Trace filter can be applied only to Trace objects!')

        # photon packets that exceed the trace length
        too_long_mask = trace_obj.overflow
        # start with all the photon packets
        valid_mask = np.ones_like(too_long_mask)

        var_mask = np.zeros_like(valid_mask)

        if self._x is not None:
            for x_range in self._x:
                np.logical_or(
                    self._x_filter(
                        x_range, trace_obj, valid_mask),
                    var_mask, out=var_mask
                )
            np.logical_and(var_mask, valid_mask, out=valid_mask)

        if self._y is not None:
            var_mask.fill(0)
            for y_range in self._y:
                np.logical_or(
                    self._y_filter(
                        y_range, trace_obj, valid_mask),
                    var_mask, out=var_mask
                )
            np.logical_and(var_mask, valid_mask, out=valid_mask)

        if self._z is not None:
            var_mask.fill(0)
            for z_range in self._z:
                np.logical_or(
                    self._z_filter(
                        z_range, trace_obj, valid_mask),
                    var_mask, out=var_mask
                )
            np.logical_and(var_mask, valid_mask, out=valid_mask)

        if self._pz is not None:
            var_mask.fill(0)
            for pz_range in self._pz:
                np.logical_or(
                    self._pz_filter(
                        pz_range, trace_obj, valid_mask),
                    var_mask, out=var_mask
                )
            np.logical_and(var_mask, valid_mask, out=valid_mask)

        if self._r is not None:
            var_mask.fill(0)
            for r_range in self._r:
                np.logical_or(
                    self._r_filter(
                        r_range, trace_obj, valid_mask),
                    var_mask, out=var_mask
                )
            np.logical_and(var_mask, valid_mask, out=valid_mask)

        if self._dir is not None:
            var_mask.fill(0)
            for dir_cfg in self._dir:
                np.logical_or(
                    self._dir_filter(
                        dir_cfg, trace_obj, valid_mask),
                    var_mask, out=var_mask
                )
            np.logical_and(var_mask, valid_mask, out=valid_mask)

        if self._pl is not None and trace_obj.plon:
            for pl_range in self._pl:
                np.logical_or(
                    self._pl_filter(
                        pl_range, trace_obj, valid_mask),
                    var_mask, out=var_mask
                )
            np.logical_and(var_mask, valid_mask, out=valid_mask)

        # valid but saturated
        selected_mask = np.logical_and(valid_mask, np.logical_not(too_long_mask))
        dropped_mask = np.logical_and(valid_mask, too_long_mask)
        n_dropped = np.count_nonzero(dropped_mask)

        valid_data = trace_obj.data[selected_mask, :]
        valid_n = trace_obj.n[selected_mask]
        if not update:
            out = Trace(trace_obj)
        else:
            out = trace_obj
        out.data = valid_data
        out.n = valid_n

        return out, n_dropped

    def _x_filter(self, x_range, trace_obj, input_mask):
        x = trace_obj.terminal['x']
        mask = np.logical_and(x >= x_range[0], x <= x_range[1])
        return np.logical_and(mask, input_mask, out=mask)

    def _y_filter(self, y_range, trace_obj, input_mask):
        y = trace_obj.terminal['y']
        mask = np.logical_and(y >= y_range[0], y <= y_range[1])
        return np.logical_and(mask, input_mask, out=mask)

    def _z_filter(self, z_range, trace_obj, input_mask):
        z = trace_obj.terminal['z']
        mask = np.logical_and(z >= z_range[0], z <= z_range[1])
        return np.logical_and(mask, input_mask, out=mask)

    def _pz_filter(self, pz_range, trace_obj, input_mask):
        pz = trace_obj.terminal['pz']
        mask = np.logical_and(pz >= pz_range[0], pz <= pz_range[1])
        return np.logical_and(mask, input_mask, out=mask)

    def _r_filter(self, r_range, trace_obj, input_mask):
        if len(r_range) > 2:
            x0, y0 = r_range[2]
        else:
            x0 = y0 = 0.0
        x = trace_obj.terminal['x'] - x0
        y = trace_obj.terminal['y'] - y0
        rr = x**2 + y**2
        mask = np.logical_and(rr >= r_range[0]**2, rr <= r_range[1]**2)
        return np.logical_and(mask, input_mask, out=mask)

    def _dir_filter(self, dir_cfg, trace_obj, input_mask):
        px = trace_obj.terminal['px']
        py = trace_obj.terminal['py']
        pz = trace_obj.terminal['pz']

        cosmin, cosmax = dir_cfg[0], dir_cfg[1]
        v = dir_cfg[2]

        cos_theta = px*v[0] + py*v[1] + pz*v[2]

        mask = np.logical_and(cos_theta >= cosmin, cos_theta <= cosmax)
        return np.logical_and(mask, input_mask, out=mask)

    def _pl_filter(self, pl_range, trace_obj, input_mask):
        pl = trace_obj.terminal['pl']
        mask = np.logical_and(pl >= pl_range[0], pl <= pl_range[1])
        return np.logical_and(mask, input_mask, out=mask)

    def _chek_range(self, param, value):
        if not isinstance(value, (list, tuple)):
            raise TypeError(
                'Filter parameter {} must be a tuple of two float '\
                'values (low, high), or a list of such tuples!'.format(param))

        if isinstance(value[0], (list, tuple)):
            result = [_range(item) for item in value]
        else:
            result = (_range(value),)

        return result

    def _chek_r(self, r: tuple):
        if not isinstance(r, (list, tuple)):
            raise TypeError(
                'Filter parameter r must be a tuple '\
                '(rmin, rmax, (x_origin, y_origin)) or a list or tuple of '\
                'such tuples!')

        if isinstance(r[0], (list, tuple)):
            result = [(float(item[0]), float(item[1]), _origin2(item[2]))
                      for item in r]
        else:
            result = ((float(r[0]), float(r[1]), _origin2(r[2])),)

        return result

    def _chek_dir(self, dir: tuple):
        if not isinstance(dir, (list, tuple)):
            raise TypeError(
                'Filter parameter dir must be a tuple '\
                '(cosmin, cosmax, (px, py, pz)) or a list or tuple of '\
                'such tuples!')

        if isinstance(dir[0], (list, tuple)):
            result = []
            for item in dir:
                p = item[2]
                l = np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
                if l == 0.0:
                    raise ValueError('Direction vector length must not be 0!')
                result.append(
                    (float(item[0]), float(item[1]), _dir3(p, 1.0/l))
                )
        else:
            p = dir[2]
            l = np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
            if l == 0.0:
                raise ValueError('Direction vector length must not be 0!')
            result = ((float(dir[0]), float(dir[1]), _dir3(p, 1.0/l)),)

        return result

    def __str__(self):
        return 'Filter(x={}, y={}, z={}, pz={}, dir={}, r={}, pl={})'.format(
            self.x, self.y, self.z, self.pz, self.dir, self.r, self.pl)

    def __repr__(self):
        return self.__str__() + \
            ' # Filter object at 0x{:>08X}.'.format(id(self))


class Trace(mcobject.McObject):
    '''
    Photon packet trace interface.
    '''
    # Trace no events - turns off the photon packet trace functionality.
    TRACE_NONE = 0
    # Trace start events - immediately after the packet is launched.
    TRACE_START = 1
    # Trace end events - just before the packet simulation is terminated.
    TRACE_END = 2
    # Trace all events including start, stop and all in between.
    TRACE_ALL = 7

    # Number of fields in a single trace entry.
    TRACE_ENTRY_LEN = 8

    @staticmethod
    def cl_type(mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Structure that is passed to the Monte carlo simulator kernel.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents a trace in the Monte Carlo kernel.

            The returned structure type implements the following fields:

            - max_events: mc_int_t 
                Maximum number of trace events per photon packet.
            - data_buffer_offset: mc_size_t
                Offset of the first element of the trace data buffer.
            - count_buffer_offset: mc_size_t
                Offset of the first element of the trace event counter buffer.
        '''
        T = mc.types
        class ClTrace(cltypes.Structure):
            _fields_ = [
                ('max_events', T.mc_int_t),
                ('data_buffer_offset', T.mc_size_t),
                ('count_buffer_offset', T.mc_size_t),
            ]

        return ClTrace

    @staticmethod
    def cl_declaration(mc: mcobject.McObject) -> str:
        '''
        Structure that defines the trace in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'struct MC_STRUCT_ATTRIBUTES McTrace{',
            '	mc_int_t max_events;',
            '	mc_size_t data_buffer_offset;',
            '	mc_size_t count_buffer_offset;',
            '};',
        ))

    @staticmethod
    def cl_implementation(mc: mcobject.McObject) -> str:
        '''
        Implementation of the trace in the Monte Carlo simulator.
        '''
        return '\n'.join((
            'void dbg_print_trace(__mc_trace_mem const McTrace *trace){',
            '	dbg_print("McTrace object:");',
            '	dbg_print_int(INDENT "max_events:", trace->max_events);',
            '	dbg_print_size_t(INDENT "data_buffer_offset:", trace->data_buffer_offset);',
            '	dbg_print_size_t(INDENT "count_buffer_offset:", trace->count_buffer_offset);',
            '};',
            '',
            'inline int mcsim_trace_event(McSim *mcsim, mc_uint_t event_count){',
            '	__mc_trace_mem const McTrace *trace = mcsim_trace(mcsim);',
            '	mc_size_t pos = mc_min(event_count, trace->max_events - 1)*TRACE_ENTRY_LEN + ',
            '		mcsim_packet_index(mcsim)*trace->max_events*TRACE_ENTRY_LEN + ',
            '		trace->data_buffer_offset;',
            '',
            '	dbg_print("Tracing photon packet event:");'
            '	dbg_print_cnt(INDENT       "packet index:", mcsim_packet_index(mcsim));',
            '	dbg_print_cnt(INDENT       "event index :", event_count);',
            '	dbg_print_point3f(INDENT   "position    :", mcsim_position(mcsim));',
            '	dbg_print_point3f(INDENT   "direction   :", mcsim_direction(mcsim));',
            '	dbg_print_float(INDENT     "weight      :", mcsim_weight(mcsim));',
            '	#if MC_TRACK_OPTICAL_PATHLENGTH',
            '		dbg_print_float(INDENT "path length :", mcsim_optical_pathlength(mcsim));',
            '	#endif',
            '',
            '	mcsim_float_buffer(mcsim)[pos++] = mcsim_position_x(mcsim);',
            '	mcsim_float_buffer(mcsim)[pos++] = mcsim_position_y(mcsim);',
            '	mcsim_float_buffer(mcsim)[pos++] = mcsim_position_z(mcsim);',
            '	mcsim_float_buffer(mcsim)[pos++] = mcsim_direction_x(mcsim);',
            '	mcsim_float_buffer(mcsim)[pos++] = mcsim_direction_y(mcsim);',
            '	mcsim_float_buffer(mcsim)[pos++] = mcsim_direction_z(mcsim);',
            '	mcsim_float_buffer(mcsim)[pos++] = mcsim_weight(mcsim);',
            '	#if MC_TRACK_OPTICAL_PATHLENGTH',
            '		mcsim_float_buffer(mcsim)[pos++] = mcsim_optical_pathlength(mcsim);',
            '	#else',
            '		mcsim_float_buffer(mcsim)[pos++] = FP_0;',
            '	#endif',
            '	return 1;',
            '};',
            '',
            'inline void mcsim_trace_complete(McSim *mcsim, mc_uint_t event_count){',
            '	mcsim_integer_buffer(mcsim)[',
            '		mcsim_trace(mcsim)->count_buffer_offset + ',
            '		mcsim_packet_index(mcsim)] = (mc_int_t)event_count;',
            '	dbg_print("Trace finished");',
            '	dbg_print_uint(INDENT "total events:", event_count);',
            '	dbg_print_cnt(INDENT "packet:", mcsim_packet_index(mcsim));',
            '};',
        ))

    def cl_options(self, mc: mcobject.McObject) -> mcoptions.RawOptions:
        return [('MC_USE_TRACE', self._options),
                ('TRACE_ENTRY_LEN', int(Trace.TRACE_ENTRY_LEN)),
                ('MC_USE_SAMPLING_VOLUME', True),
                ('MC_TRACK_OPTICAL_PATHLENGTH', bool(self.plon))]

    def __init__(self, maxlen: int or 'Trace' = 500, options: int = TRACE_ALL,
                 filter: Filter = None, plon: bool = True):
        '''
        Photon packet trace object constructor.

        Parameters
        ----------
        maxlen: int or Trace
            Maximum number of trace entries (events) per photon packet.
            If maxlen is an instance of Trace, a new copy of the trace object
            is created.
            Each trace entry is composed of seven fields that describe the
            state of the photon packet at the time of the event:

            'x' - x coordinate of the photon packet

            'y' - y coordinate of the photon packet

            'z' - z coordinate of the photon packet

            'px' - x component of the photon packet propagation direction

            'py' - y component of the photon packet propagation direction

            'pz' - z component of the photon packet propagation direction

            'w'  - photon packet weight

            'pl' - photon path length

        options: int
            Trace options must be formed by combining (logical or) the
            following flags:

                TRACE_START - trace the initial photon packet state

                TRACE_END - trace the final/terminal photon packet state

                TRACE_ALL - trace all (scattering/absorption and boundary transition) events

        filter: Filter
            Filter trace results according to the supplied filter.
            A filter class needs to define a __call__ method that returns a
            new Trace instance with filtered events.
            See the Filter class for details and example.

        plon: bool
            Turn on/off tracking of the optical pathlength. If set to nonzero,
            the simulator option MC_TRACK_OPTICAL_PATHLENGTH will be set
            turning on the tracking of optical pathlength.

        Note
        ----
        The value of maxlen parameter is automatically adjusted to:

            (a) 1 - if only the TRACE_START or TRACE_END option flags are set,

            (b) 2 - if the TRACE_START and the TRACE_END option flags are set,

            (c) left unchanged if the TRACE_ALL option flag is set.

        Use the update method to convert raw simulation output buffer into
        a numpy structured array with items 'x', 'y', 'z', 'px', 'py', 'pz',
        and 'w'.
        The array shape is (nphotons, length). To access the trace of first
        photon packet use obj.data[0]

        The number of valid trace entries for each photon packet can
        be obtained by obj.n (e.g. obj.n[0] for the first photon packet).
        '''

        if isinstance(maxlen, Trace):
            trace = maxlen
            self._maxlen = trace.maxlen
            self._data = trace.data
            self._n = trace.n
            self._options = trace.options
            self._filter = trace.filter
            self._n_dropped = trace.dropped
            self._plon = trace.plon
        else:
            self._data = None
            self._n = None
            self._options = int(options) & 7
            self._set_maxlen(maxlen)
            self._filter = filter
            self._n_dropped = 0
            self._plon = bool(plon)

        self._terminal = None
        self._overflow_mask = None

    def _update_terminal(self):
        if self._data is not None:
            arange_n = np.arange(self.n.size)
            last_indices = np.minimum(self.n - 1, self.maxlen - 1)

            self._terminal = self._data[arange_n, last_indices]

    def _update_overflow_mask(self):
        if self._data is not None:
            self._overflow_mask = self.n >= self.maxlen

    def todict(self) -> dict:
        '''
        Save the trace configuration without the accumulator data to
        a dictionary.

        Returns
        -------
        data: dict
            Trace configuration as a dictionary.
        '''
        filt = None
        if self._filter is not None:
            filt = self._filter.todict()
        return {
            'type':'Trace',
            'maxlen':self._maxlen,
            'options':self._options,
            'plon': self._plon,
            'filter':filt
        }

    @classmethod
    def fromdict(cls, data: dict) -> 'Trace':
        '''
        Create a Trace instance from a dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`Trace.todict` method.
        '''
        trace_type = data.pop('type')
        if trace_type != cls.__name__:
            raise TypeError('Expected "{}" type bot got "{}"!'.format(
                cls.__name__, trace_type))
        filter_data = data.pop('filter')
        filt = None
        if filter_data is not None:
            filt = Filter.fromdict(filter_data)

        return cls(filter=filt, **data)

    def dtype(self, mc: mcobject.McObject) -> np.dtype:
        '''
        Return a structure numpy data type that represents the trace data.

        Parameters
        ----------
        mc: mcobject.McObject
            Instance of the Monte Carlo simulator which will be used with this
            trace object.
        '''
        T = mc.types
        return np.dtype([
            ('x', T.np_float),
            ('y', T.np_float),
            ('z', T.np_float),
            ('px', T.np_float),
            ('py', T.np_float),
            ('pz', T.np_float),
            ('w', T.np_float),
            ('pl', T.np_float)
        ])

    def float_len(self, nphotons: int) -> int:
        '''
        Number of float values required to save trace data for the given number
        of photon packets.

        Parameters
        ----------
        nphotons: int
            Number of photon packets to trace
        '''
        return Trace.TRACE_ENTRY_LEN*self.maxlen*nphotons

    def _get_terminal(self) -> np.ndarray:
        if self._terminal is None:
            self._update_terminal()
        return self._terminal
    terminal = property(_get_terminal, None, None,
                        'Terminal/final state of the photon packets.')

    def _get_overflow(self) -> np.ndarray:
        if self._overflow_mask is None:
            self._update_overflow_mask()
        return self._overflow_mask
    overflow = property(_get_overflow, None, None,
                        'Binary mask of photon packets that have overflowed '
                        'the trace')

    def _get_filter(self) -> Filter:
        return self._filter
    def _set_filter(self, filter: Filter):
        if not isinstance(filter, Filter):
            raise TypeError('Expected a "Filter" instance but '
                            'got "{}"'.format(type(filter)))
        self._filter = filter
    filter = property(_get_filter, _set_filter, None, 'Trace filter.')

    def _get_plon(self) -> bool:
        return self._plon
    plon = property(_get_plon, None, None,
                    'State of the optical pathlength tracking.')


    def _get_data(self) -> np.ndarray:
        return self._data
    def _set_data(self, data: np.ndarray):
        self._data = data
    data = property(_get_data, _set_data, None, 'Trace data.')

    def _get_n(self) -> int:
        return self._n
    def _set_n(self, n: int):
        self._n = n
    n = property(_get_n, _set_n, None,
                 'Number of trace events for each photon packet.')

    def _get_options(self) -> int:
        return self._options
    options = property(_get_options, None, None, 'Compile time trace options.')

    def _get_maxlen(self) -> int:
        return self._maxlen
    def _set_maxlen(self, maxlen: int):
        if self._options == Trace.TRACE_ALL:
            pass
        elif self._options == (Trace.TRACE_START | Trace.TRACE_END):
            maxlen = 2
        else:
            maxlen = 1
        if maxlen < 1:
            raise ValueError('Maximum trace length must be at least 1!')
        self._maxlen = int(maxlen)
    maxlen = property(_get_maxlen, _set_maxlen, None, 'Maximum trace length.')

    def _get_nphotons(self):
        n = 0
        if self._n is not None:
            n = self._n.size
        return n
    nphotons = property(_get_nphotons, None, None,
                        'Number of photon packets in the trace.')

    def np_buffer(self, mc: mcobject.McObject, allocation: BufferAllocation,
                  nphotons=None, **kwargs) -> np.ndarray:
        '''
        Efficient numpy buffer allocation that can be used to directly download
        data to the Trace object.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that produced the results.
        allocation: 
            Buffer allocation made by the :py:meth:`Trace.cl_pack` method.
        nphotons: int
            The number of photon packets traced by the simulation.

        Return
        ------
        np_buffer: np.ndarray
            A numpy buffer that can be used to directly download the data
            (no copy is required in the :py:meth:`Trace.update_data` method).
        '''
        dtype_int = np.dtype(mc.types.np_int)
        dtype_float = np.dtype(mc.types.np_float)
        if allocation.dtype == dtype_int:
            return np.empty((nphotons,), dtype=dtype_int)
        elif allocation.dtype == dtype_float:
            return np.empty((nphotons, self._maxlen), dtype=self.dtype(mc))
        else:
            raise RuntimeError('Unexpected buffer allocation!')

    def n_to_cl(self, mc: mcobject.McObject) -> np.ndarray:
        '''
        Prepare the array with number of events per photon packet for
        upload to the OpenCL kernel.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that will use the data.

        Returns
        -------
        cl_data: np.ndarray
            Flat numpy data array of basic type ready for transfer to the
            OpenCL device.
        '''
        return np.asarray(self._n, dtype=mc.types.np_int)

    def data_to_cl(self, mc: mcobject.McObject) -> np.ndarray:
        '''
        Prepare the trace data array for upload to the OpenCL kernel.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that will use the data.

        Returns
        -------
        cl_data: np.ndarray
            Flat numpy data array of basic type ready for transfer to the
            OpenCL device.
        '''
        cl_dtype = self.dtype(mc)
        #if self._data.dtype == cl_dtype:
        #    return np.frombuffer(self._data, dtype=mc.types.np_float)
        #else:
        return np.frombuffer(np.asarray(self._data, dtype=cl_dtype),
                             dtype=mc.types.np_float)

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None,
                nphotons: int = None) -> cltypes.Structure:
        '''
        Fills the ctypes structure (target) with the data required by the
        Monte Carlo simulator. See the Trace.cl_type class for a detailed
        list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Ctypes structure that is filled with the source data.
        nphotons: int
            The number of photon packets that will be run by the simulator.

        Returns
        -------
        target: cltypes.Structure
            Filled ctypes structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()
        if nphotons is None:
            raise ValueError('The number of photon packets was not defined!')

        # Allocate a data buffer for the trace events.
        allocation = mc.cl_allocate_rw_float_buffer(
            self, (Trace.TRACE_ENTRY_LEN*self.maxlen*nphotons,))
        target.data_buffer_offset = allocation.offset

        # Allocate data buffer for holding the number of trace events.
        allocation = mc.cl_allocate_rw_int_buffer(self, (nphotons, ))
        target.count_buffer_offset = allocation.offset

        target.max_events = self.maxlen

        return target

    def apply_filter(self):
        if self._filter is not None:
            self._n_dropped = self._filter(self, update=True)[1]

    def _get_dropped(self) -> int:
        return self._n_dropped
    dropped = property(
        _get_dropped, None, None,
        'The number of photon packets that matched the filter, but were '\
        'dropped due to exceeding the maximum trace length.')

    def plot(self, inds: np.ndarray = None, view='3d', show: bool = True):
        '''
        Plots the paths of photon packets.

        Parameters
        ----------
        inds: np.ndarray
            A vector of indices or a logical mask for selecting the photon
            packet. If a logical mask, the shape of the mask should be
            (num_packets).
        view: str
            A string that specifies the plot view:

            * 'xy' - projection into the x-y plane.
            * 'xy' - projection into the x-z plane.
            * 'xy' - projection into the y-z plane.
            * '3d' - 3D plot (default).
        show: bool
            If True, show the plot.
        '''
        from matplotlib import pyplot as pp
        from mpl_toolkits.mplot3d import Axes3D

        if inds is None:
            inds = np.arange(self._data.shape[0])
        else:
            inds = np.asarray(inds)
            if inds.dtype == np.bool:
                inds, = np.nonzeros(inds)

        fig = None

        if view == '3d':
            fig = pp.figure()
            ax = Axes3D(fig)

            for photon in inds:
                n = min(self._n[photon], self.maxlen)
                ax.plot(self._data[photon]['x'][:n],
                        self._data[photon]['y'][:n],
                        self._data[photon]['z'][:n])
                ax.plot([self._data[photon]['x'][0]],
                        [self._data[photon]['y'][0]],
                        [self._data[photon]['z'][0]],
                        marker='.', color=(0, 1, 0))
                ax.plot([self._data[photon]['x'][n - 1]],
                        [self._data[photon]['y'][n - 1]],
                        [self._data[photon]['z'][n - 1]],
                        marker='.', color=(1, 0, 0))

            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

        elif view in('xy', 'xz', 'yz'):
            fig, ax = pp.subplots()
            for photon in inds:
                n = min(self._n[photon], self.maxlen)
                ax.plot(self._data[photon][view[0]][:n],
                        self._data[photon][view[1]][:n])
                ax.plot([self._data[photon][view[0]][0]],
                        [self._data[photon][view[1]][0]],
                        marker='.', color=(0, 1, 0))
                ax.plot([self._data[photon][view[0]][n - 1]],
                        [self._data[photon][view[1]][n - 1]],
                        marker='.', color=(1, 0, 0))

            ax.set_xlabel(view[0])
            ax.set_ylabel(view[1])

        if fig is not None:
            fig.canvas.manager.set_window_title('Trace view')
            pp.tight_layout()
            if show:
                pp.show()

    def update(self, obj: 'Trace'):
        '''
        Update the trace with data from the given trace object.

        Parameters
        ----------
        obj: Trace
            Update the trace of this instance with the data
            from Trace instance obj.
        '''
        if self.maxlen != obj.maxlen:
            raise ValueError('Cannot update the trace with trace data '
                             'of different maximum length!')

        self._terminal = None
        self._overflow_mask = None

        if self._data is not None:
            if self._data.dtype != obj.data.dtype:
                raise TypeError('Cannot update the trace with trace data of '
                                'different type!')

            if self._filter is not None:
                obj = self.filter(obj, update=False)

            self._n = np.stack([self._n, obj.n])
            self._data = np.vstack([self._data, obj.data])
        else:
            self._data = obj.data
            self._n = obj.n

    def update_data(self, mc: mcobject.McObject,
                    data: Dict[np.dtype, List[np.ndarray]],
                    nphotons: int, **kwargs):
        '''
        Update the content of raw data accumulators with simulation results.

        Parameters
        ----------
        mc: mcobject.McObject
            Simulator instance that produced the data.
        data: Dict[np.dtype, List[np.ndarray]]
            Dict of allocated buffers (this implementation
            uses only one integer and one floating point buffer).
        nphotons: int
            Number of photon packets used in the simulations.
        kwargs: dict
            Additional keyword arguments not used by this implementation.
        '''
        #int_buffers = data[np.dtype(mc.types.np_int)]
        #self._n = np.copy(int_buffers[0])

        #float_buffers = data[np.dtype(mc.types.np_float)]
        #self._data = np.frombuffer(
        #    np.copy(float_buffers[0]), dtype=self.dtype(mc))
        #self._data.shape = (nphotons, self._maxlen)
        new_data_dtype = np.dtype(mc.types.np_float)
        new_data = data[new_data_dtype][0]
        new_n = data[np.dtype(mc.types.np_int)][0]

        self._terminal = None
        self._overflow_mask = None

        if self._data is not None:
            if self._data.dtype != new_data.dtype:
                raise TypeError('Cannot update the trace with trace data of '
                                'different type!')
            if self._data.shape[1] != new_data.shape[1]:
                raise ValueError('Cannot update the trace with trace data '
                                 'of different maximum length!')

            tmp_trace = type(self)(self)
            tmp_trace.data = None
            tmp_trace.n = None
            tmp_trace.update_data(mc, data, nphotons)

            self._n = np.hstack([self._n, tmp_trace.n])
            self._data = np.vstack([self._data, tmp_trace._data])
        else:
            self._n = new_n
            self._data = new_data
            # apply filter to the trace data
            self.apply_filter()


    def __len__(self) -> int:
        '''
        Returns the number of photon packets in the trace.

        Returns
        -------
        n: int
            The number of photon packets in the trace.
        '''
        if self._n is not None:
            return self._n.size

        return 0

    def __str__(self):
        if self._options:
            if self._options == Trace.TRACE_ALL:
                flagstr = 'Trace.TRACE_ALL'
            else:
                flagstr = ''
                if self._options & Trace.TRACE_START:
                    flagstr += 'Trace.TRACE_START | '
                if self._options & Trace.TRACE_END:
                    flagstr += 'Trace.TRACE_END | '
                flagstr = flagstr[:-3]
        else:
            flagstr = 'Trace.TRACE_NONE'
        return 'Trace(maxlen={}, options={}, filter={}, plon={})'.format(
            self.maxlen, flagstr, self.filter, self.plon)

    def __repr__(self):
        return self.__str__() + \
            ' # Trace object at 0x{:>08X}.'.format(id(self))
