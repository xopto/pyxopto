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

from xopto.mccyl import mcoptions
from xopto.mccyl import mcobject
from xopto.mccyl import cltypes

OUTER = 'outer'
''' Identifier of the outer sample surface. '''

INTERNAL = 'internal'
''' Identifier of the internal sample surface. '''

NONE = 'none'
''' Identifier used for a surface that is not yet assigned. '''

class SurfaceLayoutBase(mcobject.McObject):
    '''
    Base class of all the sample surface layouts.
    '''
    def __init__(self, location: str = OUTER or INTERNAL):
        '''
        Complex surface layout at the specified sample surface.

        Parameters
        ----------
        location: OUTER or INTERNAL
            Location of the sample surface to which the layout applies.
        '''
        super().__init__()
        if location not in (OUTER, INTERNAL, NONE):
            raise ValueError(
                'Surface layout location must be "{}", "{}" or "{}"!'.format(
                    NONE, OUTER, INTERNAL)
            )
        self._location = location

    def _get_location(self) -> str:
        return self._location
    def _set_location(self, location: OUTER or INTERNAL):
        if location not in (NONE, OUTER, INTERNAL):
            raise ValueError(
                'Surface layout location must be "{}", "{}" or "{}"!'.format(
                    NONE, OUTER, INTERNAL)
            )
        if location != self._location and self._location != NONE:
            raise RuntimeError('Surface layout cannot be changed!')
        self._location = str(location)
    location = property(_get_location, _set_location, None,
                        'Location of the surface layout.')


class SurfaceLayoutOuter(SurfaceLayoutBase):
    '''
    Base class of the outer sample surface layout.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create a layout for the  outer sample surface.
        '''
        super().__init__(OUTER)


class SurfaceLayoutInternal(SurfaceLayoutBase):
    '''
    Base class of the internal sample surfaces layout.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create a layout for the internal sample surfaces.
        '''
        super().__init__(INTERNAL)


class SurfaceLayoutAny(SurfaceLayoutBase):
    '''
    Base class of surface layouts that can be used for both the
    outer or internal sample surfaces.
    '''
    def __init__(self, location: str = NONE):
        '''
        Create a layout for the outer or internal sample surfaces. By default
        the location is not initialized.
        '''
        super().__init__(location)

    @classmethod
    def fromdict(cls, data:dict) -> 'SurfaceLayoutAny':
        '''
        Create a new instance of a surface layout from a dict.
        The dict keys must match the parameter names defined by the
        constructor.
        '''
        data_ = dict(data)
        T = data_.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Cannot initialize an instance of "{:s}" from the data'
                'of instance "{}"!'.format(cls.__name__, T))

        return cls(**data_)

class SurfaceLayoutDefault(SurfaceLayoutAny):
    '''
    Default / dummy layout for the outer or internal sample surface.
    '''
    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Structure that is passed to the Monte carlo simulator kernel.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents the surface layout in
            the Monte Carlo kernel.

            The returned structure type implements the following fields:

            - dummy: mc_int_t
                Dummy field of the dummy detector.
        '''
        T = mc.types
        class ClSurfaceLayoutDefault(cltypes.Structure):
            _fields_ = [('dummy', T.mc_int_t)]
        return ClSurfaceLayoutDefault

    @staticmethod
    def cl_options(mc: mcobject.McObject) -> str:
        '''
        Raw OpenCL options of the default detector.
        '''
        return ''

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Structure that defines the dummy reflector in the Monte Carlo simulator.
        '''
        Loc = self.location.capitalize()
        return 'struct Mc{}SurfaceLayout{{mc_int_t dummy;}};'.format(Loc)

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        OpenCL implementation of the default surface layout
        '''
        loc = self.location
        Loc = loc.capitalize()
        return '\n'.join((
            'void dbg_print_{}_surface_layout('.format(loc),
            '		__mc_surface_mem const Mc{}SurfaceLayout *layout){{'.format(Loc),
            '	dbg_print("Mc{}SurfaceLayout - default layout:");'.format(Loc),
            '};',
        ))

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`SurfaceLayoutDefault.cl_type` for a detailed
        list of fields.

        Parameters
        ----------
        mc: mcobject.McObject
            Monte Carlo simulator instance.
        target: cltypes.Structure
            Structure that is filled with the source data.

        Returns
        -------
        target: cltypes.Structure
            Filled ctypes structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        target.dummy = 0

        return target

    def todict(self) -> dict:
        '''
        Export object to dict.
        '''
        return {'type': self.__class__.__name__}

    @classmethod
    def fromdict(cls, data: dict) -> 'SurfaceLayoutDefault':
        '''
        Create an instance of :py:class:`SurfaceLayoutDefault` from a
        dictionary.

        Parameters
        ----------
        data: dict
            Dictionary created by the :py:meth:`SurfaceLayoutDefault.todict`
            method.
        '''
        data_ = dict(data)
        layout_type = data_.pop('type')
        if layout_type != cls.__name__:
            raise TypeError(
                'Expected a "{}" type bot got "{}"!'.format(
                    cls.__name__, layout_type))

        return cls(**data_)


class SurfaceLayouts(mcobject.McObject):
    def cl_type(self, mc: mcobject.McObject) -> cltypes.Structure:
        '''
        Structure that is passed to the OpenCL simulator.

        Parameters
        ----------
        mc: McObject
            A Monte Carlo simulator instance.

        Returns
        -------
        struct: cltypes.Structure
            A structure type that represents surface layouts in the
            Monte Carlo kernel.

            The returned structure type implements the following fields:

            - outer: SurfaceLayoutOuter
                Layout of the outer sample surface (z = 0).
            - internal: SurfaceLayoutInternal
                Layout of the internal sample surfaces (z = 0).
        '''
        class ClSurfaceLayouts(cltypes.Structure):
            _fields_ = [
                ('outer',    self.outer.fetch_cl_type(mc)),
                ('internal', self.internal.fetch_cl_type(mc))
            ]
        return ClSurfaceLayouts

    def cl_options(self, mc: mcobject.McObject) -> mcoptions.RawOptions:
        '''
        Returns the OpenCL options of the surface layout.

        If the outer and / or internal sample surface layouts are specified
        (not the dummy default :py:class:`SurfaceLayoutDefault`),
        the corresponding OpenCL options that activate the use of surface
        layouts are set, i.e. MC_USE_OUTER_SURFACE_LAYOUT for the outer surface
        and MC_USE_INTERNAL_SURFACE_LAYOUT for the internal surface.
        '''
        out = []
        use_layouts = False
        if type(self.outer) != SurfaceLayoutDefault:
            out.append(('MC_USE_OUTER_SURFACE_LAYOUT', True))
            out.extend(self._outer.fetch_cl_options(mc))
            use_layouts = True
        if type(self.internal) != SurfaceLayoutDefault:
            out.append(('MC_USE_INTERNAL_SURFACE_LAYOUT', True))
            out.extend(self._internal.fetch_cl_options(mc))
            use_layouts = True
        if use_layouts:
            out.append(('MC_USE_SURFACE_LAYOUTS', True))
        return out

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Declarations of surface layouts in OpenCL.
        '''
        return '\n'.join((
            self._outer.fetch_cl_declaration(mc),
            'typedef struct McOuterSurfaceLayout McOuterSurfaceLayout;',
            '',
            self._internal.fetch_cl_declaration(mc),
            'typedef struct McInternalSurfaceLayout McInternalSurfaceLayout;',
            '',
            'struct MC_STRUCT_ATTRIBUTES McSurfaceLayouts{',
            '	McOuterSurfaceLayout outer;',
            '	McInternalSurfaceLayout internal;',
            '};',
            'typedef struct McSurfaceLayouts McSurfaceLayouts;',
            '',
            '#define mcsim_outer_surface_layout(psim) (&mcsim_surface_layouts(psim)->outer)',
            '#define mcsim_internal_surface_layout(psim) (&mcsim_surface_layouts(psim)->internal)',
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of surface layouts.
        '''
        return '\n'.join((
            self._outer.fetch_cl_implementation(mc),
            '',
            self._internal.fetch_cl_implementation(mc),
            '',
            'void dbg_print_surface_layouts(__mc_surface_mem const McSurfaceLayouts *layouts){',
            '	dbg_print("McSurfaceLayouts");',
            '	dbg_print_outer_surface_layout(&layouts->outer);',
            '	dbg_print_internal_surface_layout(&layouts->internal);',
            '};',
        ))


    def __init__(self, outer: SurfaceLayoutOuter or SurfaceLayoutAny = None,
                 internal: SurfaceLayoutInternal or SurfaceLayoutAny = None):
        '''
        Container of the outer and internal sample surface layouts.

        Parameters
        ----------
        outer: SurfaceLayoutOuter or SurfaceLayoutAny
            Layout of the outer sample surface.
        internal: SurfaceLayoutInternal or SurfaceLayoutAny
            Layout at the internal sample surfaces.
        '''
        if isinstance(outer, SurfaceLayouts):
            sg = outer
            outer = sg.outer
            internal = sg.internal

        if outer is None:
            outer = SurfaceLayoutDefault()
        elif not isinstance(outer, (SurfaceLayoutOuter,
                                  SurfaceLayoutAny)):
            raise TypeError('Layout of the outer sample surface must be an '
                            'instance of SurfaceLayoutOuter or '
                            'SurfaceLayoutAny!')

        if internal is None:
            internal = SurfaceLayoutDefault()
        elif not isinstance(internal, (SurfaceLayoutInternal,
                                     SurfaceLayoutAny)):
            raise TypeError('Layout of the internal sample surfaces must be an '
                            'instance of SurfaceLayoutInternal or '
                            'SurfaceLayoutAny!')

        outer.location = OUTER
        internal.location = INTERNAL

        self._outer = outer
        self._internal = internal

    def _get_outer(self) -> mcobject.McObject:
        return self._outer
    outer = property(_get_outer, None, None, 'Outer sample surface layout.')

    def _get_internal(self) -> mcobject.McObject:
        return self._internal
    internal = property(_get_internal, None, None,
                        'Layouts of the internal sample surfaces.')

    def cl_pack(self, mc: mcobject.McObject, target: cltypes.Structure = None) \
            -> cltypes.Structure:
        '''
        Fills the structure (target) with the data required by the
        Monte Carlo simulator.
        See the :py:meth:`SurfaceLayouts.cl_type` method for a detailed
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
            Filled structure received as an input argument or a new
            instance if the input argument target is None.
        '''
        if target is None:
            target_type = self.cl_type(mc)
            target = target_type()

        self.outer.cl_pack(mc, target.outer)
        self.internal.cl_pack(mc, target.internal)

        return target

    def types(self) -> tuple:
        '''
        Returns a tuple of surface layout types assigned to this instance.
        '''
        return type(self._outer), type(self._internal)

    def todict(self) -> dict:
        '''
        Export object to dict.
        '''

        return {
            'type': self.__class__.__name__,
            'outer': self._outer.todict(),
            'internal': self._internal.todict()
        }

    def __str__(self):
        return 'SurfaceLayouts(outer={}, internal={})'.format(
            self._outer, self._internal)

    def __repr__(self):
        return '{} #{}'.format(self.__str__(), id(self))
