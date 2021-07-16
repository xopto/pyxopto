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

from xopto.mcvox import mcoptions
from xopto.mcvox import mcobject
from xopto.mcvox import cltypes


TOP = 'top'
''' Identifier of the top sample surface. '''

BOTTOM = 'bottom'
''' Identifier of the bottom sample surface. '''

NONE = 'none'
''' Identifier used for a surface that is not yet assigned. '''

class SurfaceLayoutBase(mcobject.McObject):
    '''
    Base class of all the sample surface layouts.
    '''
    def __init__(self, location: str = TOP or BOTTOM):
        '''
        Complex surface layout at the specified sample surface.

        Parameters
        ----------
        location: TOP or BOTTOM
            Location of the sample surface to which the layout applies.
        '''
        super().__init__()
        if location not in (TOP, BOTTOM, NONE):
            raise ValueError(
                'Surface layout location must be "{}", "{}" or "{}"!'.format(
                    NONE, TOP, BOTTOM)
            )
        self._location = location

    def _get_location(self) -> str:
        return self._location
    def _set_location(self, location: TOP or BOTTOM):
        if location not in (NONE, TOP, BOTTOM):
            raise ValueError(
                'Surface layout location must be "{}", "{}" or "{}"!'.format(
                    NONE, TOP, BOTTOM)
            )
        if location != self._location and self._location != NONE:
            raise RuntimeError('Surface layout cannot be changed!')
        self._location = str(location)
    location = property(_get_location, _set_location, None,
                        'Location of the surface layout.')


class SurfaceLayoutTop(SurfaceLayoutBase):
    '''
    Base class of the top sample surface layout.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create a layout for the top sample surface.
        '''
        super().__init__(TOP)


class SurfaceLayoutBottom(SurfaceLayoutBase):
    '''
    Base class of the bottom sample surface layout.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create a layout for the bottom sample surface.
        '''
        super().__init__(BOTTOM)


class SurfaceLayoutAny(SurfaceLayoutBase):
    '''
    Base class of surface layouts that can be used for both the
    top or bottom sample surface.
    '''
    def __init__(self, location: str = NONE):
        '''
        Create a layout for the top or bottom sample surface. By default
        the location is not initialized.
        '''
        super().__init__(location)


class SurfaceLayoutDefault(SurfaceLayoutAny):
    '''
    Default / dummy layout for the top or bottom sample surface.
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

        Fields
        ------
        dummy: mc_int_t
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
            '		__mc_detector_mem const Mc{}SurfaceLayout *layout){{'.format(Loc),
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

            - top: SurfaceLayoutTop
                Layout of the top sample surface (z = 0).
            - bottom: SurfaceLayoutBottom
                Layout of the bottom sample surface (z = 0).
        '''
        class ClSurfaceLayouts(cltypes.Structure):
            _fields_ = [
                ('top',    self.top.fetch_cl_type(mc)),
                ('bottom', self.bottom.fetch_cl_type(mc))
            ]
        return ClSurfaceLayouts

    def cl_options(self, mc: mcobject.McObject) -> mcoptions.RawOptions:
        '''
        Returns the OpenCL options of the surface layout.

        If the top and / or bottom sample surface layouts are specified
        (not the dummy default layout:py:class:`SurfaceLayoutDefault`),
        the corresponding OpenCL options that activate the use of surface
        layouts are set, i.e. MC_USE_TOP_SURFACE_LAYOUT for the top surface
        and MC_USE_BOTTOM_SURFACE_LAYOUT for the bottom surface.
        '''
        out = []
        use_layouts = False
        if type(self.top) != SurfaceLayoutDefault:
            out.append(('MC_USE_TOP_SURFACE_LAYOUT', True))
            out.extend(self._top.fetch_cl_options(mc))
            use_layouts = True
        if type(self.bottom) != SurfaceLayoutDefault:
            out.append(('MC_USE_BOTTOM_SURFACE_LAYOUT', True))
            out.extend(self._bottom.fetch_cl_options(mc))
            use_layouts = True
        if use_layouts:
            out.append(('MC_USE_SURFACE_LAYOUTS', True))
        return out

    def cl_declaration(self, mc: mcobject.McObject) -> str:
        '''
        Declarations of surface layouts in OpenCL.
        '''
        return '\n'.join((
            self._top.fetch_cl_declaration(mc),
            'typedef struct McTopSurfaceLayout McTopSurfaceLayout;',
            '',
            self._bottom.fetch_cl_declaration(mc),
            'typedef struct McBottomSurfaceLayout McBottomSurfaceLayout;',
            '',
            'struct MC_STRUCT_ATTRIBUTES McSurfaceLayouts{',
            '	McTopSurfaceLayout top;',
            '	McBottomSurfaceLayout bottom;',
            '};',
            'typedef struct McSurfaceLayouts McSurfaceLayouts;',
            '',
            '#define mcsim_top_surface_layout(psim) (&mcsim_surface_layouts(psim)->top)',
            '#define mcsim_bottom_surface_layout(psim) (&mcsim_surface_layouts(psim)->bottom)',
        ))

    def cl_implementation(self, mc: mcobject.McObject) -> str:
        '''
        Implementation of surface layouts.
        '''
        return '\n'.join((
            self._top.fetch_cl_implementation(mc),
            '',
            self._bottom.fetch_cl_implementation(mc),
            '',
            'void dbg_print_surface_layouts(__mc_detector_mem const McSurfaceLayouts *layouts){',
            '	dbg_print("McSurfaceLayouts");',
            '	dbg_print_top_surface_layout(&layouts->top);',
            '	dbg_print_bottom_surface_layout(&layouts->bottom);',
            '};',
        ))


    def __init__(self, top: SurfaceLayoutTop or SurfaceLayoutAny = None,
                       bottom: SurfaceLayoutBottom or SurfaceLayoutAny = None):
        '''
        Container of the top and bottom sample surface layouts.

        Parameters
        ----------
        top: SurfaceLayoutTop or SurfaceLayoutAny
            Layout of the top sample surface.
        bottom: SurfaceLayoutBottom or SurfaceLayoutAny
            Layout at the bottom sample surface.
        '''
        if isinstance(top, SurfaceLayouts):
            sg = top
            top = sg.top
            bottom = sg.bottom

        if top is None:
            top = SurfaceLayoutDefault()
        elif not isinstance(top, (SurfaceLayoutTop,
                                  SurfaceLayoutAny)):
            raise TypeError('Layout of the top sample surface must be an '
                            'instance of SurfaceLayoutTop or '
                            'SurfaceLayoutAny!')

        if bottom is None:
            bottom = SurfaceLayoutDefault()
        elif not isinstance(bottom, (SurfaceLayoutBottom,
                                     SurfaceLayoutAny)):
            raise TypeError('Layout of the bottom sample surface must be an '
                            'instance of SurfaceLayoutBottom or '
                            'SurfaceLayoutAny!')

        top.location = TOP
        bottom.location = BOTTOM

        self._top = top
        self._bottom = bottom

    def _get_top(self) -> mcobject.McObject:
        return self._top
    top = property(_get_top, None, None, 'Top sample surface layout.')

    def _get_bottom(self) -> mcobject.McObject:
        return self._bottom
    bottom = property(_get_bottom, None, None, 'Bottom sample surface layout.')

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

        self.top.cl_pack(mc, target.top)
        self.bottom.cl_pack(mc, target.bottom)

        return target

    def types(self) -> tuple:
        '''
        Returns a tuple of surface layout types assigned to this instance.
        '''
        return type(self._top), type(self._bottom)
