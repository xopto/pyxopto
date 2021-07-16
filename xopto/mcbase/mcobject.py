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

from xopto.cl import cltypes

RawOption = Tuple[str, str or bool or int or float]
RawOptions = List[RawOption]

class McObject:
    '''
    Base class of all Monte Carlo objects that interface with the OpenCL
    simulator kernel.
    '''

    cl_type = None
    '''
    Definition of OpenCL data type (a basic cltypes or a Structure type)
    exposed by the class.
    This attribute can be a str or a callable that takes a simulator instance
    as the first input argument and returns a basic cltypes object or a
    Structure type.
    (can be a regular method, static method or a class method).
    Default implementation (None) - no exposed OpenCL types.
    '''

    cl_options = None
    '''
    Definition of one or simulator options exposed by the class.
    This attribute can be a list of options defined as
    [(name, value), (name, value), ...] or a callable that takes a simulator
    instance as the first input argument and returns a list of options
    (can be a regular method, static method or a class method).
    The options will be coverted to defines ("#define <name> <value>")
    Default implementation (None) - no exposed OpenCL options.
    '''

    cl_declaration = None
    '''
    Declaration of one or more OpenCL data types exposed by the class.
    This can be a str of valid OpenCL code or a callable that takes a simulator
    instance as the first input argument and returns a str
    (can be a regular method, static method or a class method).
    Default implementation (None) - no exposed OpenCL declarations.
    '''

    cl_implementation = None
    '''
    Implementation of one or more OpenCL functions exposed by
    the class.
    This can be a str of valid OpenCL code or a callable that takes a simulator
    instance as the first input argument and returns a str
    (can be a regular method, static method or a class method).
    Default implementation (None) - no exposed OpenCL implementation.
    '''

    def fetch_cl_options(self, mc: 'McObject') -> List[RawOption]:
        '''
        Fetch OpenCL options of the object. 
        OpenCl options of an object are defined through the cl_options
        attribute that can be a callable or a list of options defined as
        [(name, value), ...].
        A callable attribute will be called with a simulator instance
        as the first input argument and should return a list of options.

        Parameters
        ----------
        mc: McObject
            Simulator instance.

        Returns
        -------
        options: list[RawOption]
            A list of OpenCl options defined as [(name, value), ...].
        '''
        cl_options = getattr(self, 'cl_options', [])
        if callable(cl_options):
            cl_options = cl_options(mc)
        if cl_options is None:
            cl_options = []
        return cl_options

    def fetch_cl_type(self, mc: 'McObject') -> cltypes.Structure:
        '''
        Fetch the OpenCL type representation of this object.
        The type representation of an object is defined through the cl_type
        attribute that can be a basic cltypes or Structure, or a callable.
        A callable attribute will be called with a simulator instance
        as the first input argument and should return a data type.

        Parameters
        ----------
        mc: McObject
            Simulator instance.

        Returns
        -------
        cl_type: Structure or basic ctypes object
            OpenCL type representing this object in the OpenCL kernel.
        '''
        result = getattr(self, 'cl_type', '')
        if callable(result):
            result = result(mc)
        return result

    def fetch_cl_declaration(self, mc: 'McObject') -> str:
        '''
        Fetch OpenCL declarations of the object.
        OpenCL declaration is defined through the cl_declaration attribute that
        can be a callable or a string.
        A callable attribute will be called with a target simulator instance as
        the first input argument and should return a str (valid OpenCL code).

        Parameters
        ----------
        mc: McObject
            Simulator instance.

        Returns
        -------
        declaration: str 
            OpenCl declarations.
        '''
        result = getattr(self, 'cl_declaration', None)
        if callable(result):
            result = result(mc)
        return result

    def fetch_cl_implementation(self, mc: 'McObject') -> str:
        '''
        Fetch OpenCL implementation of the object. 
        OpenCl implementation is defined through the cl_declaration attribute 
        that can be a callable or a string.
        A callable attribute will be called with a target simulator instance as
        the first input argument and should return a str (valid OpenCL code).

        Parameters
        ----------
        mc: McObject
            Simulator instance.

        Returns
        -------
        implementation: str 
            OpenCl implementation.
        '''
        result = getattr(self, 'cl_implementation', None)
        if callable(result):
            result = result(mc)
        return result
