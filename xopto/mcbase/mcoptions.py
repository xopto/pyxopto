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

from .mcobject import McObject, RawOption, RawOptions

def make_defines(options: List[RawOption], indent: str=None) -> str:
    '''
    Create defines for the given list of options. These defines should be
    included at the beginning of the OpenCL source file.

    Parameters
    ----------
    options: List[RawOption]
        List of options defined as [(option_name, option_value), ...].
    indent: str
        Indent to use when creating define directives. 

    Returns
    -------
    defs: str
        Options as defines, one per line.
    '''
    if indent is None:
        indent = ''
    d = []
    names = []
    if options:
        for name, value in options:
            d.append(McOption.make_define(name, value))

    return '\n'.join(d)

def resolve_cl_options(*args: Tuple[RawOption] or List[RawOption]) \
        -> List[RawOption]:
    '''
    Resolve collections of OpenCL options passed as input arguments.
    Check for duplicate options with different values and raise
    ValueError if the same option is defined with different values.

    Parameters
    ----------
    args: Tuple[RawOption]
        Each positional argument must be a list of options defined as
        [(name, value), ...].

    Returns
    -------
    options: List[RawOption]
        List of resolved options defined as [(name, value), ...]. 
    '''
    options = []
    for item in args:
        options.extend(item)

    options = list(set(options))

    names = []
    for name, value in options:
        if name in names:
            value_ = options[names.index(name)][1]
            raise ValueError(
                'The same option cannot be defned with different values!'
                'Found multiple definitions of option {} with values '
                '{} and {}!'.format(name, value_, value))
        names.append(name)

    return options

def cl_value(value: str or bool or int or float) -> str:
    '''
    OpenCL representation of the option value that can be set by a define.
    Parameters
    ----------
    value: str or bool or int or float
        Python representation of the option value.

    Returns
    -------
    cl_value: str
        Representation of the option value that can be used in a
        standard C/OpenCL define.
    '''
    if isinstance(value, bool):
        value = {False:'FALSE', True:'TRUE'}.get(value)
    elif isinstance(value, int):
        value = '{:d}'.format(value)
    elif isinstance(value, float):
        value = '{:.16g}'.format(value)
    elif isinstance(value, str):
        pass
    elif value is None:
        value = ''
    else:
        raise TypeError('Option must be of str, bool, int or float type!')

    return value

class McOption(McObject):
    '''
    Base Class of Monte Carlo kernel options that are set up through defines.
    '''

    @staticmethod
    def cl_value(value: str or bool or int or float) -> str:
        '''
        OpenCL representation of the option value that can be set by a define.

        Parameters
        ----------
        value: str or bool or int or float
            Python representation of the option value.

        Returns
        -------
        cl_value: str
            Representation of the option value that can be used in a
            standard C/OpenCL define.
        '''
        if isinstance(value, bool):
            value = {False:'FALSE', True:'TRUE'}.get(value)
        elif isinstance(value, int):
            value = '{:d}'.format(value)
        elif isinstance(value, float):
            value = 'FP_LITERAL({:.16g})'.format(value)
        elif isinstance(value, str):
            pass
        elif value is None:
            value = ''
        else:
            raise TypeError('Option must be of str, bool, int or float type!')

        return value

    @classmethod
    def make_define(cls, name, value) -> str:
        '''
        Create and return a standard C/OpenCL define.

        Returns
        -------
        define: str
            A standard C/OpenCL define representing the option.
        '''
        if value is None:
            return '#define {}'.format(name)
        else:
            return '#define {} {}'.format(name, cls.cl_value(value))


    def __init__(self, name:str, value: str or bool or int or float):
        '''
        Initializes a kernel option.

        Parameters
        ----------
        name: str
            Option name
        value: str, bool, int, float
            Option value.
        '''
        self.cl_options = [(name, value)]

    def _get_name(self) -> str:
        return self.cl_options[0][0]
    name = property(_get_name, None, None, 'Option name.')

    def _get_value(self) -> str or bool or int or float:
        return self.cl_options[0][1]
    value = property(_get_value, None, None, 'Option value.')

    def __repr__(self):
        return "{}('{}', {})".format(self.__class__.__name__, 
                                   self.cl_options[0][0], self.cl_options[0][1])

    def __str__(self):
        return self.__repr__()

class McBoolOption(McOption):
    '''
    Base Class of Monte Carlo kernel options that have a boolean type.
    '''
    def __init__(self, name: str, value: bool):
        '''
        Initializes a boolean kernel option.

        Parameters
        ----------
        name: str
            Property name
        value: bool
            Boolean value of the kernel option.
        '''
        super().__init__(name, bool(value))

    def get_define(self):
        return '#define {}'.format(self.name, self.value)

class McIntOption(McOption):
    '''
    Base Class of Monte Carlo kernel options that have an integer type.
    '''
    def __init__(self, name: str, value: int):
        '''
        Initializes an integer kernel option.

        Parameters
        ----------
        name: str
            Property name
        value: int
            Integer value of the kernel option.
        '''
        super().__init__(name, int(value))

class McFloatOption(McOption):
    def __init__(self, name: str, value: float):
        '''
        Initializes a floating-point kernel option.

        Parameters
        ----------
        name: str
            Property name
        value: float
            Floating-point value of the kernel option.
        '''
        super().__init__(name, float(value))

class McTypeOption(McOption):
    '''
    Base Class of Monte Carlo kernel options that have a label type
    (no quotes in the define).
    '''
    def __init__(self, name: str, value: str):
        '''
        Initializes a label (no quotes in the define) kernel option.

        Parameters
        ----------
        name: str
            Property name
        value: str
            String value of the kernel option.
        '''
        super().__init__(name, str(value))


class McMethod(McIntOption):
    '''
    Selects the Monte Carlo simulation method:
      - Albedo Weight (0)
      - Albedo Rejection (1)
      - Microscopic Beer Lambert (2)

    A detailed description of the available methods can be found in:

    A. Sassaroli and F. Martelli
    Equivalence of four Monte Carlo methods for photon migration in turbid media
    J. Opt. Soc. Am. A / Vol. 29, No. 10 / October 2012

    Note
    ----
    Default method is Albedo Weight. The Albedo Rejection is fast but
    produces noisy results. The Microscopic Beer-Lambert is useful for
    fluence simulations with voxel size that is smaller than the mean
    free path. 
    '''

    albedo_weight = McIntOption('MC_METHOD', 0)
    aw = albedo_weight

    albedo_rejection = McIntOption('MC_METHOD', 1)
    ar = albedo_rejection

    microscopic_beer_lambert = McIntOption('MC_METHOD', 2)
    mbl = microscopic_beer_lambert

    default = albedo_weight

    def __init__(self, value: int=0):
        '''
        Initializes the Monte Carlo simulation method option.

        Parameters
        ----------
        value: int
            Allowed values are:
              - 0 for Albedo Weight
              - 1 for Albedo Rejection or,
              - 2 for Microscopic Beer-Lambert.
        '''
        if value not in (0, 1, 2):
            raise ValueError('Allowed values are 0, 1 or 2!')
        super().__init__('MC_METHOD', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McUseFluenceCache(McBoolOption):
    '''
    Turn on or off the use of OpenCL fluence cache.
    Default is off.

    Note
    ----
    Fluence cache will improve performance for large deposition voxels, where
    there is a high likelihood that consecutive weight depositions will go
    into the same accumulator. 
    '''
    on = McBoolOption('MC_USE_FLUENCE_CACHE', True)

    off = McBoolOption('MC_USE_FLUENCE_CACHE', False)

    default = off

    def __init__(self, value: bool=False):
        '''
        Initializes the fluence cache option.

        Parameters
        ----------
        value: bool
            Use True to enable fluence cache.
        '''
        super().__init__('MC_USE_FLUENCE_CACHE', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McUseNativeMath(McBoolOption):
    '''
    Turn on or off the use of OpenCL native math.
    Default is off.

    Note
    ----
    Native math usually gives some performance benefit, but might not
    be fully compliant with precisions defined by the IEEE standards. 
    '''
    on = McBoolOption('MC_USE_NATIVE_MATH', True)

    off = McBoolOption('MC_USE_NATIVE_MATH', False)

    default = off

    def __init__(self, value: bool=False):
        '''
        Initializes the native math kernel option.

        Parameters
        ----------
        value: bool
            Use True to enable native math or False to disable native math.
        '''
        super().__init__('MC_USE_NATIVE_MATH', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McIntLutMemory(McTypeOption):
    '''
    OpenCL memory type used to hold the floating-point lookup table data.
    Use one of :py:attr:`__constant` or :py:attr:`__global`.
    Default is :py:attr:`__global`.

    Note
    ----
    Selecting "__global" memory will likely lead to a significant performance
    degradation, in particular on older GPUs.
    The amount of available "__constant" memory on GPUs is typically limited
    to about 64k.
    '''
    #class global(McObject):
    #    cl_options = [('MC_INT_LUT_ARRAY_MEMORY', '__global')]

    __constant = McTypeOption('MC_INT_LUT_ARRAY_MEMORY', 'constant')

    __global = McTypeOption('MC_INT_LUT_ARRAY_MEMORY', 'global')

    default = __constant

    def __init__(self, value: str='global'):
        '''
        Initializes the integer type lookup table memory type kernel option.

        Parameters
        ----------
        value: str
            Use "global" to move integer type lookup table data to the
            __global OpenCL memory or use "constant" to move the lookup table
            data to the __constant OpenCL memory.
        '''
        value = {'global': '__global', '__global':'__global',
                 'constant':'__constant', '__constant':'__constant'}.get(value)
        if value is None:
            raise ValueError(
                'Lookup table memory type must be one of '
                '"constant" or "global", but got "{}"!'.format(value)
            )
        super().__init__('MC_INT_LUT_ARRAY_MEMORY', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McFloatLutMemory(McTypeOption):
    '''
    OpenCL memory type used to hold the floating-point lookup table data.
    Use one of "constant_mem" or "global_mem".
    Default is "global_mem".

    Note
    ----
    Selecting "global_mem" memory will likely lead to a significant performance
    degradation, in particular on older GPUs.
    The amount of available "constant_mem" memory on GPUs is typically limited
    to about 64k.
    '''
    #class global(McObject):
    #    cl_options = [('MC_FP_LUT_ARRAY_MEMORY', '__global')]

    constant_mem = McTypeOption('MC_FP_LUT_ARRAY_MEMORY', '__constant')

    global_mem = McTypeOption('MC_FP_LUT_ARRAY_MEMORY', '__global')

    default = constant_mem

    def __init__(self, value: str='global'):
        '''
        Initializes the floating-point lookup table memory type kernel option.

        Parameters
        ----------
        value: str
            Use "global_mem" to move integer lookup table data to the __global
            OpenCL memory or use "constant_mem" to move the lookup table data
            to the __constant OpenCL memory.
        '''
        value = {'global': '__global', '__global':'__global',
                 'constant':'__constant', '__constant':'__constant'}.get(value)
        if value is None:
            raise ValueError(
                'Lookup table memory type must be one of '
                '"constant" or "global", but got "{}"!'.format(value)
            )
        super().__init__('MC_FP_LUT_ARRAY_MEMORY', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McDebugMode(McBoolOption):
    '''
    Switch the kernel debug mode on or off. Use only for kernel development.
    Default is off.
    '''
    on = McBoolOption('MC_ENABLE_DEBUG', True)

    off = McBoolOption('MC_ENABLE_DEBUG', False)

    default = off

    def __init__(self, value: bool):
        '''
        Initializes the debug mode kernel option.

        Parameters
        ----------
        value: bool
            Use True to enable the debug mode or False to disable the debug mode.
        '''
        super().__init__('MC_ENABLE_DEBUG', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])

class McUseSoft64Atomics(McBoolOption):
    '''
    Force software implementation of 64-bit atomic operations.
    '''
    on = McBoolOption('MC_USE_SOFT_64_ATOMICS', True)

    off = McBoolOption('MC_USE_SOFT_64_ATOMICS', False)

    default = off

    def __init__(self, value: bool):
        '''
        Initializes the software-based 64-bit integer atomics usage option.

        Parameters
        ----------
        value: bool
            Use True to force software-implemented 64-bit integer atomics.
        '''
        super().__init__('MC_USE_SOFT_64_ATOMICS', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])

class McUseLottery(McBoolOption):
    '''
    Switch on or off photon packet termination by lottery.
    Default is on.
    '''
    on = McBoolOption('MC_USE_LOTTERY', True)

    off = McBoolOption('MC_USE_LOTTERY', False)

    default = on

    def __init__(self, value: bool):
        '''
        Initializes the lottery kernel option.

        Parameters
        ----------
        value: bool
            Use True to enable the lottery or False to disable the lottery.
        '''
        super().__init__('MC_USE_LOTTERY', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McMinimumPacketWeight(McFloatOption):
    '''
    Minimum packet weight allowed before starting packet termination or lottery.
    Default value is 1e-4.
    '''
    default = McFloatOption('MC_PACKET_WEIGHT_MIN', 1e-4)

    def __init__(self, value: float=1e-4):
        if 0.0 >= value or value > 1.0:
            raise ValueError(
                'Minimum photon packet weight must be > 0.0 and <= 1.0, '
                'but got {}!'.format(value)
            )
        super().__init__('MC_PACKET_WEIGHT_MIN', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McPacketLotteryChance(McFloatOption):
    '''
    Terminate photon packet by lottery if the value of a uniform random
    number from [0, 1] exceeds this value.
    Default value is 0.1.
    '''
    default = McFloatOption('MC_PACKET_LOTTERY_CHANCE', 0.1)

    def __init__(self, value: float=0.1):
        '''
        Initializes the lottery chance kernel option.

        Parameters
        ----------
        value: float
            Lottery survival fraction from 0.0 to 1.0.
        '''
        if  0.0 > value or value > 1.0:
            raise ValueError(
                'Lottery chance must be a value from 0.0 to 1.0, '
                'but got {}!'.format(value)
            )
        super().__init__('MC_PACKET_LOTTERY_CHANCE', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])


class McUsePackedStructures(McBoolOption):
    '''
    Turn on/off the use of tightly packed structures.
    Default value is off.

    Note
    ----
    Note that packed structures can lead to significant performance
    degradation of the MonteCarlo kernel. This option is the last resort if
    the fields of the OpenCL and host structures cannot be properly aligned.
    When declaring OPenCL or host structures always start with the
    largest data type and move towards smaller data types. Use data types
    that are of size no less than 4 bytes.
    '''
    on = McBoolOption('MC_USE_PACKED_STRUCTURES', True)

    off = McBoolOption('MC_USE_PACKED_STRUCTURES', False)

    default = off

    def __init__(self, value:bool):
        '''
        Initializes the packed structures kernel option.

        Parameters
        ----------
        value: bool
            Use True to enable packed structures or
            False to disable packed structures.
        '''
        super().__init__('MC_USE_PACKED_STRUCTURES', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])

Options = List[McOption]
