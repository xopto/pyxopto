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

from xopto.mcbase.mcoptions import *

class McMaterialMemory(McTypeOption):
    '''
    OpenCL memory type used to hold the array of materials.
    Use one of "constant_mem" or "global_mem".
    Default is "global_mem".

    Note
    ----
    Selecting "constant_mem" memory will likely lead to a significant
    performance boost, in particular on older GPUs.
    However, note that the amount of available "constant_mem" memory on GPUs is
    typically limited to about 64k.
    '''
    #class global(McObject):
    #    cl_options = [('MC_MATERIAL_ARRAY_MEMORY', '__global')]

    constant_mem = McTypeOption('MC_MATERIAL_ARRAY_MEMORY', '__constant')

    global_mem = McTypeOption('MC_MATERIAL_ARRAY_MEMORY', '__global')

    default = constant_mem

    def __init__(self, value: str='global'):
        '''
        Initializes kernel option that sets the OpenCL memory type used for
        the array of materials.

        Parameters
        ----------
        value: str
            Use "global_mem" to move the array of materials to the __global
            OpenCL memory or use "constant_mem" to move the array of materials
            to the __constant OpenCL memory.
        '''
        value = {'global': '__global', '__global':'__global',
                 'constant':'__constant', '__constant':'__constant'}.get(value)
        if value is None:
            raise ValueError(
                'Material data memory type must be one of '
                '"constant" or "global", but got "{}"!'.format(value)
            )
        super().__init__('MC_MATERIAL_ARRAY_MEMORY', value)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.cl_options[0][1])
