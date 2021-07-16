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

from xopto.mcbase import cltypes
from xopto.mcbase.mcobject import McObject

'''
Since the scattering phase function parameters are passed to the OpenCL
kernel, packing of the structure elements becomes important. If all the
members of the McPf struct are of size 4-bytes, no packing is required, since
the parent McLayer structure comprises only 4-byte single precision floats.
'''

class PfBase(McObject):

    @classmethod
    def fromdict(cls, data: dict) -> 'PfBase':
        '''
        Create a new instance of a scattering phase function from a dict.
        The dict keys must match the parameter names defined by the
        constructor.
        '''
        T = data.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Cannot initialize an instance of {:s} scattering phase '
                'function from the provided data!'.format(cls.__name__))
        return cls(**data)

    def __repr__(self):
        return self.__str__() + ' # id 0x{:>08X}.'.format(id(self))
