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

from xopto.mcml import mcobject

class Source(mcobject.McObject):
    @classmethod
    def fromdict(cls, data: dict) -> 'Source':
        '''
        Create a new instance of a photon packet source from a dict.
        The dict keys must match the parameter names defined by the
        constructor.
        '''
        data_ = dict(data)
        T = data_.pop('type')
        if T != cls.__name__:
            raise TypeError(
                'Cannot initialize an instance of {:s} photon packet source '
                'from the data of "{}"!'.format(cls.__name__, T))
        return cls(**data_)

    def __repr__(self):
        return '{} # id {}'.format(self.__str__(), id(self))
