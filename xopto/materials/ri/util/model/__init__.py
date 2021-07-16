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

from .base import Model, Scale, Normalize
from .cauchy import Cauchy
from .conrady import Conrady_1, ConradyEx_1, \
                     Conrady_2, ConradyEx_2
from .exponential import Exponential
from .herzberg import Herzberger_3_2, HerzbergerEx_3_2, \
                      Herzberger_4_2, HerzbergerEx_4_2
from .sellmeier import Sellmeier_1, SellmeierEx_1, \
                       Sellmeier_2, SellmeierEx_2, \
                       Sellmeier_3, SellmeierEx_3, \
                       Sellmeier_4, SellmeierEx_4, \
                       Sellmeier_5, SellmeierEx_5

Schott_2_4 = Cauchy
