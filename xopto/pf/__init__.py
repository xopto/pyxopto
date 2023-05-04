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

"""Scattering phase functions."""

from .pfbase import PfBase, lut_function, ilut_function, fastg, fastgs

from .rayleigh import Rayleigh
from .hg import Hg
from .mhg import MHg
from .gk import Gk
from .gk2 import Gk2
from .mgk import MGk
from .doublehg import DoubleHg
from .pc import Pc
from .mpc import MPc
from .discrete import Discrete

from .mie import Mie
from .miepd import MiePd
from .mienormal import MieNormal
from .miefractal import MieFractal
from .miepolystyrene import MiePolystyrene, MieFractalPolystyrene, MieNormalPolystyrene

from .miemixture import MieMixture

try:
    from .mieml import MieMl
    from .miemlpd import MieMlPd
    from .miemlnormal import MieMlPd
    from .miemlfractal import MieMlFractal
except ImportError:
    import xopto
    if xopto.VERBOSE > 1:
        print(
            'Scattering phase functions for multilayered spherical particles '\
            'not available! Missing scattnlay module!'
        )

from . import util
