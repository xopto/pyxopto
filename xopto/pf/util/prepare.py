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

import os

import numpy as np

from xopto.pf.util import GkMap, MHgMap, MPcMap, \
                           MiePolystyreneMap, MieFractalPolystyreneMap, \
                           GkPolygon, MHgPolygon, MPcPolygon

from xopto import DATA_PATH


def precalculate(verbose: bool = True):
    '''
    Precalculate lookup tables for Gk, MHg, MPc, MiePolystyrene and
    MieFractalPolystyreneMap scattering phase function and
    precalculate the domain boundary points for the Gk, MHg and MPc
    scattering phase functions.

    Parameters
    ----------
    verbose: bool
        Turn on verbose progress report.
    '''
    for cls in [GkMap, MHgMap, MPcMap, MiePolystyreneMap]:
        target_dir = os.path.join(cls.DEFAULT_MAP_PATH, cls.DEFAULT_MAP_FILE)
        target_dir = os.dirname(os.path.abspath(target_dir))
        try:
            os.makedirs(target_dir)
        except OSError:
            pass
        cls.precalculate(verbose=verbose)

    for cls in [GkPolygon, MHgPolygon, MPcPolygon]:
        target_dir = os.path.join(
            cls.DEFAULT_POLYGON_PATH, cls.DEFAULT_POLYGON_FILE)
        target_dir = os.dirname(os.path.abspath(target_dir))
        try:
            os.makedirs(target_dir)
        except OSError:
            pass
        cls.precalculate(verbose=verbose)

def plot():
    '''
    Plot the Gk, MHg and MPc scattering phase function properties in terms
    of g, gamma and delta parameters.
    '''
    gkmap = GkMap.fromfile()
    mhgmap = MHgMap.fromfile()
    mpcmap = MPcMap.fromfile()
    miepolymap = MiePolystyreneMap.fromfile()
    # fracmiepolymap = FractalMiePolystyreneMap.fromfile()

    # plot maps
    from matplotlib import pyplot as pp
    gkmap.showMaps()
    mhgmap.showMaps()
    mpcmap.showMaps()
    miepolymap.showMaps()

    pp.figure()
    ax = pp.subplot(111)
    ax.plot(gkmap.gamma().flatten(), gkmap.delta().flatten(),
            'r+', label=gkmap.pfTypeName())
    ax.plot(mhgmap.gamma().flatten(), mhgmap.delta().flatten(),
            'b+', label=mhgmap.pfTypeName())
    ax.plot(mpcmap.gamma().flatten(), mpcmap.delta().flatten(),
            'g+', label=mpcmap.pfTypeName())
    ax.plot(miepolymap.gamma().flatten(), miepolymap.delta().flatten(),
            'k+', label=miepolymap.pfTypeName())
    ax.legend(loc='upper left')

    pp.show()

if __name__ == '__main__':
    precalculate()
    plot()
