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

import os.path
import tempfile

import matplotlib as mpl
import matplotlib.pyplot as pp
from matplotlib.animation import FuncAnimation
import numpy as np

from xopto.mcml import mc
from xopto import make_user_dirs, USER_TMP_PATH


exportsrc = os.path.join(USER_TMP_PATH, 'mcml_fluencerz.h')
make_user_dirs()

cl_device = mc.clinfo.gpus()[0]

pf = mc.mcpf.Hg(0.8)

layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=1.0, n=1.0, mua=1.0e2, mus=50.0e2, pf=pf),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
])

source = mc.mcsource.Line()

detectors = mc.mcdetector.Detectors(
    bottom=mc.mcdetector.Total()
)
fluence = mc.mcfluence.FluenceRz(
    raxis=mc.mcfluence.Axis(0.0e-3, 1.0e-3, 101),
    zaxis=mc.mcfluence.Axis( 0.0e-3, 1.0e-3, 101),
)

nphotons = 50e6

mc_obj = mc.Mc(layers, source, detectors, fluence=fluence,
               cl_devices=cl_device)
mc_obj.rmax = 2.0e-3

_, fluence_res, detectors_res = mc_obj.run(
    nphotons, exportsrc=exportsrc, verbose=True)

fluence_res.plot()

