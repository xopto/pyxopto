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

import matplotlib.pyplot as pp
import numpy as np

from xopto.util import convolve
from xopto.mcml import mc
from xopto import make_user_dirs, USER_TMP_PATH


debug = False

exportsrc = os.path.join(USER_TMP_PATH, 'mcml.h')
make_user_dirs()

pf = mc.mcpf.Hg(0.8)
cl_device = mc.clinfo.gpus()[0]
layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=np.inf, n=1.33, mua=1e2, mus=100e2, pf=pf),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=1e2, mus=0.0, pf=pf),
])

source = mc.mcsource.Line()

pl_detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.TotalPl(mc.mcdetector.Axis(0, 10e-3, 500)),
)

detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.Total(),
)

pl_mc_obj = mc.Mc(layers, source, pl_detectors, cl_devices=cl_device,
           options=[mc.mcoptions.McDebugMode(debug)])
pl_mc_obj.rmax = 25e-3

mc_obj = mc.Mc(layers, source, detectors, cl_devices=cl_device,
           options=[mc.mcoptions.McDebugMode(debug)])
mc_obj.rmax = 25e-3

nphotons = 4e6 if not debug else 1
res = mc_obj.run(nphotons, exportsrc=exportsrc, verbose=True)[-1]
pl_res = pl_mc_obj.run(nphotons, exportsrc=exportsrc, verbose=True)[-1]

pp.figure()

pp.subplot(121)
pp.semilogy(res.top.reflectance, 'og', label='Total')
t = pl_res.top.pl/3e8
pp.semilogy(pl_res.top.reflectance.sum(0), 'xr', label='TotalPl')
pp.legend()

pp.subplot(122)
pp.semilogy(t*1e9, pl_res.top.reflectance, label='TotalPl')
pp.title('Total reflectance')
pp.xlabel('Time (ns)')

pp.show()
