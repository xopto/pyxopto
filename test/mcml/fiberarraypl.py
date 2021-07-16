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

from xopto.mcml import mc
from xopto.mcml.mcutil.fiber import MultimodeFiber, FiberLayout
from xopto import make_user_dirs, USER_TMP_PATH

exportsrc = os.path.join(USER_TMP_PATH, 'mcml.h')
make_user_dirs()

pf = mc.mcpf.Hg(0.8)
cl_device = mc.clinfo.gpus()[0]
layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf),
    mc.mclayer.Layer(d=np.inf, n=1.33, mua=1e2, mus=100e2, pf=pf),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=1e2, mus=0.0, pf=pf),
])

fiber = MultimodeFiber(200e-6, 220e-6, 1.45, 0.22)

source = mc.mcsource.UniformFiber(fiber)

fibers = [
    FiberLayout(fiber, (0.5e-3, 0.0, 0.0)),
    FiberLayout(fiber, (1.0e-3, 0.0, 0.0))
]

pl_detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.FiberArrayPl(
        fibers, plaxis=mc.mcdetector.Axis(0, 10e-3, 250)
    ),
)

detectors = mc.mcdetector.Detectors(
    top=mc.mcdetector.FiberArray(fibers),
)

pl_mc_obj = mc.Mc(layers, source, pl_detectors, cl_devices=cl_device,
           options=[mc.mcoptions.McDebugMode.off])
pl_mc_obj.rmax = 25e-3

mc_obj = mc.Mc(layers, source, detectors, cl_devices=cl_device,
           options=[mc.mcoptions.McDebugMode.off])
mc_obj.rmax = 25e-3

res = mc_obj.run(4e6, exportsrc=exportsrc, verbose=True)[-1]
pl_res = pl_mc_obj.run(4e6, exportsrc=exportsrc, verbose=True)[-1]

pp.figure()
pp.suptitle('Two 200 um fibers at 250 and 500 um SDS.')
pp.subplot(121)
pp.semilogy(res.top.reflectance, 'og', label='FiberArray')
pp.semilogy(pl_res.top.reflectance.sum(0), 'xr', label='FiberArrayPl')
pp.xlabel('Fiber index')
pp.legend()

t = pl_res.top.pl/3e8
pp.subplot(122)
pp.semilogy(t*1e9, pl_res.top.reflectance[:, 0], label='250 um')
pp.semilogy(t*1e9, pl_res.top.reflectance[:, 1], label='500 um')
pp.xlabel('Time (ns)')
pp.legend()

pp.show()
