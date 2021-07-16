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

from xopto.mcml import mc
from xopto.cl import clinfo

import numpy as np
from matplotlib import pyplot as pp


# creating the sample layer and the 2 layers of the surrounding medium
layers = mc.mclayer.Layers(
    [
        mc.mclayer.Layer(n=1.00, d=np.inf,  mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(n=1.33, d=10.0e-3, mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.8)),
        mc.mclayer.Layer(n=1.00, d=np.inf,  mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.0))
    ]
)


# creating the source of photon packets
source = mc.mcsource.Line()


# creating surface detectors of backscattered or transmitted light
detectors = mc.mcdetector.Detectors(
    top = mc.mcdetector.Radial(
        mc.mcdetector.Axis(0.0, 10.0e-3, 1000),
        cosmin=np.cos(np.deg2rad(20.0))
    ),
    bottom = mc.mcdetector.Radial(
        mc.mcdetector.Axis(0.0, 10.0e-3, 100)
    )
)


# selecting the first available OpenCL GPU device
gpu = clinfo.gpu()


# creating a Monte Carlo simulator
mc_obj = mc.Mc(layers, source, detectors, cl_devices=gpu)
mc_obj.rmax = 25.0e-3


# running the MC simulation with 1,000,000 photon packets
trace_res, fluence_res, detectors_res = mc_obj.run(10e6)


# plotting the simulation results
fig, (ax1, ax2) = pp.subplots(2, 1)

ax1.semilogy(detectors_res.top.r*1e3, detectors_res.top.reflectance)
ax1.set_xlabel('Distance from source (mm)')
ax1.set_ylabel('Reflectance')

ax2.semilogy(detectors_res.bottom.r*1e3, detectors_res.bottom.reflectance)
ax2.set_xlabel('Distance from source (mm)')
ax2.set_ylabel('Transmittance')

pp.tight_layout()
pp.show()
