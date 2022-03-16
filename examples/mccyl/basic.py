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

from xopto.mccyl import mc
from xopto.cl import clinfo

import numpy as np
from matplotlib import pyplot as pp


# creating the sample layer and the layer of the surrounding medium
layers = mc.mclayer.Layers(
    [
        mc.mclayer.Layer(n=1.00, d=np.inf,  mua=0.0e2, mus=0.0e2, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(n=1.00, d=100e-3,  mua=0.0e2, mus=0.0e2, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(n=1.45, d=10.0e-3, mua=0.0e2, mus=0.0e2, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(n=1.33, d=8.0e-3,  mua=1.0e2, mus=50.0e2, pf=mc.mcpf.Hg(0.8))
    ]
)


# creating the source of photon packets
source = mc.mcsource.GaussianBeam(100e-6)


# creating polar surface detector
detectors = mc.mcdetector.Detectors(
    outer=mc.mcdetector.FiZ(
        fiaxis=mc.mcdetector.Axis(-np.pi, np.pi, 360),
        zaxis=mc.mcdetector.Axis(-1.5e-3, 1.5e-3, 3)
    ),
)


# selecting the first available OpenCL GPU device
gpu = clinfo.gpu()


# creating a Monte Carlo simulator
mc_obj = mc.Mc(layers, source, detectors, cl_devices=gpu)
mc_obj.rmax = 1000.0e-3


# running the MC simulation with 10,000,000 photon packets
trace_res, fluence_res, detectors_res = mc_obj.run(10e6, verbose=True)


# plotting the simulation results
fig, ax = pp.subplots()

ax.semilogy(
    np.rad2deg(detectors_res.outer.fi),
    detectors_res.outer.reflectance[1]
)
ax.set_xlabel('Angle (°)')
ax.set_ylabel('Reflectance/transmittance')

pp.tight_layout()
pp.show()
