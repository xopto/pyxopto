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
        mc.mclayer.Layer(n=1.00, d=np.inf, mua=0.0e2, mus=0.0e2, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(n=1.00, d=50.0e-3, mua=0.0e2, mus=0.0e2, pf=mc.mcpf.Hg(0.0)),
        mc.mclayer.Layer(n=1.45, d=10.0e-3, mua=1.0e2, mus=0.0e2, pf=mc.mcpf.Hg(0.0)),
    ]
)


# creating the source of photon packets
source = mc.mcsource.GaussianBeam(1e-3)


# creating the trace object that can hold 128 events per each photon packet
trace = mc.mctrace.Trace(128)


# selecting the first available OpenCL GPU device
gpu = clinfo.gpu()


# creating a Monte Carlo simulator
mc_obj = mc.Mc(layers, source, None, trace, cl_devices=gpu)
mc_obj.rmax = 100.0e-3


# running the MC simulation with 100 photon packets
trace_res, fluence_res, detectors_res = mc_obj.run(100, verbose=True)


# crate a plot window and a plot axis
fig, ax = pp.subplots()

# plot the packet traces projected onto the x-y plane
trace_res.plot(view='xy', ax=ax)
# plot the glass cylindr
cylinder = pp.Circle((0.0, 0.0), 5.0e-3, color='blue', alpha=0.25,
                     fill=True, label='glass cylinder', zorder=2.0)
ax.add_patch(cylinder)
# plot the outer boundary of the sample
detector = pp.Circle((0.0, 0.0), 25e-3, color='black',
                     fill=False, label='sample boundary')
ax.add_patch(detector)

# add axis labels
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

# make aspect ratio of the x and y axis equal
ax.set_aspect('equal')

# show plot legend
pp.legend()

# display the plot window
pp.show()
