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

# This example covers the photon packet tracing in a two-layered 
# turbid medium utilizing multimode optical fibers as sources 
# and detectors

from xopto.mcml import mc
from xopto.mcml.mcutil import fiber
from xopto.cl import clinfo

import numpy as np
import matplotlib.pyplot as pp

# DEFINE A COMPUTATIONAL DEVICE
# select the first device from the provided platform
cl_device = clinfo.gpu(platform='nvidia')

# DEFINE OPTICAL PROPERTIES FOR EACH LAYER
d1 = 1e-3
d2 = 1e-3
layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),  # layer above the medium
    mc.mclayer.Layer(d=d1, n=1.33, mua=1e2, mus=125e2, pf=mc.mcpf.Hg(0.8)),
    mc.mclayer.Layer(d=d2, n=1.4, mua=2e2, mus=125e2, pf=mc.mcpf.Hg(0.9)),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),  # layer below the medium
])

# DEFINE UNIFORM MULTIMODE FIBER SOURCE
source = mc.mcsource.UniformFiber(
    fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
    ),
    position=(-440e-6, 0, 0),
    direction=(0.0, 0.0, 1.0)
)

# DEFINE MULTIMODE FIBER DETECTOR (IN THIS EXAMPLE ONLY BOTTOM IS USED)
detector_bottom = mc.mcdetector.FiberArray(
    [fiber.FiberLayout(
        fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
        ),
        position=(440e-6, 0, d1 + d2),
        direction=(0.0, 0.0, -1.0)
    ),]
)

detectors = mc.mcdetector.Detectors(
    bottom=detector_bottom
)

# DEFINE TRACE OBJECT
nphotons = 10000  # should be in the order of 10000 for recording traces
recorded_traces = 100  # desired number of recorded traces / only useful when filtering traces
number_of_events = 500  # maximal number of photon packet interaction events

# without filter
trace_no_filter = mc.mctrace.Trace(
    maxlen=number_of_events, options=mc.mctrace.Trace.TRACE_ALL,
)

# using default simulator data types
eps = mc.mctypes.McDataTypes.eps

# with filter
filter = mc.mctrace.Filter(
    z=(d1 + d2 - eps, d1 + d2 + eps), 
    pz=(np.cos(np.arcsin(0.22/1.45)), 1.00),
    r=(0, 100e-6, (440e-6, 0))
)
trace_filter = mc.mctrace.Trace(
    maxlen=number_of_events, options=mc.mctrace.Trace.TRACE_ALL, 
    filter=filter
)

# DEFINE MC OBJECT FOR MONTE CARLO SIMULATIONS
# without filter
mc_obj_no_filter = mc.Mc(
    layers=layers, 
    source=source,
    detectors=detectors,
    trace=trace_no_filter, 
    cl_devices=cl_device
)
mc_obj_no_filter.rmax = 1e-2

# with filter
mc_obj_filter = mc.Mc(
    layers=layers, 
    source=source,
    detectors=detectors,
    trace=trace_filter, 
    cl_devices=cl_device
)
mc_obj_filter.rmax = 1e-2

# RUN MONTE CARLO SIMULATIONS
# without filter
output = mc_obj_no_filter.run(nphotons, verbose=True)
trace_no_filter = output[0]

# with filter
output = mc_obj_filter.run(nphotons, verbose=True)
while len(output[0]) < recorded_traces:
    output = mc_obj_filter.run(nphotons, verbose=True, out=output)
    print('Photon packet traces detected:', len(output[0]))
trace_filter = output[0]

# VISUALIZE THE TRACES ONLY FIRST 100 TRACES
# (X-Z PLANE)
# without filter
pp.figure()
for i in range(recorded_traces):
    pp.plot(
        1e3*trace_no_filter.data['x'][i,:trace_no_filter.n[i]],
        1e3*trace_no_filter.data['z'][i,:trace_no_filter.n[i]]
    )
pp.title('Traces of all photon packets')
pp.xlabel('x (mm)')
pp.ylabel('z (mm)')
pp.xlim((-0.75, 0.75))
pp.ylim((0, 1e3*(d1 + d2)))
pp.gca().invert_yaxis()

# with filter
pp.figure()
inds = np.arange(number_of_events)
for i in range(recorded_traces):
    pp.plot(
        1e3*trace_filter.data['x'][i,:trace_filter.n[i]],
        1e3*trace_filter.data['z'][i,:trace_filter.n[i]]
    )
pp.title('Traces of filtered photon packets')
pp.xlabel('x (mm)')
pp.ylabel('z (mm)')
pp.xlim((-0.75, 0.75))
pp.ylim((0, 1e3*(d1 + d2)))
pp.gca().invert_yaxis()

pp.show()