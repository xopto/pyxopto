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

# This example covers all the neccessary steps for simulating 
# sampling volume utilizing multimode optical fibers as 
# sources and detectors. Sampling volume gives us an 
# understanding of which part of the turbid medium under 
# investigation is primarily responsible for the given 
# detected signal (see Meglinsky, I. V., and S. J. Matcher, 
# Med. Biol. Eng. Comput., 39(1), 44-50 (2001).)

from xopto.mcml import mc
from xopto.mcml.mcutil import fiber
from xopto.cl import clinfo

import numpy as np
import matplotlib.pyplot as pp

# DEFINE A COMPUTATIONAL DEVICE
# select the first device from the provided platform
cl_device = clinfo.gpu(platform='nvidia')

# DEFINE OPTICAL PROPERTIES FOR A SINGLE-LAYER MEDIUM
d = 1e-2
layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
    mc.mclayer.Layer(d=d, n=1.33, mua=2e2, mus=250e2, pf=mc.mcpf.Hg(0.90)),
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
])

# DEFINE UNIFORM MULTIMODE FIBER SOURCE
source = mc.mcsource.UniformFiber(
    fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
    ),
    direction=(0.0, 0.0, 1.0),
    position=(-220e-6, 0.0, 0.0)
)

# DEFINE MULTIMODE FIBER DETECTOR (IN THIS EXAMPLE ONLY TOP IS USED)
detector_top = mc.mcdetector.FiberArray(
    [fiber.FiberLayout(
        fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
        ),
        position=(220e-6, 0.0, 0.0),
        direction=(0.0, 0.0, 1.0)
    ),]
)

detectors = mc.mcdetector.Detectors(
    top=detector_top
)

# DEFINE THE TRACE OJECT WITH FILTER
nphotons = 100000  # should be in the order of 10000 for recording traces
recorded_traces = 1000  # desired number of recorded traces / only useful when filtering traces
number_of_events = 500  # maximal number of photon packet interaction events

trace = mc.mctrace.Trace(
    maxlen=number_of_events, options=mc.mctrace.Trace.TRACE_ALL,
    filter=mc.mctrace.Filter(
        z=(-1e-9, 1e-9), 
        pz=(-1, -np.cos(np.arcsin(0.22/1.45))), 
        r=(0, 100e-6, (220e-6, 0))
    )
)

# DEFINE THE SAMPLING VOLUME OBJECT
nx = 200
ny = 200
nz = 200
sv = mc.mcsv.SamplingVolume(
    xaxis=mc.mcsv.Axis(-0.75e-3, 0.75e-3, nx), 
    yaxis=mc.mcsv.Axis(-0.75e-3, 0.75e-3, ny),
    zaxis=mc.mcsv.Axis(0e-3, 1e-3, nz)
)

# DEFINE MC OBJECT FOR MONTE CARLO SIMULATIONS
mc_obj = mc.Mc(
    layers=layers, source=source, 
    detectors=detectors, trace=trace, 
    cl_devices=cl_device
)
mc_obj.rmax = 1e-1

# RUN MONTE CARLO SIMULATIONS AND COMPUTE THE SAMPLING VOLUME
output = None
while output is None or len(output[0]) < recorded_traces:
    output = mc_obj.run(nphotons, verbose=True, out=output)
    print('Photon packet traces collected:', len(output[0]))
trace = output[0]
sv = mc_obj.sampling_volume(trace, sv)

pp.figure()
pp.title('Sampling volume')
pp.imshow(sv.data.sum(axis=1), extent=(-0.75, 0.75, 0, 1))
pp.xlabel('x (mm)')
pp.ylabel('z (mm)')
pp.show()
