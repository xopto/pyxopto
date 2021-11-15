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

# This example shows how to simulate reflectance as acquired with optical fibers 
# and optical fiber probes. 

from xopto.mcml import mc
from xopto.mcml.mcutil import fiber
from xopto.util import convolve
from xopto.cl import clinfo

import numpy as np

# DEFINE A COMPUTATIONAL DEVICE
# select the first device from the provided platform
cl_device = clinfo.gpu(platform='nvidia')

# DEFINE OPTICAL PROPERTIES FOR A SINGLE-LAYER MEDIUM
d = 1e-2
layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
    mc.mclayer.Layer(d=d, n=1.33, mua=0.0, mus=100e2, pf=mc.mcpf.Hg(0.8)),
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
])

# DEFINE UNIFORM MULTIMODE FIBER SOURCE
source = mc.mcsource.UniformFiber(
    fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
    ),
    direction=(0.0, 0.0, 1.0),
    position=(0.0, 0.0, 0.0)
)

# DEFINE SEVERAL DETECTORS (IN THIS EXAMPLE ONLY TOP IS USED)
detector_top_fiber = mc.mcdetector.FiberArray(
    [fiber.FiberLayout(
        fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
        ),
        position=(220e-6, 0.0, 0.0),
        direction=(0.0, 0.0, 1.0)
    ),]
)

detector_top_radial = mc.mcdetector.Radial(
    raxis=mc.mcdetector.RadialAxis(
        start=0.0, stop=1e-3, n=250
    ),
    cosmin=np.cos(np.arcsin(0.22/1.45))
)

detector_top_probe = mc.mcdetector.SixAroundOne(
    fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
    ),
    position=(0.0, 0.0),
    direction=(0.0, 0.0, 1.0)
)

# DEFINE SURFACE TO REALISTIC OPTICAL PROBE CASE
surface_realistic = mc.mcsurface.SixAroundOne(
    fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
    ),
    position=(0.0, 0.0),
    direction=(0.0, 0.0, 1.0),
    spacing=220e-6,
    diameter=6e-3,
    reflectivity=0.55,
    cutout=330e-6,
    cutoutn=1.6
)

surface = mc.mcsurface.SurfaceLayouts(top=surface_realistic)

# DEFINE MC OBJECT FOR MONTE CARLO SIMULATIONS
mc_obj_fiber = mc.Mc(
    layers=layers, 
    source=source, 
    detectors=mc.mcdetector.Detectors(top=detector_top_fiber), 
    cl_devices=cl_device
)

mc_obj_radial = mc.Mc(
    layers=layers, 
    source=source, 
    detectors=mc.mcdetector.Detectors(top=detector_top_radial), 
    cl_devices=cl_device
)

mc_obj_probe = mc.Mc(
    layers=layers, 
    source=source, 
    detectors=mc.mcdetector.Detectors(top=detector_top_probe), 
    cl_devices=cl_device
)

mc_obj_probe_realistic = mc.Mc(
    layers=layers, 
    source=source,
    surface=surface,
    detectors=mc.mcdetector.Detectors(top=detector_top_probe), 
    cl_devices=cl_device
)

# set the termination radius for all Monte Carlo simulator objects
mc_obj_fiber.rmax = 1e-2
mc_obj_radial.rmax = 1e-2
mc_obj_probe.rmax = 1e-2
mc_obj_probe_realistic.rmax = 1e-2

# RUN MONTE CARLO SIMULATIONS AND COMPUTE THE SAMPLING VOLUME
nphotons = 1e8
detector_fiber = mc_obj_fiber.run(nphotons, verbose=True)[-1]
detector_radial = mc_obj_radial.run(nphotons, verbose=True)[-1]
detector_probe = mc_obj_probe.run(nphotons, verbose=True)[-1]
detector_probe_realistic = mc_obj_probe_realistic.run(
    nphotons, verbose=True)[-1]

# each returned detector object 
reflectance_fiber = detector_fiber.top.reflectance[0]  
reflectance_fiber2 = convolve.fiber_reflectance( 
    detector_radial.top.r,
    detector_radial.top.reflectance,
    sds=220e-6,
    dcore=200e-6
)[0,0]
reflectance_probe = detector_probe.top.reflectance[1:].mean()
reflectance_probe_realistic = detector_probe_realistic.top.reflectance[1:].mean()

# COMPARE OBTAINED REFLECTANCE
print('{:.4e} - Reflectance for uniform multimode fiber detector.'
    .format(reflectance_fiber))
print('{:.4e} - Reflectance for uniform mutlimode fiber detector obtained'
    ' via annular ring integration'. format(reflectance_fiber2))
print('{:.4e} - Reflectance for optical probe.'
    .format(reflectance_probe))
print('{:.4e} - Reflectance for optical probe with realistic boundary.'
    .format(reflectance_probe_realistic))

# 2.9472e-04 - Reflectance detected with uniform multimode fiber detector.
# 2.9483e-04 - Reflectance detected with uniform mutlimode fiber detector obtained via annular ring integration
# 2.9375e-04 - Reflectance detected with six-around-one optical probe with uniform medium-probe interface.
# 3.1453e-04 - Reflectance detected with optical probe with realistic boundary.
