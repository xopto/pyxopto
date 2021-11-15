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

# This example shows reflectance simulations from a suspension of 
# 1 μm polystyrene microspheres as detected via a multimode optical
# fiber position tightly next to the source fiber. The example
# gives a comparison how the selection of the scattering phase function
# (in this case Mie and Henyey-Greenstein) can affect the simulated 
# reflectance even though the scattering coefficient and anisotropy factor 
# are the same

from xopto.mcml import mc
from xopto.mcml.mcutil import fiber
from xopto.util import convolve
from xopto.pf import miepolystyrene
from xopto.materials.ri import polystyrene, water
from xopto.cl import clinfo

import numpy as np
import matplotlib.pyplot as pp

# CALCULATE THE SPECTRUM OF THE SCATTERING COEFFICIENT
# initialize the refractive index functions
ri_poly = polystyrene.Nikolov()
ri_water = water.Daimon()

# determine the target reduced scattering coefficient
# and calculate the multiplikative factor
wavelength = 600e-9
musr_target = 15e2
mie = miepolystyrene.MiePolystyrene(
    diameter=1.0e-6,
    wavelength=wavelength,
    ripolystyrene=ri_poly(wavelength),
    riliquid=ri_water(wavelength)
)

mus_target = musr_target/(1-mie.g(1))
n = mus_target/mie.scs()

# determine the scattering coefficient spectrum
wavelengths = np.arange(450e-9, 801e-9, 2e-9)
mus = np.zeros_like(wavelengths)
g = np.zeros_like(wavelengths)
for i, w in enumerate(wavelengths):

    mie = miepolystyrene.MiePolystyrene(
        diameter=1.0e-6,
        wavelength=w,
        ripolystyrene=ri_poly(w),
        riliquid=ri_water(w)
    )

    mus[i] =  n * mie.scs()
    g[i] = mie.g(1)

# plot the reduced scattering coefficient spectrum
pp.figure()
pp.title('Spectrum of the reduced scattering coefficient')
pp.plot(1e9 * wavelengths, 1e-2 * mus * (1-g))
pp.ylabel('Reduced scattering coefficient (1/cm)')
pp.xlabel('Wavelength (nm)')
pp.show()

# DEFINE A COMPUTATIONAL DEVICE
# select the first device from the provided platform
cl_device = clinfo.gpu(platform='nvidia')

# DEFINE OPTICAL PROPERTIES FOR A SINGLE-LAYER MEDIUM 
# (DYNAMICALLY CHANGED AT SIMULATION RUNS)
d = 1e-2
layers_hg = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
    mc.mclayer.Layer(d=d, n=1.33, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.8)),
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mc.mcpf.Hg(0.0)),
])

mie = miepolystyrene.MiePolystyrene(
    diameter=1.0e-6,
    wavelength=w,
    ripolystyrene=ri_poly(w),
    riliquid=ri_water(w)
)

mie_pf = mc.mcpf.Lut(*mie.mclut())
layers_mie = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mie_pf),
    mc.mclayer.Layer(d=d, n=1.33, mua=0.0, mus=0.0, pf=mie_pf),
    mc.mclayer.Layer(d=0.0, n=1.45, mua=0.0, mus=0.0, pf=mie_pf),
])

# DEFINE UNIFORM MULTIMODE FIBER SOURCE
source = mc.mcsource.UniformFiber(
    fiber=fiber.MultimodeFiber(
        dcore=200e-6, dcladding=220e-6, ncore=1.45, na=0.22
    ),
    direction=(0.0, 0.0, 1.0),
    position=(0.0, 0.0, 0.0)
)

# DEFINE DETECTOR (IN THIS EXAMPLE TOP IS USED)
detector_top_radial = mc.mcdetector.Radial(
    raxis=mc.mcdetector.RadialAxis(
        start=0.0, stop=1e-3, n=250
    ),
    cosmin=np.cos(np.arcsin(0.22/1.45))
)

# DEFINE MC OBJECT FOR MONTE CARLO SIMULATIONS
mc_obj_hg = mc.Mc(
    layers=layers_hg, 
    source=source, 
    detectors=mc.mcdetector.Detectors(top=detector_top_radial), 
    cl_devices=cl_device
)

mc_obj_mie = mc.Mc(
    layers=layers_mie, 
    source=source, 
    detectors=mc.mcdetector.Detectors(top=detector_top_radial), 
    cl_devices=cl_device
)

mc_obj_hg.rmax = 1e-2
mc_obj_mie.rmax = 1e-2

# DEFINE REFLECTANCE ARRAYS, RUN SIMULATIONS
# AND PROCESS THE REFLECTANCE
reflectance_hg = np.zeros_like(wavelengths)
reflectance_mie = np.zeros_like(wavelengths)

nphotons = 1e8
for i, w in enumerate(wavelengths):

    mc_obj_hg.layers[1].mus = mus[i]
    mc_obj_hg.layers[1].pf = mc.mcpf.Hg(g[i])

    mc_obj_mie.layers[1].mus = mus[i]

    mie = miepolystyrene.MiePolystyrene(
        diameter=1.0e-6,
        wavelength=w,
        ripolystyrene=ri_poly(w),
        riliquid=ri_water(w)
    )
    mie_pf = mc.mcpf.Lut(*mie.mclut())
    mc_obj_mie.layers[1].pf = mie_pf

    detector_hg = mc_obj_hg.run(nphotons, verbose=True)[-1]
    detector_mie = mc_obj_mie.run(nphotons, verbose=True)[-1]

    reflectance_hg[i] = convolve.fiber_reflectance( 
        detector_hg.top.r,
        detector_hg.top.reflectance,
        sds=220e-6,
        dcore=200e-6
    )[0,0]

    reflectance_mie[i] = convolve.fiber_reflectance( 
        detector_mie.top.r,
        detector_mie.top.reflectance,
        sds=220e-6,
        dcore=200e-6
    )[0,0]

# PLOT THE REFLECTANCE SPECTRUM
pp.figure()
pp.plot(1e9*wavelengths, 100*reflectance_hg, label='Henyey-Greenstein pf')
pp.plot(1e9*wavelengths, 100*reflectance_mie, label='Mie pf')
pp.legend()
pp.xlabel('Wavelength (nm)')
pp.ylabel('Reflectance (%)')
pp.show()
