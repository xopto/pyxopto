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

# This is a sample script for Monte Carlo calculations and
# visualization for reflectance as acquired with a 1 cm
# diameter opening integrating sphere from a 
# sample of human skin

from xopto.mcml import mc
from xopto.cl import clinfo
from xopto.materials.absorption import oxyhem, deoxyhem

import numpy as np
import matplotlib.pyplot as pp

# DEFINE A COMPUTATIONAL DEVICE
# select the first device from the provided platform
cl_device = clinfo.gpu(platform='nvidia')

# DEFINE RELEVANT SIMULATION PARAMETERS
nphotons = 1e6
wavelengths = np.arange(450e-9, 801e-9, 2e-9)

# DEFINE OPTICAL PROPERTIES FOR TWO-LAYERED HUMAN SKIN
# layer 1 - EPIDERMIS
d1 = 100e-6  # layer thickness in m
n1 = 1.4  # refractive index
m = 0.02  # melanin volume fraction
g1 = 0.8  # anisotropy factor constant with wavelength
pf1 = mc.mcpf.Hg(g1)  # Henyey-Greenstein scatterin phase function

# epidermis absortpion coefficient
mua1 = lambda wavelength: m * 6.6*1e13*(1e9*wavelength)**-3.33 + \
    (1-m) * 1e2*0.5*(0.244 + 85.3*np.exp(-(1e9*wavelength - 154)/66.2))

# epidermis scattering coefficient
mus1 = lambda wavelength: (2*1e7*(1e9*wavelength)**-1.5 + \
    2*1e14*(1e9*wavelength)**-4) / (1-g1)

# layer 2 - DERMIS
d2 = 10e-3  # layer thickness in m
n2 = 1.4  # refractive index
bl = 0.02  # blood volume fraction
oxy = 0.90  # oxygenation
g2 = 0.8  # anisotropy factor
pf2 = mc.mcpf.Hg(g2)  # Henyey-Greenstein scatterin phase function

# dermis absorption coefficient
mua_oxy = oxyhem.OxyHem()
mua_deoxy = deoxyhem.DeOxyHem()

mua2 = lambda wavelength: bl * (oxy * mua_oxy(wavelength, None) + \
    (1-oxy) * mua_deoxy(wavelength, None)) + \
    (1-bl) * 1e2 * (0.244 + 16.82*np.exp(-(1e9*wavelength - 400) / 80.5))

# dermis scattering coefficient
mus2 = lambda wavelength: (2*1e7*(1e9*wavelength)**-1.5 + \
    2*1e14*(1e9*wavelength)**-4) / (1-g2)

# DEFINE TWO-LAYERED SKIN MODEL LAYER STACK
layers = mc.mclayer.Layers([
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf1),  # layer above the medium
    mc.mclayer.Layer(d=d1, n=n1, mua=1.0e2, mus=1.0e2, pf=pf1),
    mc.mclayer.Layer(d=d2, n=n2, mua=1.0e2, mus=1.0e2, pf=pf2),
    mc.mclayer.Layer(d=0.0, n=1.0, mua=0.0, mus=0.0, pf=pf1),  # layer below the medium
])

# DEFINE SOURCE
source = mc.mcsource.Line(
    position=(0.0, 0.0, 0.0),
    direction=(0.0, 0.0, 1.0)
)

# DEFINE A DETECTOR FOR INTEGRATING SPHERE
sp_r = 0.5e-2  # integrating sphere opening in m
detector_top = mc.mcdetector.Radial(
    mc.mcdetector.RadialAxis(
        start=0.0,
        stop=2*sp_r,
        n=2)
)
detectors = mc.mcdetector.Detectors(
    top=detector_top
)

# DEFINE THE MC OBJECT FOR SIMULATIONS
mc_obj = mc.Mc(
    layers=layers, 
    source=source,
    detectors=detectors, 
    cl_devices=cl_device
)
mc_obj.rmax = 10e-2

# DO SIMULATIONS FOR DESIRED WAVELENGTH RANGE
reflectance_spectrum = np.zeros_like(wavelengths)
for i, w in enumerate(wavelengths):

    # for each wavelength redefine optical properties
    mc_obj.layers[1].mua = mua1(w)
    mc_obj.layers[1].mus = mus1(w)
    mc_obj.layers[2].mua = mua2(w)
    mc_obj.layers[2].mus = mus2(w)

    detector = mc_obj.run(nphotons, verbose=True, wgsize=256)[-1]
    reflectance_spectrum[i] = detector.top.reflectance[0] * np.pi * sp_r**2

# PLOT THE REFLECTANCE SPECTRUM
pp.figure()
pp.title('Reflectance spectrum of human skin')
pp.plot(1e9*wavelengths, 100*reflectance_spectrum)
pp.xlabel('Wavelength (nm)')
pp.ylabel('Reflectance (%)')
pp.show()